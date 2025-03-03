//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <ctime> 
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

////////// Other functions //////////

// Return maximum of 3 numbers
// [[Rcpp::export]]
double ret_max(double a,double b,double c){
  vec temp={a,b,c};
  return(max(temp));
}

// [[Rcpp::export]]
double ret_min(double a,double b,double c){
  vec temp={a,b,c};
  return(min(temp));
}

// [[Rcpp::export]]
int vec_to_num(vec X){
  if(X.max()>1 || X.min()<0){
    Rcpp::Rcout <<"Error with entries in the vector"<< std::endl;
    return -1;
  }
  int n=X.n_rows;
  int number=0;
  for(int i=0;i<n;i++){
    if(X(i)==1){
      number+=std::pow(2,i); 
    }
  }
  return number;
}

// [[Rcpp::export]]
vec num_to_vec(int n, int d){
  vec X(d);
  X.zeros();
  int temp;
  if(n<0 || n>=std::pow(2,d)){
    Rcpp::Rcout <<"Error, number bigger than dimension.\n Returning vector of 0s " << std::endl;
    return(X);
  }else{
    for(int i=0;i<d;i++){
      // Rcpp::Rcout <<"iteration: " <<i<< std::endl;
      temp=n % 2;
      X(i)=temp;
      if(temp==1){n = (n-1)/2;}
      if(temp==0){n=n/2;}
    }
  }
  return(X);
}


////////// Balancing functions //////////

///// List of balancing functions that apply to log-likelihoods
// [[Rcpp::export]]
double bf_sq(double x){
  return x/2;
}

// [[Rcpp::export]]
double bf_min(double x){
  double result;
  if(x<0){
    result = x;
  }else{
    result = 0;
  }
  return result;
}

//// Bounded balancing function based on sqrt 
// [[Rcpp::export]]
double bound_sq(double x, double log_gamma){
  double temp1=x/2-log_gamma;
  return ret_min(temp1,x,0);
}


///// Functions to call balancing functions inside other functions

//This function we don't export to R because it generates errors
double invoke(double x, double (*func)(double)) {
  return func(x);
}

// [[Rcpp::export]]
double bal_func(double x,String chosen){
  if (chosen == "sq") {
    return invoke(x, &bf_sq);
  } else if (chosen == "min") {
    return invoke(x, &bf_min);
  } else {
    cout << "Unknown operation!" << endl;
    Rcpp::Rcout <<"Unknown operation!" << std::endl;
    return 0; // Default return for unknown operation
  }
}


////////// loglikelihood functions //////////


// 7 modes log-likelihood, dimension 16
// [[Rcpp::export]]
double loglik(const arma::vec& X){
  int size=X.n_rows;
  if(size!=16){
   Rcpp::Rcout <<"Error:Size is not 16" << std::endl;
  return -std::numeric_limits<double>::infinity();
  }
  
  double theta=15;

  // Defined modes
    arma::vec  mod1 = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};// 16 1s
    arma::vec  mod2 = {1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0};// 8 1s
    arma::vec  mod3 = {0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1};// 8 1s
    arma::vec  mod4 = {1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0};// 8 1s
    arma::vec  mod5 = {0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1};// 8 1s
    arma::vec  mod6 = {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1};// 2 1s
    arma::vec  mod7 = {0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0};// 2 1s
  
  double loglik_comp=0;
  
  loglik_comp+=exp(-theta*sum(abs(X-mod1)));
  loglik_comp+=exp(-theta*sum(abs(X-mod2)));
  loglik_comp+=exp(-theta*sum(abs(X-mod3)));
  loglik_comp+=exp(-theta*sum(abs(X-mod4)));
  loglik_comp+=exp(-theta*sum(abs(X-mod5)));
  loglik_comp+=exp(-theta*sum(abs(X-mod6)));
  loglik_comp+=exp(-theta*sum(abs(X-mod7)));
  
  return log(loglik_comp);
}

////////// Updating functions //////////

///// Functions to update individual replicas
// [[Rcpp::export]]
List IIT_update_w(vec X, String chosen_bf, double temperature){
  int total_neighbors = X.n_rows; // total number of neighbors is p spacial
  vec probs(total_neighbors, fill::zeros); //probabilities
  
////// Compute likelihood of current state
  double logpi_current=0;
  // uvec current_coord = find(X==1);
  logpi_current = loglik(X);
  ////// Finish computing likelihood of current state
  
  // Rcpp::Rcout << "current loglik: "<< logpi_current << std::endl;
  ////// Compute weight for all neighbors
  double temporal=0;
  vec newX;
  for(int j=0; j<total_neighbors;j++){
    // Rcpp::Rcout << "Starts checking neighbors  "<< j<<std::endl; 
    newX = X;
    newX.row(j) = 1-X.row(j);
    //Rcpp::Rcout << newX << std::endl;
    // uvec coord = find(newX==1);
    //Rcpp::Rcout << coord << std::endl;
    temporal=loglik(newX)-logpi_current;
    //Apply balancing function to log probability times temperature ////
    probs(j)=bal_func(temporal*temperature, chosen_bf);
  }
  //////Choose the next neighbor
  vec u = Rcpp::runif(total_neighbors);
  vec probs_choose = log(-log(u)) - probs;
  
  //Find the index of the minimum element. 
  //This corresponds to choosing that neighbor
  int neigh_pos = (std::min_element(probs_choose.begin(), probs_choose.end()))-probs_choose.begin();
  // Rcpp::Rcout <<"probs vector: "<< probs << std::endl;
  // Rcpp::Rcout <<"chosen neighbor: "<< neigh_pos << std::endl;
  
  X.row(neigh_pos) = 1-X.row(neigh_pos); //modify the coordinate of the chosen neighbor
  List result;
  result["X"]=X; // Return new state
  result["Z"]=sum(exp(probs))/total_neighbors; //Compute Z factor
  return result;
}

// [[Rcpp::export]]
List a_IIT_update(vec X, String chosen_bf, double temperature, double log_bound, bool decrease, double prob_to_dec, double decreasing_constant){
  int total_neighbors = X.n_rows; // total number of neighbors is p spacial
  vec logprobs(total_neighbors, fill::zeros); //vector to store log-probabilities
  vec max_logprobs(total_neighbors,fill::zeros);//vector to store max-log-probabilities
  ////// Compute likelihood of current state
  double logpi_current=0;
  logpi_current = loglik(X);
  ////// Compute weight for all neighbors
  double temporal=0;
  // double temp_bound = 0;
  vec newX;
  for(int j=0; j<total_neighbors;j++){
    newX = X;
    newX.row(j) = 1-X.row(j); //Change coordinate of the state to define new neighbor
    temporal = temperature*(loglik(newX)-logpi_current);
    logprobs(j)=temporal; //Store raw log_probability
    max_logprobs(j)=abs(temporal); // Store the max log-probability, either pi_y/pi_x or pi_x/pi_y
  }// End of loop to compute raw log-probability of neighbors
  
  //Updating the log-bound
  log_bound=ret_max(bal_func(max(max_logprobs),"sq"),log_bound,0);
  
  //Arithmetic reduction of the bounding constant
  if(decrease){//If it's defined to reduce the bounding constant
    if(prob_to_dec>0){//If we consider a probability to decrease the constant
      double test_prob=0;
      if(prob_to_dec<1){
        vec ppp=Rcpp::runif(1);//Get a random number
        test_prob=ppp(0);
      }
      if(test_prob<prob_to_dec){//In case the update is accepted
        if(exp(log_bound)<1+decreasing_constant){//If the bound is less than 1+constant then we can only decrease it to 1 
          log_bound=0;
        }else{
          log_bound-=decreasing_constant; //Reduce the constant by decreasing_constant units
        }  
      }
    }
  }
  
  //After the bound has been updated we apply the bounded balancing function to the vector of log probabilities
/////////////
////IMPORTANT: we're only using bound sqrt root for the adaptive IIT
/////////////
   // Rcpp::Rcout <<"log-probs vector: \n"<< logprobs << std::endl;
  for(int j=0;j<total_neighbors;j++){
    logprobs(j)=bound_sq(logprobs(j),log_bound);
  }
  // Rcpp::Rcout <<"log-probs vector after BF: \n"<< logprobs << std::endl;
  //////Choose the next neighbor
  vec u = Rcpp::runif(total_neighbors);
  vec probs_choose = log(-log(u)) - logprobs;
  
  //Find the index of the minimum element. 
  //This corresponds to choosing that neighbor
  int neigh_pos = (std::min_element(probs_choose.begin(), probs_choose.end()))-probs_choose.begin();
  // Rcpp::Rcout <<"probs vector: "<< probs << std::endl;
  // Rcpp::Rcout <<"chosen neighbor: "<< neigh_pos << std::endl;
  
  X.row(neigh_pos) = 1-X.row(neigh_pos); //modify the coordinate of the chosen neighbor
  List ret;
  ret["X"]=X;
  ret["Z"]=sum(exp(logprobs))/total_neighbors; // Compute Z factor with uniform proposal distribution
  ret["logbound"]=log_bound;
  return ret;
}



// // [[Rcpp::export]]
// void IPT_update(mat& logprob_matrix, mat& states_matrix, vec& temp,vec& index_process, const std::vector<std::string>& bal_function){
//   int rows=logprob_matrix.n_rows;
//   int cols=logprob_matrix.n_cols;
//   int total_neighbors=states_matrix.n_rows;
//   double temporal=0;
//   double logpi_current;
// ////  Rcpp::Rcout << "Rows: " << rows <<", cols:"<<cols<<", tot_neigh: "<<total_neighbors<< std::endl;
//   std::string current_bal_fun;
//   double current_temp;
//   mat random_matrix=mat(rows,cols,fill::randu);
//   // Rcpp::Rcout << "Random matrix: \n" << random_matrix << std::endl;
//   random_matrix=log(-log(random_matrix));
//   //Rcpp::Rcout << "Transformed random matrix: \n" << random_matrix << std::endl;
//   
//   mat result_matrix = random_matrix - logprob_matrix;
//   // Rcpp::Rcout << "Result matrix: \n" << result_matrix << std::endl;
//   // uword index=result_matrix.index_min();
//   // Rcpp::Rcout << "Chosen index: " << index << std::endl;
//   uvec index=ind2sub(size(rows,cols),result_matrix.index_min());
//   uword neighbor=index(0);//Randomly chosen neighbor
//   uword replica=index(1);//Randomly chosen replica
// //If the last row is chosen, it's a replica swap
//   if(neighbor==rows-1){
//     //Update index process
//     int mod=index_process.max()+1;// Define modulus for temperature
//     vec epsilon_indic(cols);
//     epsilon_indic.fill(0);
//     uvec replica_to_swap_v;
//     uword replica_to_swap;
//     
//     int temp_to_swap=((static_cast<int>(index_process(replica))+1)%mod);
//     replica_to_swap_v=find(index_process==temp_to_swap);
//     replica_to_swap=replica_to_swap_v(0);
//     Rcpp::Rcout << "Swapping replica: " <<replica<<" with temperature: "<<index_process(replica)<<"\n with replica: "<<replica_to_swap<<" with temp: "<<temp_to_swap<< std::endl;
//     
//     epsilon_indic(replica)=temp_to_swap-index_process(replica);
//     epsilon_indic(replica_to_swap)=index_process(replica)-index_process(replica_to_swap);
//     
//     // Rcpp::Rcout << "epsilon indic: \n"<<epsilon_indic<< std::endl;
//     // Rcpp::Rcout << "index_process before: \n"<<index_process<< std::endl;
//     index_process+=epsilon_indic;
//     // Rcpp::Rcout << "resulting_swap: \n"<<index_process<< std::endl;
// //// Finish index process swap
// //Update weights of the involved neighbors.
//     vec X(total_neighbors);
//     vec newX(total_neighbors);
//     
//     X=states_matrix.col(replica); // Current state of the replica updating
//     // Rcpp::Rcout << "First X to update: \n"<<X<< std::endl;
//     current_temp=temp(index_process(replica));//Use the new temperature
//     current_bal_fun=bal_function[index_process(replica)];
//     logpi_current=loglik(X);
//     ////// Compute weight for all neighbors
//     temporal=0;
//     for(int j=0; j<total_neighbors;j++){
//       // Rcpp::Rcout << "Starts checking neighbors  "<< j<<std::endl; 
//       newX = X;
//       newX.row(j) = 1-X.row(j);
//       temporal=loglik(newX)-logpi_current;
//       //Apply balancing function to log probability times temperature ////
//       logprob_matrix(j,replica)=bal_func(temporal*current_temp, current_bal_fun)-log(total_neighbors);
//     }
//     // Rcpp::Rcout << "W Matrix after first update: \n"<<logprob_matrix<< std::endl;
//     X=states_matrix.col(replica_to_swap); // Current state of the replica updating
//     // Rcpp::Rcout << "second X to update: \n"<<X<< std::endl;
//     current_temp=temp(index_process(replica_to_swap));
//     current_bal_fun=bal_function[index_process(replica_to_swap)];
//     logpi_current=loglik(X);
//     ////// Compute weight for all neighbors
//     temporal=0;
//     for(int j=0; j<total_neighbors;j++){
//       // Rcpp::Rcout << "Starts checking neighbors  "<< j<<std::endl; 
//       newX = X;
//       newX.row(j) = 1-X.row(j);
//       temporal=loglik(newX)-logpi_current;
//       //Apply balancing function to log probability times temperature ////
//       logprob_matrix(j,replica_to_swap)=bal_func(temporal*current_temp, current_bal_fun)-log(total_neighbors);
//     }
//     // Rcpp::Rcout << "W Matrix after second update: \n"<<logprob_matrix<< std::endl;
//     
//     
// // To make things easy update all temperature weights
// uvec temporal_vec;
// // Rcpp::Rcout << "Temperature vector: \n"<<temp<< std::endl;
// for(int r=0;r<cols;r++){
//   temp_to_swap=((static_cast<int>(index_process(r))+1)%mod); //Identify the +1 in the index process
//   // Rcpp::Rcout <<"index_process(r)="<<index_process(r)<< " temp to swap:"<<temp_to_swap<< std::endl;
//   replica_to_swap_v=find(index_process==temp_to_swap);//Find the index of the replica corresponding to that entry in the index process
//   replica_to_swap=replica_to_swap_v(0);
//   // Rcpp::Rcout << "swapping temp: "<<temp(index_process(r))<<"("<<index_process(r)<<")"<<" replica: "<<r<<" with temp: "<<temp(temp_to_swap)<<" in replica: "<<replica_to_swap<< std::endl;
//   
//   // Rcpp::Rcout << "temp_current: "<< temp(index_process(r)) <<std::endl;
//   // Rcpp::Rcout << "temp_to swap: "<< temp(temp_to_swap) <<std::endl;
//   // Rcpp::Rcout << "Replica current: "<< r <<"Likelihood:"<<loglik(states_matrix.col(r))<<"Contains vector:\n"<<states_matrix.col(r)<<std::endl;
//   // Rcpp::Rcout << "Replica to swap: "<< replica_to_swap <<"Likelihood:"<<loglik(states_matrix.col(replica_to_swap))<<"Contains vector:\n"<<states_matrix.col(replica_to_swap)<<std::endl;
//   temporal= (temp(temp_to_swap)-temp(index_process(r)))*(loglik(states_matrix.col(r))-loglik(states_matrix.col(replica_to_swap)));
//   // Rcpp::Rcout << "Difference of likelihoods: "<<(loglik(states_matrix.col(r))-loglik(states_matrix.col(replica_to_swap)))<<std::endl;
//   // Rcpp::Rcout << "Difference of temps: "<<temp(temp_to_swap)-temp(index_process(r))<<std::endl;
//   // Rcpp::Rcout << "log(cols): "<<log(cols)<<std::endl;
//   logprob_matrix(rows-1,r)=bal_func(temporal,bal_function[index_process(r)])-log(cols);
// }
// 
//   }else{//If it's not a replica swap
//     Rcpp::Rcout << "Changing neighbor: " << neighbor <<" of replica:"<<replica<< std::endl;
//     //Swap the chosen coordinate of the chosen replica
//     states_matrix(neighbor,replica)=1-states_matrix(neighbor,replica);
// ////Then compute the new weights for all neighbors
//     vec X=states_matrix.col(replica); // Current state of the replica updating
//     current_temp=temp(index_process(replica));
//     current_bal_fun=bal_function[index_process(replica)];
//     logpi_current=loglik(X);
//     ////// Compute weight for all neighbors
//     // Rcpp::Rcout << "Flag 1" << std::endl;
//     vec newX(total_neighbors);
//     for(int j=0; j<total_neighbors;j++){
//       // Rcpp::Rcout << "Starts checking neighbors  "<< j<<std::endl; 
//       newX = X;
//       newX.row(j) = 1-X.row(j);
//       temporal=loglik(newX)-logpi_current;
//       //Apply balancing function to log probability times temperature ////
//       logprob_matrix(j,replica)=bal_func(temporal*current_temp, current_bal_fun)-log(total_neighbors);
//     }
//     // Rcpp::Rcout << "Flag 2" << std::endl;
// ////Compute the new weights for the possible replica swaps
//     int mod=index_process.max()+1;// Define modulus for temperature
//     int temp_to_swap=((static_cast<int>(index_process(replica))+1)%mod);
//     uvec replica_to_swap_v=find(index_process==temp_to_swap);
//     uword replica_right=replica_to_swap_v(0);
//     // Rcpp::Rcout << "Flag 3" << std::endl;
//     temp_to_swap=((static_cast<int>(index_process(replica))-1)%mod);
//     if(temp_to_swap<0){temp_to_swap+=mod;}//Avoid having nefative indexes
//     // Rcpp::Rcout << "Flag 3.1" << std::endl;
//     replica_to_swap_v=find(index_process==temp_to_swap);
//     // Rcpp::Rcout << "replica_to_swap_v" <<replica_to_swap_v<<"temp to swap: "<<temp_to_swap<< std::endl;
//     // Rcpp::Rcout << "Flag 3.2" << std::endl;
//     uword replica_left=replica_to_swap_v(0);
//     // Rcpp::Rcout << "Replica to swap:  "<< replica<<"\n replica_right: "<<replica_right<<"\n replica_left: "<<replica_left<<std::endl; 
//     // Rcpp::Rcout << "Flag 4" << std::endl;
//     //Update replica swap weight of current replica
//     temporal=(current_temp-temp(replica_right))*(loglik(states_matrix.col(replica_right)) - loglik(X));
//     logprob_matrix(rows-1,replica)=bal_func(temporal,current_bal_fun)-log(cols);
//     temporal=(temp(replica_left)-current_temp)*(loglik(X)-loglik(states_matrix.col(replica_left)));
//     logprob_matrix(rows-1,replica_left)=bal_func(temporal,current_bal_fun)-log(cols);
//     // Rcpp::Rcout << "Flag 5" << std::endl;
// }
//   // List ret;
//   // ret["weights"]=logprob_matrix;
//   // ret["states"]=states_matrix; 
//   // ret["ip"]=index_process;
//   // return ret;  
// }

////////// Code for Parallel Tempering simulations //////////

// [[Rcpp::export]]
List PT_IIT_sim(int p,int startsim,int endsim, int numiter,int iterswap,int burn_in, vec temp, const std::vector<std::string>& bal_function, bool bias_fix, int initial_state){
//// Initialize variables to use in the code
  int T=temp.n_rows; // Count number of temperatures
  double J=double(T)-1;//Number of temperatures minus 1, used in swap loops
  int total_sim = (endsim-startsim+1); //Count total number of simulations
  int total_swaps=trunc(numiter/iterswap);
  List output; // To store output of the update function
  double Z; // To store Z factor of update function
  int swap_count; //to keep track of even-odd swaps
  double current_temp; // to temporarily store the temperature
//// Initialize arrays to store information
  mat X(p,T); // To store the current state of the joint chain, as many rows as neighbors, as many columns as temperatures
  vec index_process(T);   //Initialize index process vector
  mat ind_pro_hist(total_swaps*total_sim+1,T); //To store evolution of index process
  int max_num=pow(2,p);
  vec pi_est(max_num); //Vector to store the estimated weight for each state
  mat full_pi_est(max_num,total_sim); //Matrix to store the estimated weight considering all simulations
  vec first_visit(max_num);//Vector to store the first visit to each state
  mat full_first_visit(max_num,total_sim); //Matrix to store first visits considering all simulations
  vec swap_total(J);
  vec swap_success(J);
  mat swap_rate(total_sim,J);
  ////Variables to update index process
  vec epsilon_indic(T); //Vector that indicates if the entry of the index process is proposed to change
  vec prop_swap(T); //vector to indicate a proposed swap
  vec do_swap(T); //vector to indicate actually perform a swap
  vec resulting_swap(T);//Vector to update the index process
  vec ppp; //To store the probability of replica swap
//// Variables to perform the replica swap
  double swap_prob;
  vec Xtemp_from(p);
  vec Xtemp_to(p);
  // Variables to count modes visited
  // mat modes_visited(numiter * total_sim,T);//Matrix to store the modes visited and temperature
  std::vector<double> time_taken(total_sim); // vector to store the seconds each process took
  
//// Start the loop for all simulations
  for(int s=0;s<total_sim;s++){
    for(int i=0;i<T;i++){ // Reset index process vector at the start of each simulation
      index_process.row(i)=i;
    }
    ind_pro_hist.row(0)=index_process.t(); // First entry of the index process
    swap_count=0; //Reset swap count
    // X.zeros();//Reset the starting point of all chains
    vec initialX=num_to_vec(initial_state,p);
    for(int c=0;c<T;c++){
      X.col(c)=initialX;
    }
    pi_est.zeros(); // Reset the estimated distribution
    first_visit.zeros(); //Reset the vector of first visits
    swap_total.zeros();
    swap_success.zeros();
    
    //// Start loop for burn_in period
    for(int i=0;i<burn_in;i++){
      if (i % 100 == 1) {Rcpp::Rcout << "PT-IIT - Simulation: " << s+startsim << " Burn_in period, iteration: " << i << std::endl;}
      for(int replica=0;replica<T;replica++){//For loop for replica update
        current_temp=temp(index_process(replica));
        output=IIT_update_w(X.col(replica),bal_function[index_process(replica)],current_temp);
        X.col(replica)=vec(output(0)); //Update current state of the chain
      }
      //End replica update in burn-in period
      
      //In the burn-in period we don't modify the index process
      //Just directly swap replicas
      //This way when the simulation starts after burn-in period starts we have the original starting index process
      //Start replica swap in burn-in period
      if ((i+1) % iterswap == 0){
        swap_count+=1;
        int starting=swap_count%2;
        for(int t=starting;t<J;t+=2){// For loop that runs over temperature indexes to swap
          
          epsilon_indic.elem(find(index_process==t)).ones();
          prop_swap.elem(find(index_process==t)).ones(); //we swap temperature t
          prop_swap.elem(find(index_process==t+1)).ones(); //With t+1
          //Compute swap probability
          Xtemp_from=X.col(t);
          Xtemp_to=X.col(t+1);
          
          
          //// Optional Computation of Z factors to correct bias
          double Z_fact_correc=1;//to temporarily store Z_factor correction
          if(bias_fix){
            double Z_temp11;
            double Z_temp12;
            double Z_temp21;
            double Z_temp22;
            
            output=IIT_update_w(Xtemp_from,bal_function[t],temp(t));
            Z_temp11=output(1);
            output=IIT_update_w(Xtemp_to,bal_function[t+1],temp(t+1));
            Z_temp22=output(1);
            output=IIT_update_w(Xtemp_from,bal_function[t+1],temp(t+1));
            Z_temp12=output(1);
            output=IIT_update_w(Xtemp_to,bal_function[t],temp(t));
            Z_temp21=output(1);
            
            Z_fact_correc=Z_temp12*Z_temp21/(Z_temp11*Z_temp22);
          }
          //// Computing swap probability
          swap_prob=(temp(t)-temp(t+1))*(loglik(Xtemp_to) - loglik(Xtemp_from)); 
          swap_prob=Z_fact_correc*exp(swap_prob);
          // Rcpp::Rcout <<"Swap prob "<< swap_prob << std::endl;
          ppp=Rcpp::runif(1);
          if(ppp(0)<swap_prob){//In case the swap is accepted
            //Swap vectors in the matrix X
            X.col(t+1)=Xtemp_from;
            X.col(t)=Xtemp_to;
          }
        }
      }
    }// Finish burn-in period
    swap_count=0; //Reset swap count
    
    std::clock_t start = std::clock(); // Start timer for simulation s
//// Start the loop for all iterations in simulation s
    for(int i=0;i<numiter;i++){
      // Rcpp::Rcout <<"Inside iteration loop"<< i << std::endl;
      if (i % 1000 == 1) {Rcpp::Rcout << "Simulation: " << s+startsim << " Iteration: " << i << std::endl;}
      // Rcpp::Rcout << "Simulation: " << s+startsim << " Iteration: " << i << std::endl;
      for(int replica=0;replica<T;replica++){//For loop for replicas
        current_temp=temp(index_process(replica));// Extract temperature of the replica
 // Rcpp::Rcout <<"Inside replica loop, with replica: "<< replica << std::endl;
        //Depending on the chosen method
//// Update each replica independently
            output=IIT_update_w(X.col(replica),bal_function[index_process(replica)],current_temp);
//// Store Z factor of replica with temperature 1
            if(current_temp==1){ // For the original temperature replica
// Rcpp::Rcout << "Storing weight in simulation: " << s+startsim << " Iteration: " << i << std::endl;
              Z = output(1); //Extract the Z-factor
// Rcpp::Rcout << "Printing Z: " << Z << std::endl;
              int state=vec_to_num(X.col(replica));
// Rcpp::Rcout << "Printing state: " << state << std::endl;
// pi_est(state)++;//Count how many times each state was visited
              pi_est(state)+=(1/Z);//Add weight
              if(first_visit(state)==0){
                first_visit(state)=i;//Store the first time the state is visited 
              }
// Rcpp::Rcout << "Printing pi_est: " << pi_est << std::endl;
// Rcpp::Rcout << "All good with Storing weight in simulation: " << s+startsim << " Iteration: " << i << std::endl;
            }
            X.col(replica)=vec(output(0)); //Update current state of the chain
      }//End loop to update replicas

//// Start replica swap process
      if ((i+1) % iterswap == 0){
        swap_count+=1;//Increase the count of swaps
        // Rcpp::Rcout << "Trying swap: " << swap_count << std::endl;
        epsilon_indic.fill(-1); //Epsilon indic starts as -1
        prop_swap.zeros();
        do_swap.zeros();
        //Try a replica swap every iterswap number of iterations
        //We're doing non-reversible parallel tempering
        int starting=swap_count%2; // Detect if it's even or odd
        // Rcpp::Rcout <<"Trying replica swap "<<swap_count<<" start: "<<starting <<" at iteration: "<< i << std::endl;
        for(int t=starting;t<J;t+=2){// For loop that runs over temperature indexes to swap
          swap_total(t)+=1;//Increase the count of tried swaps
          epsilon_indic.elem(find(index_process==t)).ones();
          prop_swap.elem(find(index_process==t)).ones(); //we swap temperature t
          prop_swap.elem(find(index_process==t+1)).ones(); //With t+1
          //Compute swap probability
          Xtemp_from=X.cols(find(index_process==t));
          Xtemp_to=X.cols(find(index_process==t+1));

          
  //// Optional Computation of Z factors to correct bias
  double Z_fact_correc=1;//to temporarily store Z_factor correction
  if(bias_fix){
    double Z_temp11;
    double Z_temp12;
    double Z_temp21;
    double Z_temp22;
    
    output=IIT_update_w(Xtemp_from,bal_function[t],temp(t));
    Z_temp11=output(1);
    output=IIT_update_w(Xtemp_to,bal_function[t+1],temp(t+1));
    Z_temp22=output(1);
    output=IIT_update_w(Xtemp_from,bal_function[t+1],temp(t+1));
    Z_temp12=output(1);
    output=IIT_update_w(Xtemp_to,bal_function[t],temp(t));
    Z_temp21=output(1);
    
    Z_fact_correc=Z_temp12*Z_temp21/(Z_temp11*Z_temp22);
    // Rcpp::Rcout <<"Z bias correction "<< Z_fact_correc << std::endl;
  }
//// Computing swap probability
          swap_prob=(temp(t)-temp(t+1))*(loglik(Xtemp_to) - loglik(Xtemp_from)); 
          swap_prob=Z_fact_correc*exp(swap_prob);
// Rcpp::Rcout <<"Swap prob "<< swap_prob << std::endl;
          ppp=Rcpp::runif(1);
          if(ppp(0)<swap_prob){//In case the swap is accepted
            swap_success(t)+=1;//Increase the number of successful swaps of temp t
// Rcpp::Rcout <<"Accepted swap " << std::endl;
            do_swap.elem(find(index_process==t)).ones();
            do_swap.elem(find(index_process==t+1)).ones();
          }
        }
        resulting_swap=epsilon_indic % prop_swap % do_swap;
        // Rcpp::Rcout <<"Resulting swap "<< resulting_swap << std::endl;
        index_process+=resulting_swap;
        // Rcpp::Rcout <<"New index process:\n "<< index_process << std::endl;
        ind_pro_hist.row((s*total_swaps)+swap_count)=index_process.t();
        // Rcpp::Rcout <<"Store index process " << std::endl;
      }//End of replica swap process
    }// End loop of iterations
    std::clock_t end = std::clock(); // Stop timer
    // Calculate the time taken in seconds
    double duration = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    time_taken[s] = duration;
// Store result of the simulation
full_pi_est.col(s)=pi_est;
full_first_visit.col(s)=first_visit;
vec temp_rate=swap_success / swap_total;
swap_rate.row(s)=temp_rate.t();
    // Rcpp::Rcout <<"Final state "<< X << std::endl;
  }//End loop simulations
  List ret;
  ret["est_pi"]=full_pi_est;
  ret["ip"]=ind_pro_hist;
  ret["visits"]=full_first_visit;
  ret["swap_rate"]=swap_rate;
  ret["time_taken"]=time_taken;
  return ret;
}


// [[Rcpp::export]]
List PT_a_IIT_sim(int p,int startsim,int endsim, int total_swaps,int sample_inter_swap,int burn_in, vec temp, const std::vector<std::string>& bal_function, int initial_state, bool bound_reduction, double decreasing_constant,std::string reduc_model){
  //// Initialize variables to use in the code
  int T=temp.n_rows; // Count number of temperatures
  vec log_bound_vector(T); // vector to store a log-bound for each replica
  double J=double(T)-1;//Number of temperatures minus 1, used in swap loops
  int total_sim = (endsim-startsim+1); //Count total number of simulations
  List output; // To store output of the update function
  double Z; // To store Z factor of update function
  int swap_count; //to keep track of even-odd swaps
  double current_temp; // to temporarily store the temperature
  double current_log_bound; //to temporarily store the log-bound
  int new_samples;//temporal variable to store weight with multiplicity list
  //// Initialize arrays to store information
  mat X(p,T); // To store the current state of the joint chain, as many rows as neighbors, as many columns as temperatures
  vec index_process(T);   //Initialize index process vector
  mat ind_pro_hist(total_swaps*total_sim+1,T); //To store evolution of index process
  int max_num=pow(2,p);
  // Rcpp::Rcout << "max num: " << max_num << std::endl;  
  vec pi_est(max_num); //Vector to store the estimated weight for each state
  mat full_pi_est(max_num,total_sim);
  vec first_visit(max_num);//Vector to store the first visit to each state
  mat full_first_visit(max_num,total_sim); //Matrix to store first visits considering all simulations
  vec swap_total(J,fill::ones);
  swap_total*=total_swaps;//We always have the same number of total swaps
  vec swap_success(J);
  mat swap_rate(total_sim,J);
  cube total_iterations(total_swaps,T,total_sim,fill::zeros);//To store iterations needed in between swaps
  ////Variables to update index process
  vec epsilon_indic(T); //Vector that indicates if the entry of the index process is proposed to change
  vec prop_swap(T); //vector to indicate a proposed swap
  vec do_swap(T); //vector to indicate actually perform a swap
  vec resulting_swap(T);//Vector to update the index process
  vec ppp; //To store the probability of replica swap
  //// Variables to perform the replica swap
  double swap_prob;
  vec Xtemp_from(p);
  vec Xtemp_to(p);
  // Variables to count modes visited
  // mat modes_visited(numiter * total_sim,T);//Matrix to store the modes visited and temperature
  std::vector<double> time_taken(total_sim); // vector to store the seconds each process took
  
  // Probability to update
  bool update_prob=false;
  double prob_to_dec=0;
  double percentage_start=0.05;
  double percentage_end=0.70;
  int total_replica_iterations=sample_inter_swap*total_swaps;
  int sample_iterations_count;
    if(reduc_model=="always"){prob_to_dec=1;}
    if(reduc_model=="never"){prob_to_dec=0;}
    if(reduc_model=="iterations"){update_prob=true;}
  
  //// Start the loop for all simulations
  for(int s=0;s<total_sim;s++){
    for(int i=0;i<T;i++){ // Reset index process vector at the start of each simulation
      index_process.row(i)=i;
    }
    ind_pro_hist.row(0)=index_process.t(); // First entry of the index process
    swap_count=0; //Reset swap count
    // X.zeros();//Reset the starting point of all chains
    vec initialX=num_to_vec(initial_state,p);
    for(int c=0;c<T;c++){
      X.col(c)=initialX;
    }
    pi_est.zeros(); // Reset the estimated distribution
    first_visit.zeros(); //Reset the vector of first visits
    log_bound_vector.zeros();//Reset log-bounds, all log-bounds start at 0
    swap_success.zeros();
    //Reset the probability to reduce the bounding constant
    if(reduc_model=="iterations"){update_prob=true;prob_to_dec=1;} //Reset the bool to update probability
    sample_iterations_count=0; // Reset the counting of iterations (or samples)
    ////Start the loop for burn-in period
    int track_burn_in=0;
    while(track_burn_in<burn_in){
      for(int replica=0;replica<T;replica++){//For loop for replica update in the burn-in
        int samples_replica=0;
        while(samples_replica<sample_inter_swap){//Loop to create samples for each replica until we reach the defined threshold
          current_temp=temp(replica);// Extract temperature of the replica
          current_log_bound=log_bound_vector(replica);// Extract log-bound of the corresponding temperature
          output=a_IIT_update(X.col(replica),bal_function[index_process(replica)],current_temp,current_log_bound,bound_reduction,prob_to_dec,decreasing_constant);
          //// Compute weight
          Z = output(1); //Extract the Z-factor
          new_samples=1+R::rgeom(Z);
          if(new_samples<1){
            Rcpp::Rcout <<"Error: geometric in "<< "simulation: " << s+startsim << " Burn-in period after " << track_burn_in <<"simulations,  temp:"<<current_temp<< std::endl;
            Rcpp::Rcout <<"new_samples= "<<new_samples<< ", Z=" << Z << " log-bound= " << current_log_bound << std::endl;
            new_samples=sample_inter_swap;
          }
          if(samples_replica+new_samples>sample_inter_swap){//If we're going to surpass the required number of samples
            new_samples = sample_inter_swap-samples_replica;//Wee force to stop at sample_inter_swap
          }
          samples_replica+=new_samples; // Update number of samples obtained from the replica
          X.col(replica)=vec(output(0)); //Update current state of the chain
          log_bound_vector(index_process(replica))=output(2); //Update log-bound 
        }
      }//End loop to update replicas in the burn-in
      
      //// Start replica swap process
      
      swap_count+=1;//Increase the count of swaps
      //Try a replica swap after reaching sample_inter_swap in each replica
      //We're doing non-reversible parallel tempering
      int starting=swap_count%2; // Detect if it's even or odd
      // Rcpp::Rcout <<"Trying replica swap "<<swap_count<<" start: "<<starting <<" at iteration: "<< i << std::endl;
      for(int t=starting;t<J;t+=2){// For loop that runs over temperature indexes to swap
        Xtemp_from=X.col(t);
        Xtemp_to=X.col(t+1);
        
        //// Computing swap probability
        swap_prob=(temp(t)-temp(t+1))*(loglik(Xtemp_to) - loglik(Xtemp_from)); 
        swap_prob=exp(swap_prob);
        // Rcpp::Rcout <<"Swap prob "<< swap_prob << std::endl;
        ppp=Rcpp::runif(1);
        if(ppp(0)<swap_prob){//In case the swap is accepted
          //Swap vectors in the matrix X
          X.col(t+1)=Xtemp_from;
          X.col(t)=Xtemp_to;
        }
      }
      track_burn_in+=sample_inter_swap;
      Rcpp::Rcout << "PT A-IITm - Simulation: " << s+startsim << ". Done " << track_burn_in <<" samples in burn-in period"<< std::endl;
    }
    ////Finish the loop for burn-in period
    swap_count=0; //Reset swap count
    
    std::clock_t start = std::clock(); // Start timer for simulation s
    //// Start the loop for all iterations in simulation s
    for(int i=0;i<total_swaps;i++){
      // Rcpp::Rcout <<"Inside iteration loop"<< i << std::endl;
      if (i % 10 == 1) {Rcpp::Rcout << "PT A-IITm - Simulation: " << s+startsim << " Swap: " << i <<" Prob_decrease_bound: " << prob_to_dec << std::endl;}
//<< " log-bound:\n " << log_bound_vector
    // Rcpp::Rcout << "Simulation: " << s+startsim << " Iteration: " << i << std::endl;
      for(int replica=0;replica<T;replica++){//For loop for replicas
        int samples_replica=0;
        while(samples_replica<sample_inter_swap){//Loop to create samples for each replica until we reach the defined threshold
          total_iterations(i,index_process(replica),s)+=1;//increase the number of iterations
          current_temp=temp(index_process(replica));// Extract temperature of the replica
          current_log_bound=log_bound_vector(index_process(replica));// Extract log-bound of the corresponding temperature
          
          if(update_prob){//Check if we need to update the probability
            if(current_temp==1){//The original replica defines the speed to modify the probabilty to decrease bounding constant
              sample_iterations_count+=samples_replica; //Add the number of iterations (or samples) from the previous step
              if(sample_iterations_count>total_replica_iterations*percentage_start){//Check if we start decreasing the probability
                if(sample_iterations_count>total_replica_iterations*percentage_end){//Check if we stop decreasing the probability
                  prob_to_dec=0;
                  update_prob=false;
                }else{//In case we haven't finished updating the probability
                  //Probability is proportional 
                  double progress=sample_iterations_count/total_replica_iterations;
                  prob_to_dec=1+((percentage_start-progress)/(percentage_end-percentage_start));
                }
              }
            }
          }
          
          
          output=a_IIT_update(X.col(replica),bal_function[index_process(replica)],current_temp,current_log_bound,bound_reduction,prob_to_dec,decreasing_constant);
          
          //// Compute weight
          Z = output(1); //Extract the Z-factor
          new_samples=1+R::rgeom(Z);
          if(new_samples<1){
            Rcpp::Rcout <<"Error with geometric in "<< "Simulation: " << s+startsim << " Swap: " << i <<" temperature:"<<current_temp<< std::endl;
          }
          if(samples_replica+new_samples>sample_inter_swap){//If we're going to surpass the required number of samples
            new_samples = sample_inter_swap-samples_replica;//Wee force to stop at sample_inter_swap
          }
          samples_replica+=new_samples; // Update number of samples obtained from the replica
//// Store weight of replica with temperature 1
          if(current_temp==1){ // For the original temperature replica
            int state=vec_to_num(X.col(replica));
            pi_est(state)+=new_samples;//Add weight
//// Check if it's the first time that replica with temperature 1 visits this state
            if(first_visit(state)==0){
              mat current_slice=total_iterations.slice(s);//Extract current slice
              //Store the first time the state is visited 
              first_visit(state)=sum(current_slice.col(index_process(replica)));
            }
          }
          X.col(replica)=vec(output(0)); //Update current state of the chain
          log_bound_vector(index_process(replica))=output(2); //Update log-bound 
        }
      }//End loop to update replicas
      
//// Start replica swap process

        swap_count+=1;//Increase the count of swaps
        // Rcpp::Rcout << "Trying swap: " << swap_count << std::endl;
        epsilon_indic.fill(-1); //Epsilon indic starts as -1
        prop_swap.zeros();
        do_swap.zeros();
        //Try a replica swap after reaching sample_inter_swap in each replica
        //We're doing non-reversible parallel tempering
        int starting=swap_count%2; // Detect if it's even or odd
        // Rcpp::Rcout <<"Trying replica swap "<<swap_count<<" start: "<<starting <<" at iteration: "<< i << std::endl;
        for(int t=starting;t<J;t+=2){// For loop that runs over temperature indexes to swap
          
          epsilon_indic.elem(find(index_process==t)).ones();
          prop_swap.elem(find(index_process==t)).ones(); //we swap temperature t
          prop_swap.elem(find(index_process==t+1)).ones(); //With t+1
          //Compute swap probability
          Xtemp_from=X.cols(find(index_process==t));
          Xtemp_to=X.cols(find(index_process==t+1));
          
          //// Computing swap probability
          swap_prob=(temp(t)-temp(t+1))*(loglik(Xtemp_to) - loglik(Xtemp_from)); 
          swap_prob=exp(swap_prob);
          // Rcpp::Rcout <<"Swap prob "<< swap_prob << std::endl;
          ppp=Rcpp::runif(1);
          if(ppp(0)<swap_prob){//In case the swap is accepted
            swap_success(t)+=1;//Increase the number of successful swaps of temp t
            do_swap.elem(find(index_process==t)).ones();
            do_swap.elem(find(index_process==t+1)).ones();
          }
        }
        resulting_swap=epsilon_indic % prop_swap % do_swap;
        // Rcpp::Rcout <<"Resulting swap "<< resulting_swap << std::endl;
        index_process+=resulting_swap;
        // Rcpp::Rcout <<"New index process:\n "<< index_process << std::endl;
        ind_pro_hist.row((s*total_swaps)+swap_count)=index_process.t();
        // Rcpp::Rcout <<"Store index process " << std::endl;
////End of replica swap process
    }// End loop of iterations
    std::clock_t end = std::clock(); // Stop timer
    // Calculate the time taken in seconds
    double duration = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    time_taken[s] = duration;
    // Store result of the simulation
    full_pi_est.col(s)=pi_est;
    full_first_visit.col(s)=first_visit;
    vec temp_rate=swap_success / swap_total;
    swap_rate.row(s)=temp_rate.t();
    // Rcpp::Rcout <<"Final state "<< X << std::endl;
  }//End loop simulations
  List ret;
  ret["est_pi"]=full_pi_est;
  ret["ip"]=ind_pro_hist;
  ret["visits"]=full_first_visit;
  ret["swap_rate"]=swap_rate;
  ret["total_iter"]=total_iterations;
  ret["time_taken"]=time_taken;
  return ret;
}

// [[Rcpp::export]]
List PT_a_IIT_sim_RF(int p,int startsim,int endsim, int numiter,int iterswap,int burn_in, vec temp, const std::vector<std::string>& bal_function, bool bias_fix, int initial_state, bool bound_reduction, double decreasing_constant,std::string reduc_model){
  //// Initialize variables to use in the code
  int T=temp.n_rows; // Count number of temperatures
  vec log_bound_vector(T); // vector to store a log-bound for each replica
  double J=double(T)-1;//Number of temperatures minus 1, used in swap loops
  int total_sim = (endsim-startsim+1); //Count total number of simulations
  int total_swaps=trunc(numiter/iterswap);
  List output; // To store output of the update function
  double Z; // To store Z factor of update function
  int swap_count; //to keep track of even-odd swaps
  double current_temp; // to temporarily store the temperature
  double current_log_bound; //to temporarily store the log-bound
  //// Initialize arrays to store information
  mat X(p,T); // To store the current state of the joint chain, as many rows as neighbors, as many columns as temperatures
  vec index_process(T);   //Initialize index process vector
  mat ind_pro_hist(total_swaps*total_sim+1,T); //To store evolution of index process
  int max_num=pow(2,p);
  vec pi_est(max_num); //Vector to store the estimated weight for each state
  mat full_pi_est(max_num,total_sim); //Matrix to store the estimated weight considering all simulations
  vec first_visit(max_num);//Vector to store the first visit to each state
  mat full_first_visit(max_num,total_sim); //Matrix to store first visits considering all simulations
  vec swap_total(J);
  vec swap_success(J);
  mat swap_rate(total_sim,J);
  ////Variables to update index process
  vec epsilon_indic(T); //Vector that indicates if the entry of the index process is proposed to change
  vec prop_swap(T); //vector to indicate a proposed swap
  vec do_swap(T); //vector to indicate actually perform a swap
  vec resulting_swap(T);//Vector to update the index process
  vec ppp; //To store the probability of replica swap
  //// Variables to perform the replica swap
  double swap_prob;
  vec Xtemp_from(p);
  vec Xtemp_to(p);
  // Variables to count modes visited
  // mat modes_visited(numiter * total_sim,T);//Matrix to store the modes visited and temperature
  std::vector<double> time_taken(total_sim); // vector to store the seconds each process took
  
  // Probability to update
  bool update_prob=false;
  double prob_to_dec=0;
  double percentage_start=0.05;
  double percentage_end=0.70;
  int total_replica_iterations=numiter;
  int sample_iterations_count;
  if(reduc_model=="always"){prob_to_dec=1;}
  if(reduc_model=="never"){prob_to_dec=0;}
  if(reduc_model=="iterations"){update_prob=true;}
  //// Start the loop for all simulations
  for(int s=0;s<total_sim;s++){
    for(int i=0;i<T;i++){ // Reset index process vector at the start of each simulation
      index_process.row(i)=i;
    }
    ind_pro_hist.row(0)=index_process.t(); // First entry of the index process
    swap_count=0; //Reset swap count
    // X.zeros();//Reset the starting point of all chains
    vec initialX=num_to_vec(initial_state,p);
    for(int c=0;c<T;c++){
      X.col(c)=initialX;
    }
    pi_est.zeros(); // Reset the estimated distribution
    first_visit.zeros(); //Reset the vector of first visits
    log_bound_vector.zeros();//Reset log-bounds, all log-bounds start at 0
    swap_total.zeros();
    swap_success.zeros();
    
    //Reset the probability to reduce the bounding constant
    if(reduc_model=="iterations"){update_prob=true;prob_to_dec=1;} //Reset the bool to update probability
    sample_iterations_count=0; // Reset the counting of iterations (or samples)
    //// Start loop for burn_in period
    for(int i=0;i<burn_in;i++){
      if (i % 100 == 1) {Rcpp::Rcout << "PT A-IITw - Simulation: " << s+startsim << " Burn_in period, iteration: " << i << std::endl;}
      for(int replica=0;replica<T;replica++){//For loop for replica update
        current_temp=temp(index_process(replica));
        current_log_bound=log_bound_vector(replica);// Extract log-bound of the corresponding temperature
        output=a_IIT_update(X.col(replica),bal_function[index_process(replica)],current_temp,current_log_bound,bound_reduction,prob_to_dec,decreasing_constant);
        X.col(replica)=vec(output(0)); //Update current state of the chain
        log_bound_vector(index_process(replica))=output(2); //Update log-bound 
      }
      //End replica update in burn-in period
      
      //In the burn-in period we don't modify the index process
      //Just directly swap replicas
      //This way when the simulation starts after burn-in period starts we have the original starting index process
      //Start replica swap in burn-in period
      if ((i+1) % iterswap == 0){
        swap_count+=1;
        int starting=swap_count%2;
        for(int t=starting;t<J;t+=2){// For loop that runs over temperature indexes to swap
          
          // epsilon_indic.elem(find(index_process==t)).ones();
          // prop_swap.elem(find(index_process==t)).ones(); //we swap temperature t
          // prop_swap.elem(find(index_process==t+1)).ones(); //With t+1
          //Compute swap probability
          Xtemp_from=X.col(t);
          Xtemp_to=X.col(t+1);
          
          
          //// Optional Computation of Z factors to correct bias
          double Z_fact_correc=1;//to temporarily store Z_factor correction
          if(bias_fix){
            double Z_temp11;
            double Z_temp12;
            double Z_temp21;
            double Z_temp22;
            
            output=a_IIT_update(Xtemp_from,bal_function[t],temp(t),log_bound_vector(t),false,prob_to_dec,decreasing_constant);
            Z_temp11=output(1);
              output=a_IIT_update(Xtemp_to,bal_function[t+1],temp(t+1),log_bound_vector(t+1),false,prob_to_dec,decreasing_constant);
            Z_temp22=output(1);
            output=a_IIT_update(Xtemp_from,bal_function[t+1],temp(t+1),log_bound_vector(t+1),false,prob_to_dec,decreasing_constant);
            Z_temp12=output(1);
            output=a_IIT_update(Xtemp_to,bal_function[t],temp(t),log_bound_vector(t),false,prob_to_dec,decreasing_constant);
            Z_temp21=output(1);
            
            Z_fact_correc=Z_temp12*Z_temp21/(Z_temp11*Z_temp22);
          }
          //// Computing swap probability
          swap_prob=(temp(t)-temp(t+1))*(loglik(Xtemp_to) - loglik(Xtemp_from)); 
          swap_prob=Z_fact_correc*exp(swap_prob);
          // Rcpp::Rcout <<"Swap prob "<< swap_prob << std::endl;
          ppp=Rcpp::runif(1);
          if(ppp(0)<swap_prob){//In case the swap is accepted
            //Swap vectors in the matrix X
            X.col(t+1)=Xtemp_from;
            X.col(t)=Xtemp_to;
          }
        }
      }
    }// Finish burn-in period
    swap_count=0; //Reset swap count
    
    std::clock_t start = std::clock(); // Start timer for simulation s
    //// Start the loop for all iterations in simulation s
    for(int i=0;i<numiter;i++){
      // Rcpp::Rcout <<"Inside iteration loop"<< i << std::endl;
      if (i % 1000 == 1) {Rcpp::Rcout << "PT A-IITw - Simulation: " << s+startsim << " Iteration: " << i <<"...Prob_decrease_bound: "<<prob_to_dec<< std::endl;}
      // Rcpp::Rcout << "Simulation: " << s+startsim << " Iteration: " << i << std::endl;
      for(int replica=0;replica<T;replica++){//For loop for replicas
        current_temp=temp(index_process(replica));// Extract temperature of the replica
        current_log_bound=log_bound_vector(index_process(replica));// Extract log-bound of the corresponding temperature
        
        if(update_prob){//Check if we need to update the probability
          if(current_temp==1){//The original replica defines the speed to modify the probabilty to decrease bounding constant
            sample_iterations_count=i; //Add the number of iterations (or samples) from the previous step
            if(sample_iterations_count>total_replica_iterations*percentage_start){//Check if we start decreasing the probability
              if(sample_iterations_count>total_replica_iterations*percentage_end){//Check if we stop decreasing the probability
                prob_to_dec=0;
                update_prob=false;
              }else{//In case we haven't finished updating the probability
                //Probability is proportional 
                double progress=sample_iterations_count/total_replica_iterations;
                prob_to_dec=1+((percentage_start-progress)/(percentage_end-percentage_start));
              }
            }
          }
        }
        output=a_IIT_update(X.col(replica),bal_function[index_process(replica)],current_temp,current_log_bound,bound_reduction,prob_to_dec,decreasing_constant);
        //// Store Z factor of replica with temperature 1
        if(current_temp==1){ // For the original temperature replica
          // Rcpp::Rcout << "Storing weight in simulation: " << s+startsim << " Iteration: " << i << std::endl;
          Z = output(1); //Extract the Z-factor
          // Rcpp::Rcout << "Printing Z: " << Z << std::endl;
          int state=vec_to_num(X.col(replica));
          // Rcpp::Rcout << "Printing state: " << state << std::endl;
          // pi_est(state)++;//Count how many times each state was visited
          pi_est(state)+=(1/Z);//Add weight
          if(first_visit(state)==0){
            first_visit(state)=i;//Store the first time the state is visited 
          }
          // Rcpp::Rcout << "Printing pi_est: " << pi_est << std::endl;
          // Rcpp::Rcout << "All good with Storing weight in simulation: " << s+startsim << " Iteration: " << i << std::endl;
        }
        X.col(replica)=vec(output(0)); //Update current state of the chain
        log_bound_vector(index_process(replica))=output(2); //Update log-bound 
      }//End loop to update replicas
      
      //// Start replica swap process
      if ((i+1) % iterswap == 0){
        swap_count+=1;//Increase the count of swaps
        // Rcpp::Rcout << "Trying swap: " << swap_count << std::endl;
        epsilon_indic.fill(-1); //Epsilon indic starts as -1
        prop_swap.zeros();
        do_swap.zeros();
        //Try a replica swap every iterswap number of iterations
        //We're doing non-reversible parallel tempering
        int starting=swap_count%2; // Detect if it's even or odd
        // Rcpp::Rcout <<"Trying replica swap "<<swap_count<<" start: "<<starting <<" at iteration: "<< i << std::endl;
        for(int t=starting;t<J;t+=2){// For loop that runs over temperature indexes to swap
          swap_total(t)+=1;//Increase the count of tried swaps
          epsilon_indic.elem(find(index_process==t)).ones();
          prop_swap.elem(find(index_process==t)).ones(); //we swap temperature t
          prop_swap.elem(find(index_process==t+1)).ones(); //With t+1
          //Compute swap probability
          Xtemp_from=X.cols(find(index_process==t));
          Xtemp_to=X.cols(find(index_process==t+1));
          
          
          //// Optional Computation of Z factors to correct bias
          double Z_fact_correc=1;//to temporarily store Z_factor correction
          if(bias_fix){
            double Z_temp11;
            double Z_temp12;
            double Z_temp21;
            double Z_temp22;
            

            output=a_IIT_update(Xtemp_from,bal_function[t],temp(t),log_bound_vector(t),false,prob_to_dec,decreasing_constant);
            Z_temp11=output(1);
            output=a_IIT_update(Xtemp_to,bal_function[t+1],temp(t+1),log_bound_vector(t+1),false,prob_to_dec,decreasing_constant);
            Z_temp22=output(1);
            output=a_IIT_update(Xtemp_from,bal_function[t+1],temp(t+1),log_bound_vector(t+1),false,prob_to_dec,decreasing_constant);
            Z_temp12=output(1);
            output=a_IIT_update(Xtemp_to,bal_function[t],temp(t),log_bound_vector(t),false,prob_to_dec,decreasing_constant);
            Z_temp21=output(1);
            
            Z_fact_correc=Z_temp12*Z_temp21/(Z_temp11*Z_temp22);
          }
          //// Computing swap probability
          swap_prob=(temp(t)-temp(t+1))*(loglik(Xtemp_to) - loglik(Xtemp_from)); 
          swap_prob=Z_fact_correc*exp(swap_prob);
          // Rcpp::Rcout <<"Swap prob "<< swap_prob << std::endl;
          ppp=Rcpp::runif(1);
          if(ppp(0)<swap_prob){//In case the swap is accepted
            swap_success(t)+=1;//Increase the number of successful swaps of temp t
            // Rcpp::Rcout <<"Accepted swap " << std::endl;
            do_swap.elem(find(index_process==t)).ones();
            do_swap.elem(find(index_process==t+1)).ones();
          }
        }
        resulting_swap=epsilon_indic % prop_swap % do_swap;
        // Rcpp::Rcout <<"Resulting swap "<< resulting_swap << std::endl;
        index_process+=resulting_swap;
        // Rcpp::Rcout <<"New index process:\n "<< index_process << std::endl;
        ind_pro_hist.row((s*total_swaps)+swap_count)=index_process.t();
        // Rcpp::Rcout <<"Store index process " << std::endl;
      }//End of replica swap process
    }// End loop of iterations
    std::clock_t end = std::clock(); // Stop timer
    // Calculate the time taken in seconds
    double duration = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    time_taken[s] = duration;
    // Store result of the simulation
    full_pi_est.col(s)=pi_est;
    full_first_visit.col(s)=first_visit;
    vec temp_rate=swap_success / swap_total;
    swap_rate.row(s)=temp_rate.t();
    // Rcpp::Rcout <<"Final state "<< X << std::endl;
  }//End loop simulations
  List ret;
  ret["est_pi"]=full_pi_est;
  ret["ip"]=ind_pro_hist;
  ret["visits"]=full_first_visit;
  ret["swap_rate"]=swap_rate;
  ret["time_taken"]=time_taken;
  return ret;
}

// Rcpp::Rcout <<"log-probs vector: \n"<< logprobs << std::endl;

