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

// [[Rcpp::export]]
double tvd_compute(vec dist1, vec dist2){
  if(dist1.n_elem!=dist2.n_elem){Rcpp::Rcout <<"Vectors are not the same size" << std::endl; return -1;}
  if(sum(dist1)!=1){dist1 = dist1/sum(dist1);}
  if(sum(dist2)!=1){dist2 = dist2/sum(dist2);}
  
  double sum_diff = sum(abs(dist1-dist2));
  
  return 0.5*sum_diff;
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
    Rcpp::Rcout <<"Name of the balancing function does exist!" << std::endl;
    return 0; // Default return for unknown operation
  }
}


////////// loglikelihood functions //////////


// bi-modal log-likelihood, dimension 16
// [[Rcpp::export]]
double loglik(const arma::vec& X){
  int size=X.n_rows;
  if(size!=16){
    Rcpp::Rcout <<"Error:Size is not 16" << std::endl;
    return -std::numeric_limits<double>::infinity();
  }
  
  double theta1=6;
  double theta2=6;
  // Defined modes
  // arma::vec  mod1 = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};// 16 1s
  arma::vec  mod2 = {1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0};// 8 1s
  arma::vec  mod3 = {0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1};// 8 1s
  // arma::vec  mod4 = {1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0};// 8 1s
  // arma::vec  mod5 = {0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1};// 8 1s
  // arma::vec  mod6 = {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1};// 2 1s
  // arma::vec  mod7 = {0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0};// 2 1s
  
  double loglik_comp=0;
  
  // loglik_comp+=exp(-theta*sum(abs(X-mod1)));
  loglik_comp+=exp(-theta1*sum(abs(X-mod2)));
  loglik_comp+=exp(-theta2*sum(abs(X-mod3)));
  // loglik_comp+=exp(-theta*sum(abs(X-mod4)));
  // loglik_comp+=exp(-theta*sum(abs(X-mod5)));
  // loglik_comp+=exp(-theta*sum(abs(X-mod6)));
  // loglik_comp+=exp(-theta*sum(abs(X-mod7)));
  
  return log(loglik_comp);
}

// [[Rcpp::export]]
vec compute_true_dist(int p){
  double theta=15;
  vec temporal_state(p);
  vec true_dist(pow(2,p));
  for(int i=0; i<pow(2,p); i++){
    temporal_state = num_to_vec(i,p);
    true_dist(i)=loglik(temporal_state);
  }
  return exp(true_dist)/sum(exp(true_dist));
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
List a_IIT_update(vec X, String chosen_bf, double temperature, double log_bound, bool update, double prob_to_dec, double decreasing_constant, double max_logbound_found){
  const double threshold = 1e-5;//Threshold for updating maximum bound. It has to change at least this
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
  // Rcpp::Rcout <<"Max log probs: "<< bal_func(max(max_logprobs),"sq")<<std::endl;
  //Updating the log-bound
  
  
  
  if(update){//If it's defined to keep updating the bounding constant
    // Rcpp::Rcout <<"Temp: "<< temperature<<" C_max_bound: "<<max_logbound_found<<" new max bound: "<<bal_func(max(max_logprobs),"sq")<<std::endl;
    double checking_max_logprob=bal_func(max(max_logprobs),"sq");
    if(max_logbound_found<checking_max_logprob && (checking_max_logprob-max_logbound_found)>=threshold){
      // Rcpp::Rcout <<"Diff en log bound: "<< checking_max_logprob-max_logbound_found<<std::endl;
      // Rcpp::Rcout <<"Previous max log bound: "<< max_logbound_found*1000000<<std::endl;
      max_logbound_found=checking_max_logprob;
      // Rcpp::Rcout <<"New max log bound: "<< max_logbound_found*1000000<<std::endl;
      log_bound=max_logbound_found;//First update according to the maximum 
    }
    //Then try Arithmetic reduction of the bounding constant
    if(prob_to_dec>0){//If we consider a probability to decrease the constant
      double test_prob=0;
      if(prob_to_dec<1){
        vec ppp=Rcpp::runif(1);//Get a random number
        test_prob=ppp(0);
      }
      if(test_prob<prob_to_dec){//In case the update is accepted
        double temporal_log_b=log_bound/temperature;
        double delta_bound = decreasing_constant/exp(temporal_log_b);
        //log(a-b)=log(a)+log(1-b/a) ≈ log(a) - b/a if b/a is very small
        //a=exp(log_bound), b=delta_bound
        // Rcpp::Rcout<<"Temp: "<< temperature<< " Decreasing delta: "<< delta_bound<<" C_bound: "<<exp(log_bound)<<" C_log_bound: "<<log_bound<<" temporal_log_b: "<<temporal_log_b<<std::endl;
        // if(delta_bound<.06){
        //   log_bound=temperature*(temporal_log_b - delta_bound);
        // }else{
        //   log_bound = temperature*log1p(exp(temporal_log_b)-(decreasing_constant)); //Reduce constant
        //   // log_bound = temperature*log(exp(temporal_log_b)-(decreasing_constant)); //Reduce constant
        // }
        //Second option to update the log_bound
        
        log_bound = temperature*(temporal_log_b + log1p(-delta_bound));
        
        // Rcpp::Rcout <<"New log_b: "<< log_bound<<std::endl;
        // if(temperature==0.05){
        //   Rcpp::Rcout <<"Decreasing delta: "<< delta_bound<<" Current bound: "<<exp(log_bound)<<" new bound: "<<log_bound<<std::endl;
        //   Rcpp::Rcout <<"Decreasing log-bound to "<< log_bound<<std::endl;
        // }
        if(log_bound<0){log_bound=0;} //Minimum bound is 1, minimum log_bound is 0
        
      }else{
        // Rcpp::Rcout <<"Rejected a bounding constant decrease"<< std::endl;
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
  ret["max_logbound"]=max_logbound_found;
  return ret;
}


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
  vec first_time(max_num);//Vector to store the time of first visit to each state
  mat full_first_time(max_num,total_sim); //Matrix to store times of first visits considering all simulations
  std::clock_t visit_time; // To store time of visiting each state
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
  
  vec true_distribution = compute_true_dist(p); // Generate the true target distribution for this problem
  int tvd_measurements;//Define number of times to check the TVD
  int measure_this_swap;
  int tvd_swap_count;
  uword tvd_swap_index;
  if(total_swaps>=20){
    tvd_measurements=20;
    measure_this_swap=trunc(total_swaps/tvd_measurements); //define after how many swaps there will be a measurement
  }else{
    tvd_measurements=total_swaps;
    measure_this_swap=1; //define after how many swaps there will be a measurement
  }//Check that we have enough swaps to measure
  mat tvd_report(total_sim,tvd_measurements); //Create a matrix to store the tvd measurements
  mat tvd_time_report(total_sim,tvd_measurements); //Create a matrix to store the time of tvd measurements    
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
    tvd_swap_count=measure_this_swap;//Re-start the counter to check after which swaps to measure TVD
    tvd_swap_index=0;
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
            visit_time = std::clock(); //Time to visit this state
            first_time(state)=static_cast<double>(visit_time - start) / CLOCKS_PER_SEC;
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
      //// Include the measurement of TVD after some swaps      
      if(swap_count == tvd_swap_count && tvd_swap_index<tvd_measurements){
        Rcpp::Rcout <<"Measuring TVD in swap: "<< swap_count <<" out of "<<total_swaps<<", measured every: "<<measure_this_swap<<" swaps "<< std::endl;
        tvd_report(s,tvd_swap_index)=tvd_compute(true_distribution, pi_est);
        std::clock_t time_tvd_measurement = std::clock(); // Stop timer
        // Calculate the time taken in seconds
        double dur_tvd = static_cast<double>(time_tvd_measurement - start) / CLOCKS_PER_SEC;
        tvd_time_report(s,tvd_swap_index)=dur_tvd;
        tvd_swap_count+=measure_this_swap;
        tvd_swap_index++;
      }      
    }// End loop of iterations
    std::clock_t end = std::clock(); // Stop timer
    // Calculate the time taken in seconds
    double duration = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    time_taken[s] = duration;
    // Store result of the simulation
    full_pi_est.col(s)=pi_est;
    full_first_visit.col(s)=first_visit;
    full_first_time.col(s)=first_time;
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
  ret["time_visit"]=full_first_time;
  ret["tvd_report"]=tvd_report;
  ret["tvd_time_report"]=tvd_time_report;
  return ret;
}


// [[Rcpp::export]]
List PT_a_IIT_sim(int p,int startsim,int endsim, int total_swaps,int sample_inter_swap,int burn_in, vec temp, const std::vector<std::string>& bal_function, int initial_state, double decreasing_constant,std::string reduc_model){
  //// Initialize variables to use in the code
  int T=temp.n_rows; // Count number of temperatures
  vec log_bound_vector(T); // vector to store a log-bound for each replica
  vec max_log_bound_vector(T); // vector to store a the maximum log-bound for each replica
  double J=double(T)-1;//Number of temperatures minus 1, used in swap loops
  int total_sim = (endsim-startsim+1); //Count total number of simulations
  List output; // To store output of the update function
  double Z; // To store Z factor of update function
  int swap_count; //to keep track of even-odd swaps
  double current_temp; // to temporarily store the temperature
  double current_log_bound; //to temporarily store the log-bound
  int new_samples=0;//temporal variable to store weight with multiplicity list
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
  vec first_time(max_num);//Vector to store the time of first visit to each state
  mat full_first_time(max_num,total_sim); //Matrix to store times of first visits considering all simulations
  std::clock_t visit_time; // To store time of visiting each state
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
  bool update_prob=false; //To define if the probability to decrease the constant should decrease or not
  bool update_constant=true; //In case we want to stop the adapting process at some point
  double prob_to_dec=0;
  double percentage_start=0.05;
  double percentage_end=0.70;
  int total_replica_iterations=sample_inter_swap*total_swaps;
  int sample_iterations_count;
  if(reduc_model=="always"){prob_to_dec=1;}
  if(reduc_model=="never"){prob_to_dec=0;}
  if(reduc_model=="iterations"){update_prob=true;}
  
  vec true_distribution = compute_true_dist(p); // Generate the true target distribution for this problem
  int tvd_measurements;//Define number of times to check the TVD
  int measure_this_swap;
  int tvd_swap_count;
  uword tvd_swap_index;
  if(total_swaps>=20){
    tvd_measurements=20;
    measure_this_swap=trunc(total_swaps/tvd_measurements); //define after how many swaps there will be a measurement
  }else{
    tvd_measurements=total_swaps;
    measure_this_swap=1; //define after how many swaps there will be a measurement
  }//Check that we have enough swaps to measure
  mat tvd_report(total_sim,tvd_measurements); //Create a matrix to store the tvd measurements
  mat tvd_time_report(total_sim,tvd_measurements); //Create a matrix to store the time of tvd measurements    
  
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
    max_log_bound_vector.zeros();//Reset max log-bounds, all log-bounds start at 0
    swap_success.zeros();
    //Reset the probability to reduce the bounding constant
    if(reduc_model=="iterations"){update_prob=true;prob_to_dec=1;} //Reset the bool to update probability
    sample_iterations_count=0; // Reset the counting of iterations (or samples)
    tvd_swap_count=measure_this_swap;//Re-start the counter to check after which swaps to measure TVD
    tvd_swap_index=0;
    ////Start the loop for burn-in period
    int track_burn_in=0;
    while(track_burn_in<burn_in){
      for(int replica=0;replica<T;replica++){//For loop for replica update in the burn-in
        int samples_replica=0;
        while(samples_replica<sample_inter_swap){//Loop to create samples for each replica until we reach the defined threshold
          current_temp=temp(replica);// Extract temperature of the replica
          current_log_bound=log_bound_vector(replica);// Extract log-bound of the corresponding temperature
          //// In burn-in we update (increase) the constant but we don't decrease it.
          output=a_IIT_update(X.col(replica),bal_function[index_process(replica)],current_temp,current_log_bound,true,0,0,max_log_bound_vector(replica));
          //During burn-in:
          ////// Update = true, we always update the constant
          ////// prob_to_dec=0, we never decrease the constant 
          ////// decreasing constant=0, we decrease by 0 (redundancy)
          ////// we keep track of the max log bound found
          
          // Compute weight
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
          max_log_bound_vector(index_process(replica))=output(3); //Update MAX log-bound 
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
    // Rcpp::Rcout <<"After burn-in log-bound vector:\n "<< log_bound_vector << std::endl;
    // Rcpp::Rcout <<"After burn-in  MAX log-bound vector:\n "<< max_log_bound_vector << std::endl;
    std::clock_t start = std::clock(); // Start timer for simulation s
    //// Start the loop for all iterations in simulation s
    for(int i=0;i<total_swaps;i++){
      // Rcpp::Rcout <<"Inside iteration loop"<< i << std::endl;
      if (i % 10 == 1) {Rcpp::Rcout << "PT A-IITm - Simulation: " << s+startsim << " Swap: " << i <<" Prob_decrease_bound: " << prob_to_dec << std::endl;}
      // if (i % 10 == 1) {Rcpp::Rcout <<"Current log_bound vector :\n"<< log_bound_vector<< std::endl;}
      // Rcpp::Rcout <<"Current log_bound vector :\n"<< log_bound_vector<< std::endl;
      //   bool check_bool= log_bound_vector(J)==0;
      // Rcpp::Rcout <<"Check if bound is 0 already :"<< check_bool<< std::endl;
      //<< " log-bound:\n " << log_bound_vector
      // Rcpp::Rcout << "Simulation: " << s+startsim << " Iteration: " << i << std::endl;
      // Rcpp::Rcout <<"Iter "<<i << std::endl;
      for(int replica=0;replica<T;replica++){//For loop for replicas
        int samples_replica=0;
        while(samples_replica<sample_inter_swap){//Loop to create samples for each replica until we reach the defined threshold
          total_iterations(i,index_process(replica),s)+=1;//increase the number of iterations
          current_temp=temp(index_process(replica));// Extract temperature of the replica
          current_log_bound=log_bound_vector(index_process(replica));// Extract log-bound of the corresponding temperature
          ///// Process to update probability of decreasing the bounding constant
          // Rcpp::Rcout <<"Temp: "<<current_temp<<" Current log_bound: "<< current_log_bound<<" C_max_log_bound: "<<max_log_bound_vector(index_process(replica))<< std::endl;
          output=a_IIT_update(X.col(replica),bal_function[index_process(replica)],current_temp,current_log_bound,update_constant,prob_to_dec,decreasing_constant,max_log_bound_vector(index_process(replica)));
          
          //// Compute weight
          Z = output(1); //Extract the Z-factor
          new_samples=1+R::rgeom(Z);
          if(new_samples<1){
            Rcpp::Rcout <<"Error with geometric in "<< "Simulation: " << s+startsim << " Swap: " << i <<" temperature:"<<current_temp<< std::endl;
            new_samples=sample_inter_swap;
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
              visit_time = std::clock(); //Time to visit this state
              first_time(state)=static_cast<double>(visit_time - start) / CLOCKS_PER_SEC;
            }
            if(update_prob){//Check if we need to update the probability
              if(current_temp==1){//The original replica defines the speed to modify the probability to decrease bounding constant
                sample_iterations_count+=new_samples; //Add the number of iterations (or samples) from the previous step
                // Rcpp::Rcout <<" New samples: "<<new_samples<<" sample_iterations_count= "<<sample_iterations_count<< std::endl;
                // Rcpp::Rcout <<"Update prob samples: "<< sample_iterations_count <<" total iterations: "<<total_replica_iterations<< std::endl;
                if(sample_iterations_count>(total_replica_iterations*percentage_start)){//Check if we start decreasing the probability
                  if(sample_iterations_count>(total_replica_iterations*percentage_end)){//Check if we stop decreasing the probability
                    Rcpp::Rcout <<"Stop decreasing bounds. Last bound vector:\n"<< log_bound_vector<< std::endl;
                    prob_to_dec=0;
                    update_prob=false;
                  }else{//In case we haven't finished updating the probability
                    //Probability is proportional 
                    double progress = static_cast<double>(sample_iterations_count) / total_replica_iterations;
                    // Rcpp::Rcout <<""<< sample_iterations_count <<" / "<<total_replica_iterations<<"="<<progress<< std::endl;
                    // Rcpp::Rcout <<"progress: "<< progress << std::endl;
                    prob_to_dec=1+((percentage_start-progress)/(percentage_end-percentage_start));
                    // Rcpp::Rcout <<"New prob: "<< prob_to_dec << std::endl;
                  }
                }
              }
            } 
          }
          X.col(replica)=vec(output(0)); //Update current state of the chain
          log_bound_vector(index_process(replica))=output(2); //Update log-bound 
          max_log_bound_vector(index_process(replica))=output(3); //Update maximum log-bound found
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
      
      //// Include the measurement of TVD after some swaps      
      if(swap_count == tvd_swap_count && tvd_swap_index<tvd_measurements){
        Rcpp::Rcout <<"Measuring TVD in swap: "<< swap_count <<" out of "<<total_swaps<<", measured every: "<<measure_this_swap<<" swaps "<< std::endl;
        tvd_report(s,tvd_swap_index)=tvd_compute(true_distribution, pi_est);
        std::clock_t time_tvd_measurement = std::clock(); // Stop timer
        // Calculate the time taken in seconds
        double dur_tvd = static_cast<double>(time_tvd_measurement - start) / CLOCKS_PER_SEC;
        tvd_time_report(s,tvd_swap_index)=dur_tvd;
        tvd_swap_count+=measure_this_swap;
        tvd_swap_index++;
      }
    }// End loop of iterations
    std::clock_t end = std::clock(); // Stop timer
    // Calculate the time taken in seconds
    double duration = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    time_taken[s] = duration;
    // Store result of the simulation
    full_pi_est.col(s)=pi_est;
    full_first_visit.col(s)=first_visit;
    full_first_time.col(s)=first_time;
    vec temp_rate=swap_success / swap_total;
    swap_rate.row(s)=temp_rate.t();
    Rcpp::Rcout <<"Final log-bound vector:\n "<< log_bound_vector << std::endl;
    Rcpp::Rcout <<"MAX log-bound vector:\n "<< max_log_bound_vector << std::endl;
  }//End loop simulations
  List ret;
  ret["est_pi"]=full_pi_est;
  ret["ip"]=ind_pro_hist;
  ret["visits"]=full_first_visit;
  ret["swap_rate"]=swap_rate;
  ret["total_iter"]=total_iterations;
  ret["time_taken"]=time_taken;
  ret["time_visit"]=full_first_time;
  ret["max_bounds"]=max_log_bound_vector;
  ret["final_bounds"]=log_bound_vector;
  ret["tvd_report"]=tvd_report;
  ret["tvd_time_report"]=tvd_time_report;
  return ret;
}

// [[Rcpp::export]]
List PT_a_IIT_sim_RF(int p,int startsim,int endsim, int numiter,int iterswap,int burn_in, vec temp, const std::vector<std::string>& bal_function, bool bias_fix, int initial_state, double decreasing_constant,std::string reduc_model){
  //// Initialize variables to use in the code
  int T=temp.n_rows; // Count number of temperatures
  vec log_bound_vector(T); // vector to store a log-bound for each replica
  vec max_log_bound_vector(T); // vector to store the MAX log-bound found for each replica
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
  vec first_time(max_num);//Vector to store the time of first visit to each state
  mat full_first_time(max_num,total_sim); //Matrix to store times of first visits considering all simulations
  std::clock_t visit_time; // To store time of visiting each state
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
  bool update_constant=true;
  double prob_to_dec=0;
  double percentage_start=0.05;
  double percentage_end=0.70;
  int total_replica_iterations=numiter;
  int sample_iterations_count;
  if(reduc_model=="always"){prob_to_dec=1;}
  if(reduc_model=="never"){prob_to_dec=0;}
  if(reduc_model=="iterations"){update_prob=true;}
  
  vec true_distribution = compute_true_dist(p); // Generate the true target distribution for this problem
  int tvd_measurements;//Define number of times to check the TVD
  int measure_this_swap;
  int tvd_swap_count;
  uword tvd_swap_index;
  if(total_swaps>=20){
    tvd_measurements=20;
    measure_this_swap=trunc(total_swaps/tvd_measurements); //define after how many swaps there will be a measurement
  }else{
    tvd_measurements=total_swaps;
    measure_this_swap=1; //define after how many swaps there will be a measurement
  }//Check that we have enough swaps to measure
  mat tvd_report(total_sim,tvd_measurements); //Create a matrix to store the tvd measurements
  mat tvd_time_report(total_sim,tvd_measurements); //Create a matrix to store the time of tvd measurements    
  
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
    
    tvd_swap_count=measure_this_swap;//Re-start the counter to check after which swaps to measure TVD
    tvd_swap_index=0;
    //// Start loop for burn_in period
    for(int i=0;i<burn_in;i++){
      if (i % 100 == 1) {Rcpp::Rcout << "PT A-IITw - Simulation: " << s+startsim << " Burn_in period, iteration: " << i << std::endl;}
      for(int replica=0;replica<T;replica++){//For loop for replica update
        current_temp=temp(index_process(replica));
        current_log_bound=log_bound_vector(replica);// Extract log-bound of the corresponding temperature
        //// During burn-in we update the constant (increase) but we don't decrease it.
        output=a_IIT_update(X.col(replica),bal_function[index_process(replica)],current_temp,current_log_bound,true,0,0,max_log_bound_vector(replica));
        //During burn-in:
        ////// Update = true, we always update the constant
        ////// prob_to_dec=0, we never decrease the constant 
        ////// decreasing constant=0, we decrease by 0 (redundancy)
        ////// we keep track of the max log bound found
        X.col(replica)=vec(output(0)); //Update current state of the chain
        log_bound_vector(index_process(replica))=output(2); //Update log-bound 
        max_log_bound_vector(index_process(replica))=output(3); //Update log-bound 
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
            // For replica swaps we don't update the bounding constant
            output=a_IIT_update(Xtemp_from,bal_function[t],temp(t),log_bound_vector(t),false,0,0,max_log_bound_vector(t));
            Z_temp11=output(1);
            output=a_IIT_update(Xtemp_to,bal_function[t+1],temp(t+1),log_bound_vector(t+1),false,0,0,max_log_bound_vector(t+1));
            Z_temp22=output(1);
            output=a_IIT_update(Xtemp_from,bal_function[t+1],temp(t+1),log_bound_vector(t+1),false,0,0,max_log_bound_vector(t+1));
            Z_temp12=output(1);
            output=a_IIT_update(Xtemp_to,bal_function[t],temp(t),log_bound_vector(t),false,0,0,max_log_bound_vector(t));
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
        
        output=a_IIT_update(X.col(replica),bal_function[index_process(replica)],current_temp,current_log_bound,update_constant,prob_to_dec,decreasing_constant,max_log_bound_vector(index_process(replica)));
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
            visit_time = std::clock(); //Time to visit this state
            first_time(state)=static_cast<double>(visit_time - start) / CLOCKS_PER_SEC;
          }
          // Rcpp::Rcout << "Printing pi_est: " << pi_est << std::endl;
          // Rcpp::Rcout << "All good with Storing weight in simulation: " << s+startsim << " Iteration: " << i << std::endl;
          if(update_prob){//Check if we need to update the probability
            if(current_temp==1){//The original replica defines the speed to modify the probability to decrease bounding constant
              sample_iterations_count+=1; //Add the number of iterations (or samples) from the previous step
              // Rcpp::Rcout <<" New samples: "<<new_samples<<" sample_iterations_count= "<<sample_iterations_count<< std::endl;
              // Rcpp::Rcout <<"Update prob samples: "<< sample_iterations_count <<" total iterations: "<<total_replica_iterations<< std::endl;
              if(sample_iterations_count>(total_replica_iterations*percentage_start)){//Check if we start decreasing the probability
                if(sample_iterations_count>(total_replica_iterations*percentage_end)){//Check if we stop decreasing the probability
                  prob_to_dec=0;
                  update_prob=false;
                }else{//In case we haven't finished updating the probability
                  //Probability is proportional 
                  double progress = static_cast<double>(sample_iterations_count) / total_replica_iterations;
                  // Rcpp::Rcout <<""<< sample_iterations_count <<" / "<<total_replica_iterations<<"="<<progress<< std::endl;
                  // Rcpp::Rcout <<"progress: "<< progress << std::endl;
                  prob_to_dec=1+((percentage_start-progress)/(percentage_end-percentage_start));
                  // Rcpp::Rcout <<"New prob: "<< prob_to_dec << std::endl;
                }
              }
            }
          } 
        }
        X.col(replica)=vec(output(0)); //Update current state of the chain
        log_bound_vector(index_process(replica))=output(2); //Update log-bound 
        max_log_bound_vector(index_process(replica))=output(3); //Update log-bound 
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
            
            //// In replica swaps we don't update the bounding constant
            output=a_IIT_update(Xtemp_from,bal_function[t],temp(t),log_bound_vector(t),false,0,0,max_log_bound_vector(t));
            Z_temp11=output(1);
            output=a_IIT_update(Xtemp_to,bal_function[t+1],temp(t+1),log_bound_vector(t+1),false,0,0,max_log_bound_vector(t+1));
            Z_temp22=output(1);
            output=a_IIT_update(Xtemp_from,bal_function[t+1],temp(t+1),log_bound_vector(t+1),false,0,0,log_bound_vector(t+1));
            Z_temp12=output(1);
            output=a_IIT_update(Xtemp_to,bal_function[t],temp(t),log_bound_vector(t),false,0,0,log_bound_vector(t));
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
      //// Include the measurement of TVD after some swaps      
      if(swap_count == tvd_swap_count && tvd_swap_index<tvd_measurements){
        Rcpp::Rcout <<"Measuring TVD in swap: "<< swap_count <<" out of "<<total_swaps<<", measured every: "<<measure_this_swap<<" swaps "<< std::endl;
        tvd_report(s,tvd_swap_index)=tvd_compute(true_distribution, pi_est);
        std::clock_t time_tvd_measurement = std::clock(); // Stop timer
        // Calculate the time taken in seconds
        double dur_tvd = static_cast<double>(time_tvd_measurement - start) / CLOCKS_PER_SEC;
        tvd_time_report(s,tvd_swap_index)=dur_tvd;
        tvd_swap_count+=measure_this_swap;
        tvd_swap_index++;
      }     
    }// End loop of iterations
    std::clock_t end = std::clock(); // Stop timer
    // Calculate the time taken in seconds
    double duration = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    time_taken[s] = duration;
    // Store result of the simulation
    full_pi_est.col(s)=pi_est;
    full_first_visit.col(s)=first_visit;
    full_first_time.col(s)=first_time;
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
  ret["time_visit"]=full_first_time;
  ret["max_bounds"]=max_log_bound_vector;
  ret["final_bounds"]=log_bound_vector;
  ret["tvd_report"]=tvd_report;
  ret["tvd_time_report"]=tvd_time_report;
  return ret;
}

// Rcpp::Rcout <<"log-probs vector: \n"<< logprobs << std::endl;

///// 
