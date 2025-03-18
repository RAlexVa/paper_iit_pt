//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <fstream>
#include <ctime> 
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
// #include <iostream>
// #include <fstream>
// #include <armadillo>


////////// Other functions //////////

// Return maximum of 3 numbers
// [[Rcpp::export]]
double ret_max(double a,double b,double c){
  vec temp={a,b,c};
  return(max(temp));
}
// Transform vector representation of base 2 into decimal representation
// [[Rcpp::export]]
int vec_to_num(vec X){
  int n=X.n_rows;
  int number=0;
  for(int i=0;i<n;i++){
    if(X(i)==1){
      number+=std::pow(2,i); 
    }
  }
  return number;
}
// Transform number in decimal base to vector in base 2 representation
// [[Rcpp::export]]
vec num_to_vec(int n, int d){
  vec X(d);
  X.zeros();
  int temp;
  if((n<0) | (n>=std::pow(2,d))){
    Rcpp::Rcout <<"Error, number bigger than dimension " << std::endl;
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
// Take coordinates and create a binary vector with 1 in the defined coordinates

// [[Rcpp::export]]
vec createBinaryVector(const std::vector<int>& coordinates, int size) {
  vec binaryVector(size, fill::zeros); // Initialize vector with 0s
  for (int coord : coordinates) {
    if (coord >= 0 && coord < size) {
      binaryVector(coord) = 1; // Set the specified coordinates to 1
    }
  }
  return binaryVector;
}

// [[Rcpp::export]]
mat initializeMatrix(const std::vector<int>& coordinates, int size, int num_replicas) {
  // vec createBinaryVector(coordinates,size);
  mat ini_M(size,num_replicas);
  for (int coord : coordinates) {
    if (coord >= 0 && coord < size) {
      ini_M.row(coord).fill(1);
    }
  }
  return ini_M;
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


////////// Creating sparse matrix from file //////////
arma::sp_mat readSparseMatrix(const std::string& filename){
  // Open the file
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open file: " + filename);
  }
  
  // Read the first line to get the matrix dimensions and number of non-zero entries
  size_t n_rows, n_nonzeros;
  file >> n_rows >> n_nonzeros;
  
  // Create a sparse matrix
  arma::sp_mat matrix(n_rows, n_rows);
  arma::vec diag_vec(n_rows, fill::zeros);
  // Read the subsequent lines for the non-zero entries
  size_t row, col;
  double value;
  for (size_t i = 0; i < n_nonzeros; ++i) {
    file >> row >> col >> value;
    // Adjust for 1-based indexing in the file
    matrix(row - 1, col - 1) = -value;
    matrix(col - 1, row - 1) = -value;
    
    // Add for the diagonal
    diag_vec(row-1)+=1;
    diag_vec(col-1)+=1;
  }
  file.close();
  
  matrix.diag()=diag_vec;
  // Rcpp::Rcout <<"Matrix\n " <<matrix<< std::endl;
  return matrix;
}


// [[Rcpp::export]]
double loglik(const arma::vec& X,const arma::mat& M){
  // double theta1=10;
  double theta=200;
  // double theta2=200;
  if(M.n_cols!=2){
    Rcpp::Rcout << "Error matrix has more than 2 columns: " << std::endl;
    return(-10000);
  }else{
    vec mod1=M.col(0);
    vec mod2=M.col(1);
    double dif1=sum(abs(X-mod1));
    double dif2=sum(abs(X-mod2));
    double loglik_computed;
    if(dif2<=dif1){
      loglik_computed = -dif2/theta + log1p(exp((dif2-dif1)/theta));
    }
    if(dif2>dif1){
      loglik_computed = -dif1/theta + log1p(exp((dif1-dif2)/theta));
    }
    
    
    // loglik_computed = exp(-(sum(abs(X-mod1))/theta1))+exp(-(sum(abs(X-mod2))/theta2));
    // Rcpp::Rcout << "Primera parte: " <<-(sum(abs(X-mod1))/theta1)<< std::endl;
    // Rcpp::Rcout << "EXP Primera parte: " <<exp(-(sum(abs(X-mod1))/theta1))<< std::endl;
    // Rcpp::Rcout << "Segunda parte: " <<-(sum(abs(X-mod2))/theta2)<< std::endl;
    // Rcpp::Rcout << "EXP Segunda parte: " <<exp(-(sum(abs(X-mod2))/theta2))<< std::endl;
    return(loglik_computed);
  }
  
}




////////// Updating functions //////////


// [[Rcpp::export]]
List IIT_update_w(vec X,const arma::mat& M,String chosen_bf, double temperature){
  int total_neighbors = X.n_rows; // total number of neighbors is p spacial
  vec probs(total_neighbors, fill::zeros); //probabilities
  
  ////// Compute likelihood of current state
  double logpi_current=0;
  // uvec current_coord = find(X==1);
  logpi_current = loglik(X,M);
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
    temporal=loglik(newX,M)-logpi_current;
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
List a_IIT_update(vec X,const arma::mat& M, String chosen_bf, double temperature, double log_bound){
  int total_neighbors = X.n_rows; // total number of neighbors is p spacial
  vec probs(total_neighbors, fill::zeros); //probabilities
  ////// Compute likelihood of current state
  double logpi_current=0;
  logpi_current = loglik(X,M);
  ////// Compute weight for all neighbors
  double temporal=0;
  double temp_bound = 0;
  vec newX;
  for(int j=0; j<total_neighbors;j++){
    // Rcpp::Rcout << "Starts checking neighbors  "<< j<<std::endl; 
    newX = X;
    newX.row(j) = 1-X.row(j);
    //Rcpp::Rcout << newX << std::endl;
    temporal= bal_func(temperature*(loglik(newX,M)-logpi_current), chosen_bf);
    //Apply balancing function to log probability times temperature ////
    probs(j)=temporal;
    // Update bound if needed
    temp_bound=bal_func(temperature*(logpi_current-loglik(newX,M)), chosen_bf); //apply bf to the highest of pix-piy or piy-pix
    log_bound=ret_max(temporal,temp_bound,log_bound);
  }
  probs = probs - log_bound; // Apply bound to log-probabilities
  //////Choose the next neighbor
  vec u = Rcpp::runif(total_neighbors);
  vec probs_choose = log(-log(u)) - probs;
  
  //Find the index of the minimum element. 
  //This corresponds to choosing that neighbor
  int neigh_pos = (std::min_element(probs_choose.begin(), probs_choose.end()))-probs_choose.begin();
  // Rcpp::Rcout <<"probs vector: "<< probs << std::endl;
  // Rcpp::Rcout <<"chosen neighbor: "<< neigh_pos << std::endl;
  
  X.row(neigh_pos) = 1-X.row(neigh_pos); //modify the coordinate of the chosen neighbor
  List ret;
  ret["X"]=X;
  ret["Z"]=sum(exp(probs))/total_neighbors; // Compute Z factor with uniform proposal distribution
  ret["logbound"]=log_bound;
  return ret;
}

////////// Code for Parallel Tempering simulations //////////

// [[Rcpp::export]]
List PT_IIT_sim(int p,int startsim,int endsim, int numiter, int iterswap,int burn_in, vec temp, const std::vector<std::string>& bal_function, bool bias_fix,const std::string& filename,int num_states_visited,const std::vector<int>& starting_coord){
  //// Initialize variables to use in the code
  int T=temp.n_rows; // Count number of temperatures
  double J=double(T)-1;//Number of temperatures minus 1, used in swap loops
  int total_sim = (endsim-startsim+1); //Count total number of simulations
  int total_swaps=trunc(numiter/iterswap);
  List output; // To store output of the update function
  // double Z; // To store Z factor of update function
  int swap_count; //to keep track of even-odd swaps
  double current_temp; // to temporarily store the temperature
  //// Initialize arrays to store information
  mat X(p,T); // To store the current state of the joint chain, as many rows as neighbors, as many columns as temperatures
  vec index_process(T);   //Initialize index process vector
  mat ind_pro_hist(total_swaps*total_sim+1,T); //To store evolution of index process
  // int max_num=pow(2,p);
  
  
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
  // Variables to register visits to high probability states
  //Number of states to keep track
  mat iter_to_visit(num_states_visited,total_sim);
  mat loglikelihood_visited(num_states_visited,total_sim);
  loglikelihood_visited.fill(-10);//Initialize a very negative loglikelihood
  cube states_visited(p,num_states_visited,total_sim,fill::zeros);
  Rcpp::Rcout <<"p= "<<p<<" num_states_visited= "<<num_states_visited<<" total_sim= "<<total_sim<< std::endl;
  double temporal_loglik;
  uword found_min; // to find the minimum
  // mat modes_visited(numiter * total_sim,T);//Matrix to store the modes visited and temperature
  //// Define modes
  mat Q_matrix(p,2);
  for(int i=0;i<p;i++){
    if(i%2==0){Q_matrix(i,1)=1;}
    if(i%2==1){Q_matrix(i,0)=1;}
  }
  // Rcpp::Rcout << "First rows Q_matrix: " << Q_matrix.rows(0,5) << std::endl;
  // Rcpp::Rcout << "Last rows Q_matrix: " << Q_matrix.rows(p-6,p-1) << std::endl;
  
  std::vector<double> time_taken(total_sim); // vector to store the seconds each process took
  //// Start the loop for all simulations
  for(int s=0;s<total_sim;s++){
    for(int i=0;i<T;i++){ // Reset index process vector at the start of each simulation
      index_process.row(i)=i;
    }
    ind_pro_hist.row(0)=index_process.t(); // First entry of the index process
    swap_count=0; //Reset swap count
    X=initializeMatrix(starting_coord,p,T);//Reset the starting point of all chains
    swap_total.zeros();
    swap_success.zeros();
    //// Start loop for burn_in period
    for(int i=0;i<burn_in;i++){
      if (i % 100 == 1) {Rcpp::Rcout << "Simulation: " << s+startsim << " Burn_in period, iteration: " << i << std::endl;}
      for(int replica=0;replica<T;replica++){//For loop for replica update
        current_temp=temp(index_process(replica));
        output=IIT_update_w(X.col(replica),Q_matrix,bal_function[index_process(replica)],current_temp);
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
            
            output=IIT_update_w(Xtemp_from,Q_matrix,bal_function[t],temp(t));
            Z_temp11=output(1);
            output=IIT_update_w(Xtemp_to,Q_matrix,bal_function[t+1],temp(t+1));
            Z_temp22=output(1);
            output=IIT_update_w(Xtemp_from,Q_matrix,bal_function[t+1],temp(t+1));
            Z_temp12=output(1);
            output=IIT_update_w(Xtemp_to,Q_matrix,bal_function[t],temp(t));
            Z_temp21=output(1);
            
            Z_fact_correc=Z_temp12*Z_temp21/(Z_temp11*Z_temp22);
          }
          //// Computing swap probability
          swap_prob=(temp(t)-temp(t+1))*(loglik(Xtemp_to,Q_matrix) - loglik(Xtemp_from,Q_matrix)); 
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
        output=IIT_update_w(X.col(replica),Q_matrix,bal_function[index_process(replica)],current_temp);
        //// Store Z factor of replica with temperature 1
        if(current_temp==1){ // For the original temperature replica
          // Rcpp::Rcout << "Starts update of visited states" << std::endl;
          vec curr_loglik_visited=loglikelihood_visited.col(s);
          found_min=curr_loglik_visited.index_min();
          temporal_loglik=loglik(X.col(replica),Q_matrix);
          if(curr_loglik_visited(found_min)<temporal_loglik){
            Rcpp::Rcout << "Found big likelihood" <<exp(temporal_loglik)<<" in index: "<<found_min<< std::endl;
            // Rcpp::Rcout << "Updates state!\n with likelihood " <<curr_loglik_visited(found_min)<<" to loglik: "<<temporal_loglik<<" in poisition "<<found_min<< std::endl;
            loglikelihood_visited(found_min,s)=temporal_loglik;//Record new loglikelihood
            // Rcpp::Rcout << "Stores likelihood" << std::endl;
            iter_to_visit(found_min,s)=i;//Record iterations taken to visit the state
            // Rcpp::Rcout << "Stores iterations" << std::endl;
            //Record the new state
            for(int c=0;c<p;c++){
              // Rcpp::Rcout << "Stores entry "  <<c<<" of the cube"<< std::endl;
              states_visited(c,found_min,s)=X(c,replica);
            }
          }
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
            
            output=IIT_update_w(Xtemp_from,Q_matrix,bal_function[t],temp(t));
            Z_temp11=output(1);
            output=IIT_update_w(Xtemp_to,Q_matrix,bal_function[t+1],temp(t+1));
            Z_temp22=output(1);
            output=IIT_update_w(Xtemp_from,Q_matrix,bal_function[t+1],temp(t+1));
            Z_temp12=output(1);
            output=IIT_update_w(Xtemp_to,Q_matrix,bal_function[t],temp(t));
            Z_temp21=output(1);
            
            Z_fact_correc=Z_temp12*Z_temp21/(Z_temp11*Z_temp22);
            // Rcpp::Rcout <<"Z bias correction "<< Z_fact_correc << std::endl;
          }
          //// Computing swap probability
          swap_prob=(temp(t)-temp(t+1))*(loglik(Xtemp_to,Q_matrix) - loglik(Xtemp_from,Q_matrix)); 
          swap_prob=Z_fact_correc*exp(swap_prob);
          // Rcpp::Rcout <<"Swap prob "<< swap_prob << std::endl;
          ppp=Rcpp::runif(1);
          
          ///// A lot of comments           
          // Rcpp::Rcout <<"loglik from: "<< loglik(Xtemp_from,Q_matrix)<<"temp: "<<temp(t)<< std::endl;
          // Rcpp::Rcout <<"loglik to: "<< loglik(Xtemp_to,Q_matrix)<<"temp: "<<temp(t+1)<< std::endl;
          // Rcpp::Rcout <<"Difference in loglik: "<< (loglik(Xtemp_to,Q_matrix) - loglik(Xtemp_from,Q_matrix))<<std::endl;
          // Rcpp::Rcout <<"Difference in temp: "<< (temp(t)-temp(t+1))<<std::endl;
          // Rcpp::Rcout <<"Multiplying numbers from above: "<< (temp(t)-temp(t+1))*(loglik(Xtemp_to,Q_matrix) - loglik(Xtemp_from,Q_matrix))<<std::endl;
          // Rcpp::Rcout <<"Z correction "<< Z_fact_correc<< std::endl;
          // Rcpp::Rcout <<"Swap prob "<< swap_prob <<" RN: "<<ppp<< std::endl;
          
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
    vec temp_rate=swap_success / swap_total;
    swap_rate.row(s)=temp_rate.t();
    // Rcpp::Rcout <<"Final state "<< X << std::endl;
  }//End loop simulations
  List ret;
  
  ret["ip"]=ind_pro_hist;
  ret["swap_rate"]=swap_rate;
  ret["states"]=states_visited;
  ret["loglik_visited"]=loglikelihood_visited;
  ret["iter_visit"]=iter_to_visit;
  ret["time_taken"]=time_taken;
  return ret;
}

// [[Rcpp::export]]
List PT_a_IIT_sim(int p,int startsim,int endsim, int total_swaps,int sample_inter_swap,int burn_in, vec temp, const std::vector<std::string>& bal_function,const std::string& filename,int num_states_visited,const std::vector<int>& starting_coord){
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
  // int max_num=pow(2,p);
  // Rcpp::Rcout << "max num: " << max_num << std::endl;  
  
  
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
  // Variables to register visits to high probability states
  //Number of states to keep track
  mat iter_to_visit(num_states_visited,total_sim);
  mat loglikelihood_visited(num_states_visited,total_sim);
  loglikelihood_visited.fill(-10);//Initialize a very negative loglikelihood
  cube states_visited(p,num_states_visited,total_sim,fill::zeros);
  Rcpp::Rcout <<"p= "<<p<<" num_states_visited= "<<num_states_visited<<" total_sim= "<<total_sim<< std::endl;
  double temporal_loglik;
  uword found_min; // to find the minimum
  //// Define modes
  mat Q_matrix(p,2);
  for(int i=0;i<p;i++){
    if(i%2==0){Q_matrix(i,1)=1;}
    if(i%2==1){Q_matrix(i,0)=1;}
  }
  // Rcpp::Rcout << "First rows Q_matrix: " << Q_matrix.rows(0,5) << std::endl;
  // Rcpp::Rcout << "Last rows Q_matrix: " << Q_matrix.rows(p-6,p-1) << std::endl;
  std::vector<double> time_taken(total_sim); // vector to store the seconds each process took
  //// Start the loop for all simulations
  for(int s=0;s<total_sim;s++){
    for(int i=0;i<T;i++){ // Reset index process vector at the start of each simulation
      index_process.row(i)=i;
    }
    ind_pro_hist.row(0)=index_process.t(); // First entry of the index process
    swap_count=0; //Reset swap count
    X=initializeMatrix(starting_coord,p,T);//Reset the starting point of all chains
    
    
    log_bound_vector.zeros();//Reset log-bounds, all log-bounds start at 0
    swap_success.zeros();
    ////Start the loop for burn-in period
    int track_burn_in=0;
    while(track_burn_in<burn_in){
      for(int replica=0;replica<T;replica++){//For loop for replica update in the burn-in
        int samples_replica=0;
        while(samples_replica<sample_inter_swap){//Loop to create samples for each replica until we reach the defined threshold
          current_temp=temp(replica);// Extract temperature of the replica
          current_log_bound=log_bound_vector(replica);// Extract log-bound of the corresponding temperature
          output=a_IIT_update(X.col(replica),Q_matrix,bal_function[index_process(replica)],current_temp,current_log_bound);
          
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
        swap_prob=(temp(t)-temp(t+1))*(loglik(Xtemp_to,Q_matrix) - loglik(Xtemp_from,Q_matrix)); 
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
      Rcpp::Rcout << "Simulation: " << s+startsim << ". Done " << track_burn_in <<" samples in burn-in period"<< std::endl;
    }
    ////Finish the loop for burn-in period
    swap_count=0; //Reset swap count
    std::clock_t start = std::clock(); // Start timer for simulation s
    //// Start the loop for all iterations in simulation s
    for(int i=0;i<total_swaps;i++){
      // Rcpp::Rcout <<"Inside iteration loop"<< i << std::endl;
      if (i % 10 == 1) {Rcpp::Rcout << "Simulation: " << s+startsim << " Swap: " << i << std::endl;}
      // Rcpp::Rcout << "Simulation: " << s+startsim << " Iteration: " << i << std::endl;
      for(int replica=0;replica<T;replica++){//For loop for replicas
        // Rcpp::Rcout << "Sim: " << s+startsim << " Swap: " << i <<"replica: "<<replica<< std::endl;
        int samples_replica=0;
        while(samples_replica<sample_inter_swap){//Loop to create samples for each replica until we reach the defined threshold
          total_iterations(i,index_process(replica),s)+=1;//increase the number of iterations
          current_temp=temp(index_process(replica));// Extract temperature of the replica
          current_log_bound=log_bound_vector(index_process(replica));// Extract log-bound of the corresponding temperature
          output=a_IIT_update(X.col(replica),Q_matrix,bal_function[index_process(replica)],current_temp,current_log_bound);
          
          
          //// Compute weight
          Z = output(1); //Extract the Z-factor
          new_samples=1+R::rgeom(Z);
          if(new_samples<1){
            Rcpp::Rcout <<"Error: geometric in "<< "simulation: " << s+startsim << " Swap: " << i <<" temperature:"<<current_temp<< std::endl;
            Rcpp::Rcout <<"new_samples= "<<new_samples<< ", Z=" << Z << " log-bound= " << current_log_bound << std::endl;
            new_samples=sample_inter_swap;
          }
          if(samples_replica+new_samples>sample_inter_swap){//If we're going to surpass the required number of samples
            new_samples = sample_inter_swap-samples_replica;//Wee force to stop at sample_inter_swap
          }
          samples_replica+=new_samples; // Update number of samples obtained from the replica
          //// Store weight of replica with temperature 1
          if(current_temp==1){ // For the original temperature replica
            // Rcpp::Rcout << "Starts update of visited states" << std::endl;
            vec curr_loglik_visited=loglikelihood_visited.col(s);
            found_min=curr_loglik_visited.index_min();
            temporal_loglik=loglik(X.col(replica),Q_matrix);
            if(curr_loglik_visited(found_min)<temporal_loglik){
              Rcpp::Rcout << "Found big likelihood" <<exp(temporal_loglik)<<" in index: "<<found_min<< std::endl;
              // Rcpp::Rcout << "Updates state!\n with likelihood " <<curr_loglik_visited(found_min)<<" to loglik: "<<temporal_loglik<<" in poisition "<<found_min<< std::endl;
              loglikelihood_visited(found_min,s)=temporal_loglik;//Record new loglikelihood
              // Rcpp::Rcout << "Stores likelihood" << std::endl;
              mat current_slice=total_iterations.slice(s);//Extract current slice
              //Record iterations taken to visit the state
              iter_to_visit(found_min,s)=sum(current_slice.col(index_process(replica)));
              //Record the new state
              for(int c=0;c<p;c++){
                // Rcpp::Rcout << "Stores entry "  <<c<<" of the cube"<< std::endl;
                states_visited(c,found_min,s)=X(c,replica);
              }
            }
          }
          
          
          X.col(replica)=vec(output(0)); //Update current state of the chain
          log_bound_vector(index_process(replica))=output(2); //Update log-bound 
          // Rcpp::Rcout << "New neighbor: " << find(X.col(replica)==1) << std::endl;
          
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
        swap_prob=(temp(t)-temp(t+1))*(loglik(Xtemp_to,Q_matrix) - loglik(Xtemp_from,Q_matrix)); 
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
    vec temp_rate=swap_success / swap_total;
    swap_rate.row(s)=temp_rate.t();
    // Rcpp::Rcout <<"Final state "<< X << std::endl;
  }//End loop simulations
  List ret;
  
  ret["ip"]=ind_pro_hist;
  ret["swap_rate"]=swap_rate;
  ret["states"]=states_visited;
  ret["loglik_visited"]=loglikelihood_visited;
  ret["iter_visit"]=iter_to_visit;
  ret["total_iter"]=total_iterations;
  ret["time_taken"]=time_taken;
  return ret;
}

// [[Rcpp::export]]
List PT_a_IIT_sim_RF(int p,int startsim,int endsim, int numiter, int iterswap,int burn_in, vec temp, const std::vector<std::string>& bal_function, bool bias_fix,const std::string& filename,int num_states_visited,const std::vector<int>& starting_coord){
  //// Initialize variables to use in the code
  int T=temp.n_rows; // Count number of temperatures
  vec log_bound_vector(T); // vector to store a log-bound for each replica
  double J=double(T)-1;//Number of temperatures minus 1, used in swap loops
  int total_sim = (endsim-startsim+1); //Count total number of simulations
  int total_swaps=trunc(numiter/iterswap);
  List output; // To store output of the update function
  // double Z; // To store Z factor of update function
  int swap_count; //to keep track of even-odd swaps
  double current_temp; // to temporarily store the temperature
  double current_log_bound; //to temporarily store the log-bound
  //// Initialize arrays to store information
  mat X(p,T); // To store the current state of the joint chain, as many rows as neighbors, as many columns as temperatures
  vec index_process(T);   //Initialize index process vector
  mat ind_pro_hist(total_swaps*total_sim+1,T); //To store evolution of index process
  // int max_num=pow(2,p);
  
  
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
  // Variables to register visits to high probability states
  //Number of states to keep track
  mat iter_to_visit(num_states_visited,total_sim);
  mat loglikelihood_visited(num_states_visited,total_sim);
  loglikelihood_visited.fill(-10);//Initialize a very negative loglikelihood
  cube states_visited(p,num_states_visited,total_sim,fill::zeros);
  Rcpp::Rcout <<"p= "<<p<<" num_states_visited= "<<num_states_visited<<" total_sim= "<<total_sim<< std::endl;
  double temporal_loglik;
  uword found_min; // to find the minimum
  // mat modes_visited(numiter * total_sim,T);//Matrix to store the modes visited and temperature
  //// Define modes
  mat Q_matrix(p,2);
  for(int i=0;i<p;i++){
    if(i%2==0){Q_matrix(i,1)=1;}
    if(i%2==1){Q_matrix(i,0)=1;}
  }
  // Rcpp::Rcout << "First rows Q_matrix: " << Q_matrix.rows(0,5) << std::endl;
  // Rcpp::Rcout << "Last rows Q_matrix: " << Q_matrix.rows(p-6,p-1) << std::endl;
  std::vector<double> time_taken(total_sim); // vector to store the seconds each process took
  //// Start the loop for all simulations
  for(int s=0;s<total_sim;s++){
    for(int i=0;i<T;i++){ // Reset index process vector at the start of each simulation
      index_process.row(i)=i;
    }
    ind_pro_hist.row(0)=index_process.t(); // First entry of the index process
    swap_count=0; //Reset swap count
    X=initializeMatrix(starting_coord,p,T);//Reset the starting point of all chains
    
    
    swap_total.zeros();
    swap_success.zeros();
    log_bound_vector.zeros();//Reset log-bounds, all log-bounds start at 0
    //// Start loop for burn_in period
    for(int i=0;i<burn_in;i++){
      if (i % 100 == 1) {Rcpp::Rcout << "Simulation: " << s+startsim << " Burn_in period, iteration: " << i << std::endl;}
      for(int replica=0;replica<T;replica++){//For loop for replica update
        current_temp=temp(index_process(replica));
        current_log_bound=log_bound_vector(replica);// Extract log-bound of the corresponding temperature
        // output=IIT_update_w(X.col(replica),Q_matrix,bal_function[index_process(replica)],current_temp);
        output=a_IIT_update(X.col(replica),Q_matrix,bal_function[index_process(replica)],current_temp,current_log_bound);
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
            
            // output=IIT_update_w(Xtemp_from,Q_matrix,bal_function[t],temp(t));
            output=a_IIT_update(Xtemp_from,Q_matrix,bal_function[t],temp(t),log_bound_vector(t));
            Z_temp11=output(1);
            // output=IIT_update_w(Xtemp_to,Q_matrix,bal_function[t+1],temp(t+1));
            output=a_IIT_update(Xtemp_to,Q_matrix,bal_function[t+1],temp(t+1),log_bound_vector(t+1));
            Z_temp22=output(1);
            // output=IIT_update_w(Xtemp_from,Q_matrix,bal_function[t+1],temp(t+1));
            output=a_IIT_update(Xtemp_from,Q_matrix,bal_function[t+1],temp(t+1),log_bound_vector(t+1));
            Z_temp12=output(1);
            // output=IIT_update_w(Xtemp_to,Q_matrix,bal_function[t],temp(t));
            output=a_IIT_update(Xtemp_to,Q_matrix,bal_function[t],temp(t),log_bound_vector(t));
            Z_temp21=output(1);
            
            Z_fact_correc=Z_temp12*Z_temp21/(Z_temp11*Z_temp22);
          }
          //// Computing swap probability
          swap_prob=(temp(t)-temp(t+1))*(loglik(Xtemp_to,Q_matrix) - loglik(Xtemp_from,Q_matrix)); 
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
        current_log_bound=log_bound_vector(index_process(replica));// Extract log-bound of the corresponding temperature
        //Depending on the chosen method
        //// Update each replica independently
        // output=IIT_update_w(X.col(replica),Q_matrix,bal_function[index_process(replica)],current_temp);
        output=a_IIT_update(X.col(replica),Q_matrix,bal_function[index_process(replica)],current_temp,current_log_bound);
        //// Store Z factor of replica with temperature 1
        if(current_temp==1){ // For the original temperature replica
          // Rcpp::Rcout << "Starts update of visited states" << std::endl;
          vec curr_loglik_visited=loglikelihood_visited.col(s);
          found_min=curr_loglik_visited.index_min();
          temporal_loglik=loglik(X.col(replica),Q_matrix);
          if(curr_loglik_visited(found_min)<temporal_loglik){
            // Rcpp::Rcout << "Updates state!\n with likelihood " <<curr_loglik_visited(found_min)<<" to loglik: "<<temporal_loglik<<" in poisition "<<found_min<< std::endl;
            loglikelihood_visited(found_min,s)=temporal_loglik;//Record new loglikelihood
            // Rcpp::Rcout << "Stores likelihood" << std::endl;
            iter_to_visit(found_min,s)=i;//Record iterations taken to visit the state
            // Rcpp::Rcout << "Stores iterations" << std::endl;
            //Record the new state
            for(int c=0;c<p;c++){
              // Rcpp::Rcout << "Stores entry "  <<c<<" of the cube"<< std::endl;
              states_visited(c,found_min,s)=X(c,replica);
            }
          }
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
            
            // output=IIT_update_w(Xtemp_from,Q_matrix,bal_function[t],temp(t));
            output=a_IIT_update(Xtemp_from,Q_matrix,bal_function[t],temp(t),log_bound_vector(t));
            Z_temp11=output(1);
            // output=IIT_update_w(Xtemp_to,Q_matrix,bal_function[t+1],temp(t+1));
            output=a_IIT_update(Xtemp_to,Q_matrix,bal_function[t+1],temp(t+1),log_bound_vector(t+1));
            Z_temp22=output(1);
            // output=IIT_update_w(Xtemp_from,Q_matrix,bal_function[t+1],temp(t+1));
            output=a_IIT_update(Xtemp_from,Q_matrix,bal_function[t+1],temp(t+1),log_bound_vector(t+1));
            Z_temp12=output(1);
            // output=IIT_update_w(Xtemp_to,Q_matrix,bal_function[t],temp(t));
            output=a_IIT_update(Xtemp_to,Q_matrix,bal_function[t],temp(t),log_bound_vector(t));
            Z_temp21=output(1);
            
            Z_fact_correc=Z_temp12*Z_temp21/(Z_temp11*Z_temp22);
            // Rcpp::Rcout <<"Z bias correction "<< Z_fact_correc << std::endl;
          }
          //// Computing swap probability
          swap_prob=(temp(t)-temp(t+1))*(loglik(Xtemp_to,Q_matrix) - loglik(Xtemp_from,Q_matrix)); 
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
    vec temp_rate=swap_success / swap_total;
    swap_rate.row(s)=temp_rate.t();
    // Rcpp::Rcout <<"Final state "<< X << std::endl;
  }//End loop simulations
  List ret;
  
  ret["ip"]=ind_pro_hist;
  ret["swap_rate"]=swap_rate;
  ret["states"]=states_visited;
  ret["loglik_visited"]=loglikelihood_visited;
  ret["iter_visit"]=iter_to_visit;
  ret["time_taken"]=time_taken;
  return ret;
}


