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

// [[Rcpp::export]]
double ret_min(double a,double b,double c){
  vec temp={a,b,c};
  return(min(temp));
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

// [[Rcpp::export]]
mat initializeRandom(const int& num_rows,const int& num_cols, const double& prob) {

  // Initialize a matrix with random values between 0 and 1
  arma::mat A = arma::randu<arma::mat>(num_rows, num_cols);
  
  // Threshold the random values to 0 or 1
  A = arma::conv_to<arma::mat>::from(A > prob);
  
  return(A);
}

// Compare two vectors
// [[Rcpp::export]]
bool CompareVectors(const vec& v1, const vec& v2) {
  // Check if the vectors have the same size
  if (v1.n_elem != v2.n_elem) {
    return false;
  }
  
  // Compare each element of the vectors
  for (uword i = 0; i < v1.n_elem; ++i) {
    if (v1(i) != v2(i)) {
      return false;
    }
  }
  // If all elements are the same, return true
  return true;
}

// [[Rcpp::export]]
int L1_distance(const vec& v1, const vec& v2) {
  // Check if the vectors have the same size
  if (v1.n_elem != v2.n_elem) {
    Rcpp::Rcout <<"Vectors are not the same size" << std::endl;
    return -1;
  }else{
    int dist=sum(abs(v1-v2));
    return dist;
  }
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
  double theta=3;
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
    // if(dif2<=dif1){
    //   loglik_computed = -dif2/theta + log1p(exp((dif2-dif1)/theta));
    // }
    // if(dif2>dif1){
    //   loglik_computed = -dif1/theta + log1p(exp((dif1-dif2)/theta));
    // }
    
    if(dif2<=dif1){
      loglik_computed = -dif2*theta + log1p(exp((dif2-dif1)*theta));
    }
    if(dif2>dif1){
      loglik_computed = -dif1*theta + log1p(exp((dif1-dif2)*theta));
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
List a_IIT_update(vec X,const arma::mat& M, String chosen_bf, const double& temperature, double log_bound, const bool& update, double prob_to_dec, const double& decreasing_constant, double max_logbound_found){
  int total_neighbors = X.n_rows; // total number of neighbors is p spacial
  double logpi_current= loglik(X,M);// likelihood of current state
  
  const double threshold = 1e-5;//Threshold for updating maximum bound. It has to change at least this
  vec logprobs(total_neighbors, fill::zeros); //probabilities
  vec max_logprobs(total_neighbors,fill::zeros);//vector to store max-log-probabilities
  double temporal=0;
  vec newX;


  
  
  
  if(update){//If it's defined to keep updating the bounding constant
    ////// Compute weight for all neighbors
    for(int j=0; j<total_neighbors;j++){
      newX = X;
      newX.row(j) = 1-X.row(j);//Change coordinate of the state to define new neighbor
      temporal= temperature*(loglik(newX,M)-logpi_current);
      logprobs(j)=temporal;//Store raw log_probability
      max_logprobs(j)=abs(temporal); // Store the max log-probability, either pi_y/pi_x or pi_x/pi_y
    }// End of loop to compute raw log-probability of neighbors
    
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
        //a=exp(temporal_log_b), b=decreasing_constant
        if(delta_bound>=1){//If the decreasing constant is too big
          //log(1-b/a) would be undefined
          log_bound=0;//Go the minimum 
        }else{
          log_bound = temperature*(temporal_log_b + log1p(-delta_bound));
        }
        
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
    //// Apply bounding B.F. with updated constant    
    for(int j=0;j<total_neighbors;j++){
      logprobs(j)=bound_sq(logprobs(j),log_bound);
    }
  }else{//In case we don't update the bouding constant
    
    //// Apply bounded B.F. without modifying the constant
    for(int j=0;j<total_neighbors;j++){
      newX = X;
      newX.row(j) = 1-X.row(j); //Change coordinate of the state to define new neighbor
      temporal = temperature*(loglik(newX,M)-logpi_current);
      logprobs(j)=bound_sq(temporal,log_bound);
    }
  }
  
  /////////////
  ////IMPORTANT: we're only using bound sqrt root for the adaptive IIT
  /////////////

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
  loglikelihood_visited.fill(-10000);//Initialize a very negative loglikelihood
  cube states_visited(p,num_states_visited,total_sim,fill::zeros);
  Rcpp::Rcout <<"p= "<<p<<" num_states_visited= "<<num_states_visited<<" total_sim= "<<total_sim<< std::endl;
  double temporal_loglik;
  uword found_min; // to find the minimum
  // mat modes_visited(numiter * total_sim,T);//Matrix to store the modes visited and temperature
  //// Define two modes
  mat Q_matrix(p,2);
  for(int i=0;i<p;i++){
    if(i%2==0){Q_matrix(i,1)=1;}
    if(i%2==1){Q_matrix(i,0)=1;}
  }
  vec mode1=Q_matrix.col(0);
  vec mode2=Q_matrix.col(1);
  mat distance_mode1(numiter,T);//Matrix to store distance to mode 1 for each replica
  mat distance_mode2(numiter,T);//Matrix to store distance to mode 2 for each replica
  mat distance_origin(numiter,T);//Matrix to store distance to origin for each replica
  cube full_distance_mode1(numiter,T,total_sim);//Cube to store distance to mode 1 for each replica and simulation
  cube full_distance_mode2(numiter,T,total_sim);//Cube to store distance to mode 2 for each replica and simulation
  cube full_distance_origin(numiter,T,total_sim);//Cube to store distance to origin for each replica and simulation
  
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
    // X=initializeMatrix(starting_coord,p,T);//Reset the starting point of all chains
    ppp=Rcpp::runif(1);
    X=initializeRandom(p,T,ppp(0));//Randomly initialize the state of each replica.
    swap_total.zeros();
    swap_success.zeros();
    distance_mode1.fill(-1);
    distance_mode2.fill(-1);
    distance_origin.fill(-1);
    //// Start loop for burn_in period
    for(int i=0;i<burn_in;i++){
      if (i % 100 == 1) {Rcpp::Rcout << "PT-IIT Simulation: " << s+startsim << " Burn_in period, iteration: " << i << std::endl;}
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
    }
/////////////////////// Finish burn-in period
    swap_count=0; //Reset swap count
    std::clock_t start = std::clock(); // Start timer for simulation s
    //// Start the loop for all iterations in simulation s
    for(int i=0;i<numiter;i++){
      // Rcpp::Rcout <<"Inside iteration loop"<< i << std::endl;
      if (i % 1000 == 1) {Rcpp::Rcout << "PT-IIT Simulation: " << s+startsim << " Iteration: " << i << std::endl;}
      // Rcpp::Rcout << "Simulation: " << s+startsim << " Iteration: " << i << std::endl;
      for(int replica=0;replica<T;replica++){//For loop for replicas
        current_temp=temp(index_process(replica));// Extract temperature of the replica
        //Measure distance to the modes
        distance_mode1(i,index_process(replica))=L1_distance(X.col(replica),mode1);
        distance_mode2(i,index_process(replica))=L1_distance(X.col(replica),mode2);
        distance_origin(i,index_process(replica))=sum(X.col(replica));
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
            Rcpp::Rcout << "Found big likelihood: " <<exp(temporal_loglik)<<" in index: "<<found_min<< std::endl;
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
    
    full_distance_mode1.slice(s)=distance_mode1;
    full_distance_mode2.slice(s)=distance_mode2;
    full_distance_origin.slice(s)=distance_origin;
    // Rcpp::Rcout <<"Final state "<< X << std::endl;
  }//End loop simulations
  List ret;
  
  ret["ip"]=ind_pro_hist;
  ret["swap_rate"]=swap_rate;
  ret["states"]=states_visited;
  ret["loglik_visited"]=loglikelihood_visited;
  ret["iter_visit"]=iter_to_visit;
  ret["time_taken"]=time_taken;
  // ret["modes_visit"]=full_iter_visit_modes;
  ret["distance_mode1"]=full_distance_mode1;
  ret["distance_mode2"]=full_distance_mode2;
  ret["distance_origin"]=full_distance_origin;
  return ret;
}

// [[Rcpp::export]]
List PT_a_IIT_sim(int p,int startsim,int endsim, int total_swaps,int sample_inter_swap,int burn_in, vec temp, const std::vector<std::string>& bal_function,const std::string& filename,int num_states_visited,const std::vector<int>& starting_coord, double decreasing_constant,std::string reduc_model){
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
  int new_samples;//temporal variable to store weight with multiplicity list
  //// Initialize arrays to store information
  mat X(p,T); // To store the current state of the joint chain, as many rows as neighbors, as many columns as temperatures
  vec index_process(T);   //Initialize index process vector
  mat ind_pro_hist(total_swaps*total_sim+1,T); //To store evolution of index process
  // int max_num=pow(2,p);
  // Rcpp::Rcout << "max num: " << max_num << std::endl;  
  vec current_state(p);//To print when I find a very small Z factor
  
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
  loglikelihood_visited.fill(-10000);//Initialize a very negative loglikelihood
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
  vec mode1=Q_matrix.col(0);
  vec mode2=Q_matrix.col(1);
  int numiter=total_swaps*sample_inter_swap;//Compute the total number of iterations to perform
  mat distance_mode1(numiter,T);//Matrix to store distance to mode 1 for each replica
  mat distance_mode2(numiter,T);//Matrix to store distance to mode 1 for each replica
  mat distance_origin(numiter,T);//Matrix to store distance to origin for each replica
  mat temporal_mat1(sample_inter_swap,T);//Matrix to store the distances computed between each replica swap
  mat temporal_mat2(sample_inter_swap,T);//Matrix to store the distances computed between each replica swap
  mat temporal_origin(sample_inter_swap,T);//Matrix to store the distances computed between each replica swap  
  cube full_distance_mode1(numiter,T,total_sim);//Cube to store distance to mode 1 for each replica and simulation
  cube full_distance_mode2(numiter,T,total_sim);//Cube to store distance to mode 2 for each replica and simulation
  cube full_distance_origin(numiter,T,total_sim);//Cube to store distance to origin for each replica and simulation
  
  std::vector<double> time_taken(total_sim); // vector to store the seconds each process took
  // Probability to update
  bool update_prob=false; //To define if the probability to decrease the constant should decrease or not
  bool update_constant=true; //In case we want to stop the adapting process at some point
  double prob_to_dec=0;
  double percentage_start=0.05;
  double percentage_end=0.70;
  int total_replica_iterations=sample_inter_swap*total_swaps;
  int sample_iterations_count;
  if(false){}
  else if(reduc_model=="always"){prob_to_dec=1;}
  else if(reduc_model=="never"){prob_to_dec=0;}
  else if(reduc_model=="iterations"){update_prob=true;}
  else if(reduc_model=="zero"){update_constant=false;}
  else {Rcpp::Rcout <<"reduc_model= " <<reduc_model<<" is not a valid reduc_model. Default to standard"<< std::endl;
    Rcpp::Rcout <<" The standard is: Bounding constant always increases."<< std::endl;
    prob_to_dec=0;}//If we don't define a reduc_model
  
  
  
  //// Start the loop for all simulations
  for(int s=0;s<total_sim;s++){
    for(int i=0;i<T;i++){ // Reset index process vector at the start of each simulation
      index_process.row(i)=i;
    }
    ind_pro_hist.row(0)=index_process.t(); // First entry of the index process
    swap_count=0; //Reset swap count
    // X=initializeMatrix(starting_coord,p,T);//Reset the starting point of all chains
    ppp=Rcpp::runif(1);
    X=initializeRandom(p,T,ppp(0));//Randomly initialize the state of each replica.
    
    log_bound_vector.zeros();//Reset log-bounds, all log-bounds start at 0
    swap_success.zeros();
    distance_mode1.fill(-1);
    distance_mode2.fill(-1);
    distance_origin.fill(-1);
    //Reset the probability to reduce the bounding constant
    if(reduc_model=="iterations"){update_prob=true;prob_to_dec=1;} //Reset the bool to update probability
    sample_iterations_count=0; // Reset the counting of iterations (or samples)
    ////Start the loop for burn-in period
    Rcpp::Rcout << "PT A-IIT Simulation: " << s+startsim << " Starting burn-in period "<< std::endl;
    int track_burn_in=0;
    while(track_burn_in<burn_in){
      for(int replica=0;replica<T;replica++){//For loop for replica update in the burn-in
        int samples_replica=0;
        while(samples_replica<sample_inter_swap && samples_replica<burn_in){//Loop to create samples for each replica until we reach the defined threshold
          current_temp=temp(replica);// Extract temperature of the replica
          current_log_bound=log_bound_vector(replica);// Extract log-bound of the corresponding temperature
          current_state=X.col(replica);
          output=a_IIT_update(current_state,Q_matrix,bal_function[index_process(replica)],current_temp,current_log_bound,update_constant,0,0,max_log_bound_vector(replica));
          bool update_state=true;
          //During burn-in:
          ////// Update = true, we always update the constant
          ////// prob_to_dec=0, we never decrease the constant 
          ////// decreasing constant=0, we decrease by 0 (redundancy)
          ////// we keep track of the max log bound found
          //// Compute weight
          Z = output(1); //Extract the Z-factor
          new_samples=1+R::rgeom(Z);
          if(new_samples<1){
            Rcpp::Rcout <<"Error: geometric in "<< "simulation: " << s+startsim << " Burn-in period after " << track_burn_in <<"simulations,  temp:"<<current_temp<< std::endl;
            Rcpp::Rcout <<"new_samples= "<<new_samples<< ", Z=" << Z << " log-bound= " << current_log_bound << std::endl;
            //Show the current state where the low Z factor was identified
            double cuenta_unos = sum(current_state);
            uvec coord_print;
            int coord_shown;
            if(cuenta_unos>400){//If there are more 1s
              coord_print=find(current_state==0);
              coord_shown=0;
            }else{//If there are more 0s
              coord_print=find(current_state==1);
              coord_shown=1;
            }
            Rcpp::Rcout <<"Coords with  "<<coord_shown<< " are:\n" << coord_print << std::endl;
            
            new_samples=sample_inter_swap;
          }
          if((samples_replica+new_samples)>sample_inter_swap){//If we're going to surpass the required number of samples
            new_samples = sample_inter_swap-samples_replica;//We force to stop at sample_inter_swap
            bool update_state=false;
          }
          samples_replica+=new_samples; // Update number of samples obtained from the replica
          if(update_state){
            X.col(replica)=vec(output(0)); //Update current state of the chain
          }
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
      // Rcpp::Rcout << "PT A-IIT Simulation: " << s+startsim << ". Done " << track_burn_in <<" samples in burn-in period"<< std::endl;
    }
    Rcpp::Rcout <<"END of burn-in period\n log_bound_vector:\n "<< log_bound_vector << std::endl;
//////////////////////Finish the loop for burn-in period
    swap_count=0; //Reset swap count
    std::clock_t start = std::clock(); // Start timer for simulation s
    //// Start the loop for all iterations in simulation s
    for(int i=0;i<total_swaps;i++){
      // Rcpp::Rcout <<"Inside iteration loop"<< i << std::endl;
      temporal_mat1.fill(-1); //Reset the temporal matrix to store distance to modes
      temporal_mat2.fill(-1); //Reset the temporal matrix to store distance to modes
      temporal_origin.fill(-1);//Reset the temporal matrix to store distance to modes
      if (i % 100 == 1) {Rcpp::Rcout << "PT A-IIT Simulation: " << s+startsim << " Swap: " << i<<" Prob_decrease_bound: " << prob_to_dec << std::endl;
        Rcpp::Rcout <<"log_bound_vector:\n "<< log_bound_vector << std::endl;}
      // Rcpp::Rcout << "Simulation: " << s+startsim << " Iteration: " << i << std::endl;
      for(int replica=0;replica<T;replica++){//For loop for replicas
        // Rcpp::Rcout << "Sim: " << s+startsim << " Swap: " << i <<"replica: "<<replica<< std::endl;
        int samples_replica=0;
        while(samples_replica<sample_inter_swap){//Loop to create samples for each replica until we reach the defined threshold
          if(total_swaps<10){Rcpp::Rcout << "Replica: " << index_process(replica) << " Sampled: " << (static_cast<double>(samples_replica) / sample_inter_swap)<< std::endl;}
          total_iterations(i,index_process(replica),s)+=1;//increase the number of iterations
          current_temp=temp(index_process(replica));// Extract temperature of the replica
          current_log_bound=log_bound_vector(index_process(replica));// Extract log-bound of the corresponding temperature
          current_state=X.col(replica);
          output=a_IIT_update(current_state,Q_matrix,bal_function[index_process(replica)],current_temp,current_log_bound,update_constant,prob_to_dec,decreasing_constant,max_log_bound_vector(index_process(replica)));
          bool update_state=true;
//// Compute weight
          Z = output(1); //Extract the Z-factor
          new_samples=1+R::rgeom(Z);
          if(new_samples<1){
            Rcpp::Rcout <<"Error: geometric in "<< "simulation: " << s+startsim << " Swap: " << i <<" temperature:"<<current_temp<< std::endl;
            Rcpp::Rcout <<"new_samples= "<<new_samples<< ", Z=" << Z << " log-bound= " << current_log_bound << std::endl;
            //Show the current state where the low Z factor was identified
            double cuenta_unos = sum(current_state);
            uvec coord_print;
            int coord_shown;
            if(cuenta_unos>400){//If there are more 1s
              coord_print=find(current_state==0);
              coord_shown=0;
            }else{//If there are more 0s
              coord_print=find(current_state==1);
              coord_shown=1;
            }
            Rcpp::Rcout <<"Coords with  "<<coord_shown<< " are:\n" << coord_print.t() << std::endl;
            new_samples=sample_inter_swap;
          }
          if((samples_replica+new_samples)>sample_inter_swap){//If we're going to surpass the required number of samples
            new_samples = sample_inter_swap-samples_replica;//We force to stop at sample_inter_swap
            bool update_state=true;
          }
          
          //// Store weight of replica with temperature 1
          if(current_temp==1){ // For the original temperature replica
            // Rcpp::Rcout << "Starts update of visited states" << std::endl;
            vec curr_loglik_visited=loglikelihood_visited.col(s);
            found_min=curr_loglik_visited.index_min();
            temporal_loglik=loglik(X.col(replica),Q_matrix);
            if(curr_loglik_visited(found_min)<temporal_loglik){
              Rcpp::Rcout << "Found big likelihood " <<exp(temporal_loglik)<<" in index: "<<found_min<< std::endl;
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
          ///// Measure distance to modes
          double temp_dist1=L1_distance(X.col(replica),mode1);
          double temp_dist2=L1_distance(X.col(replica),mode2);
          double temp_origin=sum(X.col(replica));
          // Rcpp::Rcout <<"Before storing in temporal matrix" << std::endl;
          ///// Store the measured distance in the temporal matrix
          for(int j=samples_replica;j<(samples_replica+new_samples);j++){
            temporal_mat1(j,index_process(replica))=temp_dist1;
            temporal_mat2(j,index_process(replica))=temp_dist2;
            temporal_origin(j,index_process(replica))=temp_origin;
          }
          // Rcpp::Rcout <<"Succesfully stored in temporal matrix" << std::endl;
          ///// Updating before the next iteration of the loop        
          samples_replica+=new_samples; // Update number of samples obtained from the replica
          if(update_state){
            X.col(replica)=vec(output(0)); //Update current state of the chain
          }
          log_bound_vector(index_process(replica))=output(2); //Update log-bound
          max_log_bound_vector(index_process(replica))=output(3); //Update maximum log-bound found
        }//End loop to update a single replica
      }//End loop to update all replicas
      
///// Store the temporal matrix into the big matrix
      // Rcpp::Rcout <<"Storing distance in big matrix" << std::endl;
      distance_mode1.rows(i*sample_inter_swap,((i+1)*sample_inter_swap)-1)=temporal_mat1;
      distance_mode2.rows(i*sample_inter_swap,((i+1)*sample_inter_swap)-1)=temporal_mat2;
      distance_origin.rows(i*sample_inter_swap,((i+1)*sample_inter_swap)-1)=temporal_origin;
      // Rcpp::Rcout <<"Succesfully stored in big matrix" << std::endl;
      //// Start replica swap process
      // Rcpp::Rcout <<"log_bound_vector:\n "<< log_bound_vector << std::endl;    
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
    full_distance_mode1.slice(s)=distance_mode1;//Store distance to mode 1
    full_distance_mode2.slice(s)=distance_mode2;//Store distance to mode 2
    full_distance_origin.slice(s)=distance_origin;//Store distance to origin
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
  // ret["modes_visit"]=full_iter_visit_modes;
  ret["distance_mode1"]=full_distance_mode1;
  ret["distance_mode2"]=full_distance_mode2;
  ret["distance_origin"]=full_distance_origin;
  ret["max_bounds"]=max_log_bound_vector;
  ret["final_bounds"]=log_bound_vector;
  return ret;
}

// [[Rcpp::export]]
List PT_a_IIT_sim_RF(int p,int startsim,int endsim, int numiter, int iterswap,int burn_in, vec temp, const std::vector<std::string>& bal_function, bool bias_fix,const std::string& filename,int num_states_visited,const std::vector<int>& starting_coord, double decreasing_constant,std::string reduc_model){
  //// Initialize variables to use in the code
  int T=temp.n_rows; // Count number of temperatures
  vec log_bound_vector(T); // vector to store a log-bound for each replica
  vec max_log_bound_vector(T); // vector to store the MAX log-bound found for each replica
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
  loglikelihood_visited.fill(-10000);//Initialize a very negative loglikelihood
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
  vec mode1=Q_matrix.col(0);
  vec mode2=Q_matrix.col(1);
  mat distance_mode1(numiter,T);//Matrix to store distance to mode 1 for each replica
  mat distance_mode2(numiter,T);//Matrix to store distance to mode 1 for each replica
  mat distance_origin(numiter,T);//Matrix to store distance to mode 1 for each replica
  cube full_distance_mode1(numiter,T,total_sim);//Cube to store distance to mode 1 for each replica and simulation
  cube full_distance_mode2(numiter,T,total_sim);//Cube to store distance to mode 2 for each replica and simulation
  cube full_distance_origin(numiter,T,total_sim);//Cube to store distance to mode 2 for each replica and simulation

  // Rcpp::Rcout << "First rows Q_matrix: " << Q_matrix.rows(0,5) << std::endl;
  // Rcpp::Rcout << "Last rows Q_matrix: " << Q_matrix.rows(p-6,p-1) << std::endl;
  std::vector<double> time_taken(total_sim); // vector to store the seconds each process took
  // Probability to update
  bool update_prob=false;
  bool update_constant=true;
  double prob_to_dec=0;
  double percentage_start=0.05;
  double percentage_end=0.70;
  int total_replica_iterations=numiter;
  int sample_iterations_count;
  if(false){}
  else if(reduc_model=="always"){prob_to_dec=1;}
  else if(reduc_model=="never"){prob_to_dec=0;}
  else if(reduc_model=="iterations"){update_prob=true;}
  else if(reduc_model=="zero"){update_constant=false;}
  else {Rcpp::Rcout <<"reduc_model= " <<reduc_model<<" is not a valid reduc_model. Default to standard"<< std::endl;
    Rcpp::Rcout <<" The standard is: Bounding constant always increases."<< std::endl;
    prob_to_dec=0;}//If we don't define a reduc_model
  //// Start the loop for all simulations
  for(int s=0;s<total_sim;s++){
    for(int i=0;i<T;i++){ // Reset index process vector at the start of each simulation
      index_process.row(i)=i;
    }
    ind_pro_hist.row(0)=index_process.t(); // First entry of the index process
    swap_count=0; //Reset swap count
    // X=initializeMatrix(starting_coord,p,T);//Reset the starting point of all chains
    ppp=Rcpp::runif(1);
    X=initializeRandom(p,T,ppp(0));//Randomly initialize the state of each replica.
    
    swap_total.zeros();
    swap_success.zeros();
    log_bound_vector.zeros();//Reset log-bounds, all log-bounds start at 0
    distance_mode1.fill(-1);
    distance_mode2.fill(-1);
    distance_origin.fill(-1);
    
    //Reset the probability to reduce the bounding constant
    if(reduc_model=="iterations"){update_prob=true;prob_to_dec=1;} //Reset the bool to update probability
    sample_iterations_count=0; // Reset the counting of iterations (or samples)
    
    //// Start loop for burn_in period
    for(int i=0;i<burn_in;i++){
      if (i % 100 == 1) {Rcpp::Rcout << "Simulation: " << s+startsim << " Burn_in period, iteration: " << i << std::endl;}
      for(int replica=0;replica<T;replica++){//For loop for replica update
        current_temp=temp(index_process(replica));
        current_log_bound=log_bound_vector(replica);// Extract log-bound of the corresponding temperature
        
        output=a_IIT_update(X.col(replica),Q_matrix,bal_function[index_process(replica)],current_temp,current_log_bound,update_constant,0,0,max_log_bound_vector(replica));
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
            // For replica swaps we don't update the bounding constant
            output=a_IIT_update(Xtemp_from,Q_matrix,bal_function[t],temp(t),log_bound_vector(t),false,0,0,max_log_bound_vector(t));
            Z_temp11=output(1);
            output=a_IIT_update(Xtemp_to,Q_matrix,bal_function[t+1],temp(t+1),log_bound_vector(t+1),false,0,0,max_log_bound_vector(t+1));
            Z_temp22=output(1);
            output=a_IIT_update(Xtemp_from,Q_matrix,bal_function[t+1],temp(t+1),log_bound_vector(t+1),false,0,0,max_log_bound_vector(t+1));
            Z_temp12=output(1);
            output=a_IIT_update(Xtemp_to,Q_matrix,bal_function[t],temp(t),log_bound_vector(t),false,0,0,max_log_bound_vector(t));
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
        //Measure distance to the modes
        distance_mode1(i,index_process(replica))=L1_distance(X.col(replica),mode1);
        distance_mode2(i,index_process(replica))=L1_distance(X.col(replica),mode2);
        distance_origin(i,index_process(replica))=sum(X.col(replica));
        //Depending on the chosen method
        //// Update each replica independently
        // output=IIT_update_w(X.col(replica),Q_matrix,bal_function[index_process(replica)],current_temp);
        output=a_IIT_update(X.col(replica),Q_matrix,bal_function[index_process(replica)],current_temp,current_log_bound,update_constant,prob_to_dec,decreasing_constant,max_log_bound_vector(index_process(replica)));
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
                  Rcpp::Rcout <<"New prob: "<< prob_to_dec << std::endl;
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
            output=a_IIT_update(Xtemp_from,Q_matrix,bal_function[t],temp(t),log_bound_vector(t),false,0,0,max_log_bound_vector(t));
            Z_temp11=output(1);
            output=a_IIT_update(Xtemp_to,Q_matrix,bal_function[t+1],temp(t+1),log_bound_vector(t+1),false,0,0,max_log_bound_vector(t+1));
            Z_temp22=output(1);
            output=a_IIT_update(Xtemp_from,Q_matrix,bal_function[t+1],temp(t+1),log_bound_vector(t+1),false,0,0,max_log_bound_vector(t+1));
            Z_temp12=output(1);
            output=a_IIT_update(Xtemp_to,Q_matrix,bal_function[t],temp(t),log_bound_vector(t),false,0,0,max_log_bound_vector(t));
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
    full_distance_mode1.slice(s)=distance_mode1;
    full_distance_mode2.slice(s)=distance_mode2;
    full_distance_origin.slice(s)=distance_origin;
    // Rcpp::Rcout <<"Final state "<< X << std::endl;
  }//End loop simulations
  List ret;
  
  ret["ip"]=ind_pro_hist;
  ret["swap_rate"]=swap_rate;
  ret["states"]=states_visited;
  ret["loglik_visited"]=loglikelihood_visited;
  ret["iter_visit"]=iter_to_visit;
  ret["time_taken"]=time_taken;
  // ret["modes_visit"]=full_iter_visit_modes;
  ret["distance_mode1"]=full_distance_mode1;
  ret["distance_mode2"]=full_distance_mode2;
  ret["distance_origin"]=full_distance_origin;
  ret["max_bounds"]=max_log_bound_vector;
  ret["final_bounds"]=log_bound_vector;
  return ret;
}


