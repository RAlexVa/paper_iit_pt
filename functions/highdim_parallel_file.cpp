// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;
using namespace arma;

#include <ctime> 
#include "workers.h"


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
arma::Mat<double> initializeRandom(const int num_rows,const int num_cols, const double prob) {
  
  // Initialize a matrix with random values between 0 and 1
  arma::mat A(num_rows, num_cols,fill::randu);
  
  // Threshold the random values to 0 or 1
  A = arma::conv_to<arma::Mat<double>>::from(A > prob);
  
  return(A);
}

////////// Creating sparse matrix from file //////////
// [[Rcpp::export]]
arma::mat readMatrixFile(const std::string& filename){
  // Open the file
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open file: " + filename);
  }
  
  // Read the first line to get the matrix dimensions and number of non-zero entries
  size_t n_rows, n_nonzeros;
  file >> n_rows >> n_nonzeros;
  
  // Create a sparse matrix
  arma::mat matrix(n_rows, n_rows);
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
int problem_dim(const std::string& filename){
  // Open the file
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open file: " + filename);
  }
  
  // Read the first line to get the matrix dimensions and number of non-zero entries
  size_t n_rows, n_nonzeros;
  file >> n_rows >> n_nonzeros;
  
  return n_rows;
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
double bf_sq(double x){
  return x/2;
}
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
double bound_sq(double x, double log_gamma){
  double temp1=x/2-log_gamma;
  return ret_min(temp1,x,0);
}
///// Functions to call balancing functions inside other functions
//This function we don't export to R because it generates errors
double invoke(double x, double (*func)(double)) {
  return func(x);
}
// double bal_func(double x,String chosen){
//   if (chosen == "sq") {
//     return invoke(x, &bf_sq);
//   } else if (chosen == "min") {
//     return invoke(x, &bf_min);
//   } else {
//     Rcpp::Rcout <<"Name of the balancing function does exist!" << std::endl;
//     return 0; // Default return for unknown operation
//   }
// }

double apply_bal_func(double x,const int chosen){
  if (chosen == 2) {
    return invoke(x, &bf_sq);
  } else if (chosen == 1) {
    return invoke(x, &bf_min);
  } else {
    Rcpp::Rcout <<"The balancing function does exist!" << std::endl;
    return -1000; // Default return for unknown operation
  }
}


// [[Rcpp::export]]
double loglik(const arma::vec& X,const arma::mat&  M){
  double temporal=arma::as_scalar(X.t() * M * X);
  if(temporal==0){return(-0.05);//This is just for the vector 0
    }else{return(log(temporal));}
}



// [[Rcpp::export]]
double loglik_R(Rcpp::NumericVector& X, NumericMatrix& M) {
  // Convert directly using RcppArmadillo
  arma::vec VEC = Rcpp::as<arma::vec>(X);
  arma::mat MAT = Rcpp::as<arma::mat>(M);
  return loglik(VEC, MAT);
}



///// Full definition of internal functions of workers 
double IIT_visit_neighbors::apply_bal_func_internal(double x,const int chosen){
  return apply_bal_func(x,chosen);
}

double IIT_visit_neighbors::loglik_internal(const arma::Col<double>& X,const arma::Mat<double>& M, const double& theta){
  // return loglik(X,M,theta);
  return loglik(X,M);
}

double IIT_visit_bounded::apply_bal_func_bounded_internal(double x,double log_bound, int bal_func){
  if(bal_func==1){
    return bf_min(x);//Apply min balancing function ignoring the log-bound
  }else   if(bal_func==2){
    return bound_sq(x,log_bound); 
  }else{
    Rcpp::Rcout <<" The balancing function is incorrect (bal_func_bounded_internal)"<< std::endl;
    
  }
  
}

double IIT_visit_bounded::loglik_internal(const arma::Col<double>& X,const arma::Mat<double>& M, const double& theta){
  // return loglik(X,M,theta);
  return loglik(X,M);
}


///// Update function that is not rejection free
// [[Rcpp::export]]
List single_step_update(NumericVector currentX, NumericMatrix Q,int p, int bal_func, double current_temp, double theta, double current_log_bound){
  NumericVector newX=clone(currentX);
  
  double current_loglik=loglik_R(currentX,Q);
  double rand_val = arma::randu(); // Random value in [0, 1)
  int random_neighbor = static_cast<int>(rand_val * p); // Scale to [0, p) and cast to int
  
  newX(random_neighbor) = 1-newX(random_neighbor); // Update coordinate
  double new_loglik=loglik_R(newX,Q);
  
  double logratio_probs = current_temp*(new_loglik - current_loglik);
  
  bool success_jump=false;
  
  // Apply the corresponding balancing function
  double jump_logprob=0;
  if(bal_func==1){//Using MIN bal_fun (Traditional Metropolis algorithm)
    jump_logprob=bf_min(logratio_probs);
  }else if(bal_func==2){//Using bounded Sqrt root 
    jump_logprob=bound_sq(logratio_probs, current_log_bound);
  }else{
    Rcpp::Rcout <<"The chosen balancing function is not correct"<< std::endl;
    jump_logprob=-10000;
  }
  double rand_jump = arma::randu();
  if(log(rand_jump)<jump_logprob){//We jumped to the new state
    success_jump=true;
  }else{//If we reject the jump
    // newX = clone(currentX);//We stay in the same state
  }
  
  List ret_single_step;
  ret_single_step["jump"]=success_jump;
  ret_single_step["coord"]=random_neighbor;
  
  return ret_single_step;
  
}

////////// Code for Parallel Tempering simulations //////////
// [[Rcpp::export]]
List PT_IIT_sim(int p,int startsim,int endsim, int numiter, int iterswap,int burn_in, vec temp, int bal_func, bool bias_fix,const std::string& filename,int num_states_visited,const std::vector<int>& starting_coord, double theta, int num_modes){
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
  // mat X(p,T); // To store the current state of the joint chain, as many rows as neighbors, as many columns as temperatures
  NumericMatrix X(p,T);
  NumericMatrix initialX(p,T);
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
  NumericVector Xtemp_from(p);//vec Xtemp_from(p);
  NumericVector Xtemp_to(p);//vec Xtemp_to(p);
  // Variables to register visits to high probability states
  // Keep track of the state with the highest likelihood for each replica
  mat iter_to_visit(T,total_sim);
  mat loglikelihood_visited(T,total_sim);
  mat time_visited(T,total_sim);
  
  loglikelihood_visited.fill(-10000);//Initialize a very negative loglikelihood
  time_visited.fill(-1);
  cube states_visited(p,T,total_sim,fill::zeros);
  Rcpp::Rcout <<"p= "<<p<<", temps= "<<T<<", total_sim= "<<total_sim<< std::endl;
  // mat modes_visited(numiter * total_sim,T);//Matrix to store the modes visited and temperature
  const std::size_t dim_size = static_cast <size_t> (p); 
  const std::size_t number_modes = static_cast <size_t> (num_modes); 
  //// Define matrix to use for the likelihood

  arma::mat Q_matrix = readMatrixFile(filename);
  NumericMatrix Q_mat_R=Rcpp::wrap(Q_matrix);//Numeric Matrix version of the Q_matrix

  bool finish_sim;
  
  std::vector<double> time_taken(total_sim); // vector to store the seconds each process took
  //// Start the loop for all simulations
  for(int s=0;s<total_sim;s++){
    for(int i=0;i<T;i++){ // Reset index process vector at the start of each simulation
      index_process.row(i)=i;
    }
    ind_pro_hist.row(0)=index_process.t(); // First entry of the index process
    swap_count=0; //Reset swap count
    // X=initializeMatrix(starting_coord,p,T);//Reset the starting point of all chains
    // double ppm=randu();
    arma::Mat<double> inter_mat(p,T);
    // inter_mat=initializeRandom(p,T,ppm);//Randomly initialize the state of each replica.
    inter_mat=initializeRandom(p,T,0.5);
    X=Rcpp::wrap(inter_mat);
    initialX=clone(X);
    swap_total.zeros();
    swap_success.zeros();
    finish_sim=false;
    
    std::clock_t start = std::clock(); // Start timer for simulation s
/////// Start loop for burn_in period
    for(int i=0;i<burn_in;i++){
      if (i % 10000 == 1) {Rcpp::Rcout << "PT-IIT Simulation: " << s+startsim << " Burn_in period, iteration: " << i << std::endl;}
      for(int replica=0;replica<T;replica++){//For loop for replica update
        current_temp=temp(index_process(replica));
        
        NumericVector output_current_X(p);
        NumericVector current_X=X(_,replica);
        
        IIT_visit_neighbors visit_current_X(current_X,
                                            Q_mat_R,
                                            bal_func,
                                            current_temp,
                                            theta,
                                            output_current_X,
                                            dim_size,
                                            dim_size);
        parallelFor(0,dim_size,visit_current_X);//Apply ParallelFor
        
        ////Sample Proportionally
        //Get random uniforms
        arma::Col<double> u_random(p,fill::randu);
        NumericVector u=Rcpp::wrap(u_random);
        //Compute the needed values
        NumericVector choose_min=log(-log(u)) - (output_current_X);
        //Find the index of the minimum entry
        GetMin min_coord(choose_min);
        parallelReduce(0,dim_size,min_coord);
        //Swap that coordinate
        X(min_coord.min_index,replica)=1-X(min_coord.min_index,replica);
        // X.col(replica)=vec(output(0)); //Update current state of the chain
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
          // 
          Xtemp_from=X(_,t);
          Xtemp_to=X(_,t+1);
          //// Optional Computation of Z factors to correct bias
          double Z_fact_correc=1;//to temporarily store Z_factor correction
          if(bias_fix){
            //Declare vector to store info of visiting neighbors
            NumericVector output_Z_v11(p); 
            NumericVector output_Z_v22(p); 
            NumericVector output_Z_v12(p); 
            NumericVector output_Z_v21(p); 
            //// Declare constructor to visit all neighbors
            IIT_visit_neighbors Z_v11(Xtemp_from,
                                      Q_mat_R,
                                      bal_func,
                                      temp(t),
                                      theta,
                                      output_Z_v11,
                                      dim_size,
                                      dim_size);
            IIT_visit_neighbors Z_v22(Xtemp_to,
                                      Q_mat_R,
                                      bal_func,
                                      temp(t+1),
                                      theta,
                                      output_Z_v22,
                                      dim_size,
                                      dim_size);
            IIT_visit_neighbors Z_v12(Xtemp_from,
                                      Q_mat_R,
                                      bal_func,
                                      temp(t+1),
                                      theta,
                                      output_Z_v12,
                                      dim_size,
                                      dim_size);
            IIT_visit_neighbors Z_v21(Xtemp_to,
                                      Q_mat_R,
                                      bal_func,
                                      temp(t),
                                      theta,
                                      output_Z_v21,
                                      dim_size,
                                      dim_size);
            //// Apply ParallelFor
            parallelFor(0,dim_size,Z_v11);
            parallelFor(0,dim_size,Z_v22);
            parallelFor(0,dim_size,Z_v12);
            parallelFor(0,dim_size,Z_v21);
            //// Declare constructor to add log-probabilities
            SumExp sum_Z11(output_Z_v11);
            SumExp sum_Z22(output_Z_v22);
            SumExp sum_Z12(output_Z_v12);
            SumExp sum_Z21(output_Z_v21);
            //// Get the sum of probabilities
            parallelReduce(0,dim_size,sum_Z11);
            parallelReduce(0,dim_size,sum_Z22);
            parallelReduce(0,dim_size,sum_Z12);
            parallelReduce(0,dim_size,sum_Z21);
            //Compute correction factor
            Z_fact_correc=(sum_Z12.Z*sum_Z21.Z)/(sum_Z11.Z*sum_Z22.Z);
          }//Finish IF bias_fix
          //// Computing swap probability
          swap_prob=(temp(t)-temp(t+1))*(loglik_R(Xtemp_to,Q_mat_R) - loglik_R(Xtemp_from,Q_mat_R)); 
          swap_prob=Z_fact_correc*exp(swap_prob);
          // Rcpp::Rcout <<"Swap prob "<< swap_prob << std::endl;
          ppp=Rcpp::runif(1);
          if(ppp(0)<swap_prob){//In case the swap is accepted
            //Swap vectors in the matrix X
            for(int coord=0;coord<p;coord++){
              X(coord,t+1)=Xtemp_from[coord];
              X(coord,t)=Xtemp_to[coord];
            }
            
          }//In case a swap is accepted
        }//Finish loop to swap many replicas
      }//Finish IF for replica swap
    }// Finish burn-in period
    
/////////////////////// Finish burn-in period
    swap_count=0; //Reset swap count
    
    
////// Start the loop for all iterations in simulation s
    for(int i=0;i<numiter;i++){
      if (i % 10000 == 1) {Rcpp::Rcout << "PT-IIT Simulation: " << s+startsim << " Iteration: " << i << std::endl;}
      // Rcpp::Rcout << "State X:\n " <<X << std::endl;
      for(int replica=0;replica<T;replica++){//For loop for replicas
        int temperature_index=index_process(replica);
        current_temp=temp(temperature_index);// Extract temperature of the replica
        NumericVector output_current_X(p);
        NumericVector current_X=X(_,replica);
        IIT_visit_neighbors visit_current_X(current_X,
                                            Q_mat_R,
                                            bal_func,
                                            current_temp,
                                            theta,
                                            output_current_X,
                                            dim_size,
                                            dim_size);
        parallelFor(0,dim_size,visit_current_X);//Apply ParallelFor
        
        ////Sample Proportionally
        //Get random uniforms
        arma::Col<double> u_random(p,fill::randu);
        NumericVector u=Rcpp::wrap(u_random);
        //Compute the needed values
        NumericVector choose_min=log(-log(u)) - (output_current_X);
        //Find the index of the minimum entry
        GetMin min_coord(choose_min);
        parallelReduce(0,dim_size,min_coord);
// Record the likelihood

double current_loglik=loglik_R(current_X,Q_mat_R);

if(current_loglik>loglikelihood_visited(temperature_index,s)){//If the current loglik is bigger
  //Store time
  std::clock_t time_find_mode = std::clock();
  double secs_find_mode = static_cast<double>(time_find_mode - start) / CLOCKS_PER_SEC;
  time_visited(temperature_index,s)=secs_find_mode;
  //Store iteration
  iter_to_visit(temperature_index,s)=i;
  //Store value of loglik
  loglikelihood_visited(temperature_index,s)=current_loglik;
  //Store state
  arma::vec current_X_arma = Rcpp::as<arma::vec>(current_X);
  states_visited.slice(s).col(temperature_index) = current_X_arma;
  if(current_loglik>log(11500)){
    Rcpp::Rcout <<"Found higher lik: "<<exp(current_loglik)<<" in iteration: "<<i<<" in temp: "<<temperature_index<<" in time: "<<secs_find_mode<< std::endl; 
  }

}

        //Swap that coordinate
        X(min_coord.min_index,replica)=1-X(min_coord.min_index,replica);
        
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
          // Xtemp_from=Rcpp::wrap(X.cols(find(index_process==t)));
          // Xtemp_to=Rcpp::wrap(X.cols(find(index_process==t+1)));
          auto find_t = std::find(index_process.begin(), index_process.end(), t);
          auto find_t_1 = std::find(index_process.begin(), index_process.end(), t+1);
          size_t index_t=std::distance(index_process.begin(), find_t);
          size_t index_t_1=std::distance(index_process.begin(), find_t_1);
          Xtemp_from=X(_,index_t);
          Xtemp_to=X(_,index_t_1);
          
          
          //// Optional Computation of Z factors to correct bias
          double Z_fact_correc=1;//to temporarily store Z_factor correction
          if(bias_fix){
            //Declare vector to store info of visiting neighbors
            NumericVector output_Z_v11(p); 
            NumericVector output_Z_v22(p); 
            NumericVector output_Z_v12(p); 
            NumericVector output_Z_v21(p); 
            //// Declare constructor to visit all neighbors
            IIT_visit_neighbors Z_v11(Xtemp_from,
                                      Q_mat_R,
                                      bal_func,
                                      temp(t),
                                      theta,
                                      output_Z_v11,
                                      dim_size,
                                      dim_size);
            IIT_visit_neighbors Z_v22(Xtemp_to,
                                      Q_mat_R,
                                      bal_func,
                                      temp(t+1),
                                      theta,
                                      output_Z_v22,
                                      dim_size,
                                      dim_size);
            IIT_visit_neighbors Z_v12(Xtemp_from,
                                      Q_mat_R,
                                      bal_func,
                                      temp(t+1),
                                      theta,
                                      output_Z_v12,
                                      dim_size,
                                      dim_size);
            IIT_visit_neighbors Z_v21(Xtemp_to,
                                      Q_mat_R,
                                      bal_func,
                                      temp(t),
                                      theta,
                                      output_Z_v21,
                                      dim_size,
                                      dim_size);
            //// Apply ParallelFor
            parallelFor(0,dim_size,Z_v11);
            parallelFor(0,dim_size,Z_v22);
            parallelFor(0,dim_size,Z_v12);
            parallelFor(0,dim_size,Z_v21);
            
            //// Declare constructor to add log-probabilities
            SumExp sum_Z11(output_Z_v11);
            SumExp sum_Z22(output_Z_v22);
            SumExp sum_Z12(output_Z_v12);
            SumExp sum_Z21(output_Z_v21);
            //// Get the sum of probabilities
            parallelReduce(0,dim_size,sum_Z11);
            parallelReduce(0,dim_size,sum_Z22);
            parallelReduce(0,dim_size,sum_Z12);
            parallelReduce(0,dim_size,sum_Z21);
            // Rcpp::Rcout <<"Finish parallel steps swap: "<<swap_count<< std::endl; 
            Z_fact_correc=(sum_Z12.Z*sum_Z21.Z)/(sum_Z11.Z*sum_Z22.Z);
            // Rcpp::Rcout <<"Correction factors: "<< sum_Z12.Z <<","<< sum_Z21.Z<<","<<sum_Z11.Z <<","<<sum_Z22.Z << std::endl;
            
          }
          //// Computing swap probability
          swap_prob=(temp(t)-temp(t+1))*(loglik_R(Xtemp_to,Q_mat_R) - loglik_R(Xtemp_from,Q_mat_R)); 
          // Rcpp::Rcout <<"Correction factor: "<< Z_fact_correc << std::endl;
          // Rcpp::Rcout <<"Swap prob: "<< swap_prob << std::endl;
          swap_prob=Z_fact_correc*exp(swap_prob);
          ppp=Rcpp::runif(1);
          
          if(ppp(0)<swap_prob){//In case the swap is accepted
            swap_success(t)+=1;//Increase the number of successful swaps of temp t
            // Rcpp::Rcout <<"Accepted swap: " <<swap_count<< std::endl;
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
    std::clock_t end_time = std::clock(); // Stop timer
    // Calculate the time taken in seconds
    double duration = static_cast<double>(end_time - start) / CLOCKS_PER_SEC;
    time_taken[s] = duration;
    // Store result of the simulation
    vec temp_rate=swap_success / swap_total;
    swap_rate.row(s)=temp_rate.t();
  }//End loop simulations
  
// Print highest likelihood visited  
  Rcpp::Rcout <<"Found liks: \n"<<exp(loglikelihood_visited )<< std::endl; 
  
  
  List ret;
  ret["ip"]=ind_pro_hist;
  ret["swap_rate"]=swap_rate;
  ret["loglik_visit"]=loglikelihood_visited;
  ret["iter_visit"]=iter_to_visit;
  ret["time_visit"]=time_visited;
  ret["state_visit"]=states_visited;
  ret["time_taken"]=time_taken;
  ret["initial_X"]=initialX;

  return ret;
}

// [[Rcpp::export]]
List PT_a_IIT_sim(int p,int startsim,int endsim, int total_swaps,int sample_inter_swap,int burn_in, vec temp, const int bal_func,const std::string& filename,int num_states_visited,const std::vector<int>& starting_coord, double decreasing_constant,std::string reduc_model, double theta, int num_modes, int temps_rf){
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
  NumericMatrix X(p,T); // To store the current state of the joint chain, as many rows as neighbors, as many columns as temperatures
  vec index_process(T);   //Initialize index process vector
  mat ind_pro_hist(total_swaps*total_sim+1,T); //To store evolution of index process
  // int max_num=pow(2,p);
  // Rcpp::Rcout << "max num: " << max_num << std::endl;
  // vec current_state(p);//To print when I find a very small Z factor
  NumericVector current_X(p);
  
  vec swap_total(J,fill::ones);
  swap_total*=total_swaps;//We always have the same number of total swaps
  int final_swap=total_swaps; //Last swap before breaking the for loop
  //By default the final swap will be the total number of specified swaps
  //If it finds the modes earlier, then it's updated before used to compute swap rate
  vec swap_success(J);
  mat swap_rate(total_sim,J);
  cube total_iterations(total_swaps,T,total_sim,fill::zeros);//To store iterations needed in between swaps
  ////Variables to update index process
  vec epsilon_indic(T); //Vector that indicates if the entry of the index process is proposed to change
  vec prop_swap(T); //vector to indicate a proposed swap
  vec do_swap(T); //vector to indicate actually perform a swap
  vec resulting_swap(T);//Vector to update the index process
  double ppp; //To store the probability of replica swap
  //// Variables to perform the replica swap
  double swap_prob;
  NumericVector Xtemp_from(p);
  NumericVector Xtemp_to(p);
  bool update_state=true;
  // Variables to register visits to high probability states
  // Keep track of the state with the highest likelihood for each replica
  mat iter_to_visit(T,total_sim);
  mat loglikelihood_visited(T,total_sim);
  mat time_visited(T,total_sim);
  
  loglikelihood_visited.fill(-10000);//Initialize a very negative loglikelihood
  time_visited.fill(-1);
  cube states_visited(p,T,total_sim,fill::zeros);
  Rcpp::Rcout <<"p= "<<p<<", temps= "<<T<<", total_sim= "<<total_sim<< std::endl;
  
  arma::mat Q_matrix = readMatrixFile(filename);

  NumericMatrix Q_mat_R=Rcpp::wrap(Q_matrix);//Numeric Matrix version of the Q_matrix
  
  const std::size_t dim_size = static_cast <size_t> (p);
  const std::size_t numer_modes = static_cast <size_t> (num_modes); 
  
  
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
  
  List single_step_output;
  
  //// Start the loop for all simulations
  for(int s=0;s<total_sim;s++){
    for(int i=0;i<T;i++){ // Reset index process vector at the start of each simulation
      index_process.row(i)=i;
    }
    ind_pro_hist.row(0)=index_process.t(); // First entry of the index process
    swap_count=0; //Reset swap count
    
    // double ppm=randu();
    arma::Mat<double> inter_mat(p,T);
    // inter_mat=initializeRandom(p,T,ppm);//Randomly initialize the state of each replica.
    inter_mat=initializeRandom(p,T,0.5);
    X=Rcpp::wrap(inter_mat);
    
    log_bound_vector.zeros();//Reset log-bounds, all log-bounds start at 0
    swap_success.zeros();
    //Reset the probability to reduce the bounding constant
    if(reduc_model=="iterations"){update_prob=true;prob_to_dec=1;} //Reset the bool to update probability
    sample_iterations_count=0; // Reset the counting of iterations (or samples)
    ////Start the loop for burn-in period
    Rcpp::Rcout << "PT A-IIT Simulation: " << s+startsim << " Starting burn-in period "<< std::endl;
    int track_burn_in=0;
    // Start tracking of time
    std::clock_t start = std::clock(); // Start timer for simulation s
    while(track_burn_in<burn_in){
      for(int replica=0;replica<T;replica++){//For loop for replica update in the burn-in
        int samples_replica=0;
        while(samples_replica<sample_inter_swap && samples_replica<burn_in){//Loop to create samples for each replica until we reach the defined threshold
          update_state=true;
          current_temp=temp(replica);// Extract temperature of the replica
          current_X=X(_,replica);
          NumericVector output_current_X(p);
          //Visit neighbors of current state
          // Rcpp::Rcout <<"Declaring vising IIT neighbors"<< std::endl;
          IIT_visit_neighbors visit_current_X(current_X,
                                              Q_mat_R,
                                              bal_func,
                                              current_temp,
                                              theta,
                                              output_current_X,
                                              dim_size,
                                              dim_size);
          // Rcpp::Rcout <<"Before parallelFor IIT neighbors"<< std::endl;
          
          parallelFor(0,dim_size,visit_current_X);//Apply ParallelFor
          // Rcpp::Rcout <<"After parallelFor IIT neighbors"<< std::endl;
          
          // Update the bounding constant
          // Rcpp::Rcout <<"Declaring getmax"<< std::endl;
          GetMax get_max(output_current_X);
          // Rcpp::Rcout <<"Before parallelReduce GetMax"<< std::endl;
          parallelReduce(0,dim_size,get_max);
          // Rcpp::Rcout <<"After parallelReduce GetMax"<< std::endl;
          //Always increase the bounding constant
          log_bound_vector(replica)=ret_max(get_max.max_value,log_bound_vector(replica),0);
          
          current_log_bound=log_bound_vector(replica);
          if(current_log_bound>700){
            Rcpp::Rcout <<"Replica:"<<replica<<" Current log-bound:"<<current_log_bound<< std::endl;
            Rcpp::Rcout <<"Current X= \n"<<current_X<< std::endl;
          }
          NumericVector bounded_vector=output_current_X - current_log_bound;
          // Rcpp::Rcout <<"Declaring getmax"<< std::endl;
          SumExp get_sum(bounded_vector);
          // Rcpp::Rcout <<"Before parallelReduce SumExp"<< std::endl;
          parallelReduce(0,dim_size,get_sum);
          // Rcpp::Rcout <<"After parallelReduce SumExp"<< std::endl;
          Z=get_sum.Z/p;//Divide over the number of neighbors since we're using uniform distribution
          new_samples=1+R::rgeom(Z);//Get multiplicity list
          if(new_samples<1){
            Rcpp::Rcout <<"Error: geometric in "<< "simulation: " << s+startsim << " Burn-in period after " << track_burn_in <<"simulations,  temp:"<<current_temp<< std::endl;
            Rcpp::Rcout <<"new_samples= "<<new_samples<< ", Z=" << Z << " log-bound= " << current_log_bound << std::endl;
            Rcpp::Rcout <<"Current X= \n"<<current_X<< std::endl;
            new_samples=sample_inter_swap;
          }
          
          if((samples_replica+new_samples)>sample_inter_swap){//If we're going to surpass the required number of samples
            new_samples = sample_inter_swap-samples_replica;//We force to stop at sample_inter_swap
            update_state=false;
          }
          samples_replica+=new_samples; // Update number of samples obtained from the replica
          if(update_state){
            ////Sample Proportionally
            //Get random uniforms
            arma::Col<double> u_random(p,fill::randu);
            NumericVector u=Rcpp::wrap(u_random);
            //Compute the needed values
            NumericVector choose_min=log(-log(u)) - (output_current_X);//We can sample from
            //Find the index of the minimum entry
            GetMin min_coord(choose_min);
            parallelReduce(0,dim_size,min_coord);
            //Swap that coordinate
            X(min_coord.min_index,replica)=1-X(min_coord.min_index,replica);
          }//End If for updating state
          
        }//End while loop to update replicas
      }//End loop to update replicas in the burn-in
      
      //// Start replica swap process for the burn-in
      
      swap_count+=1;//Increase the count of swaps
      //Try a replica swap after reaching sample_inter_swap in each replica
      //We're doing non-reversible parallel tempering
      int starting=swap_count%2; // Detect if it's even or odd
      // Rcpp::Rcout <<"Trying replica swap "<<swap_count<<" start: "<<starting <<" at iteration: "<< i << std::endl;
      for(int t=starting;t<J;t+=2){// For loop that runs over temperature indexes to swap
        Xtemp_from=X(_,t);
        Xtemp_to=X(_,t+1);
        
        //// Computing swap probability
        swap_prob=(temp(t)-temp(t+1))*(loglik_R(Xtemp_to,Q_mat_R) - loglik_R(Xtemp_from,Q_mat_R));
        swap_prob=exp(swap_prob);
        // Rcpp::Rcout <<"Swap prob "<< swap_prob << std::endl;
        ppp=randu();
        if(ppp<swap_prob){//In case the swap is accepted
          //Swap vectors in the matrix X
          for(int coord=0;coord<p;coord++){
            X(coord,t+1)=Xtemp_from[coord];
            X(coord,t)=Xtemp_to[coord];
          }
        }
      }//End for loop for swaping odd-even replicas
      track_burn_in+=sample_inter_swap;
      
    }//End while loop to track burn-in
    Rcpp::Rcout <<"END of burn-in period\n log_bound_vector:\n "<< log_bound_vector << std::endl;
    //////////////////////Finish the loop for burn-in period
    max_log_bound_vector=log_bound_vector;
    swap_count=0; //Reset swap count
    
    
    //// Start the loop for all iterations in simulation s
    for(int i=0;i<total_swaps;i++){
      if (i % 1000 == 1) {Rcpp::Rcout << "PT A-IIT Simulation: " << s+startsim << " Swap: " << i<<" Prob_decrease_bound: " << prob_to_dec << std::endl;}
      // Rcpp::Rcout <<"log_bound_vector:\n "<< log_bound_vector << std::endl;}
      for(int replica=0;replica<T;replica++){//For loop for replicas
        int samples_replica=0;
        while(samples_replica<sample_inter_swap){//Loop to create samples for each replica until we reach the defined threshold
          
          int temperature_index=index_process(replica);
          total_iterations(i,temperature_index,s)+=1;//increase the number of iterations
          current_temp=temp(temperature_index);// Extract temperature of the replica
          current_log_bound=log_bound_vector(temperature_index);// Extract log-bound of the corresponding temperature
          current_X=X(_,replica);
          bool update_state=true;
 
          // Record the likelihood
          double current_loglik=loglik_R(current_X,Q_mat_R);
          
          if(current_loglik>loglikelihood_visited(temperature_index,s)){//If the current loglik is bigger
            //Store time
            std::clock_t time_find_mode = std::clock();
            double secs_find_mode = static_cast<double>(time_find_mode - start) / CLOCKS_PER_SEC;
            time_visited(temperature_index,s)=secs_find_mode;
            //Store iteration
            iter_to_visit(temperature_index,s)=i;
            //Store value of loglik
            loglikelihood_visited(temperature_index,s)=current_loglik;
            //Store state
            arma::vec current_X_arma = Rcpp::as<arma::vec>(current_X);
            states_visited.slice(s).col(temperature_index) = current_X_arma;
            if(current_loglik>log(11500)){
              Rcpp::Rcout <<"Found higher lik: "<<exp(current_loglik)<<" in iteration: "<<i<<" in temp: "<<temperature_index<<" in time: "<<secs_find_mode<< std::endl; 
            }
          }
          
          
          if(temperature_index<temps_rf){//For the colder temperatures we use Rejection-Free
            //// Visit neighbors in parallel
            NumericVector output_current_X_bounded(p);
            IIT_visit_bounded visit_current_X_bounded(current_X,
                                                      Q_mat_R,
                                                      current_temp,
                                                      theta,
                                                      output_current_X_bounded,
                                                      dim_size,
                                                      dim_size,
                                                      current_log_bound,
                                                      bal_func);
            
            parallelFor(0,dim_size,visit_current_X_bounded);//Apply ParallelFor
            //// Add the h(piy/pix) to compute Z factor
            SumExp get_sum(output_current_X_bounded);
            parallelReduce(0,dim_size,get_sum);
            Z=get_sum.Z/p;//Divide over number of neihgbors
            //// Compute weight
            new_samples=1+R::rgeom(Z);
            if(new_samples<1){
              Rcpp::Rcout <<"Error: geometric in "<< "simulation: " << s+startsim << " Swap: " << i <<" temperature:"<<current_temp<< std::endl;
              Rcpp::Rcout <<"new_samples= "<<new_samples<< ", Z=" << Z << " log-bound= " << current_log_bound << std::endl;
              new_samples=sample_inter_swap;
            }
            if((samples_replica+new_samples)>sample_inter_swap){//If we're going to surpass the required number of samples
              new_samples = sample_inter_swap-samples_replica;//We force to stop at sample_inter_swap
              update_state=false;// We stay in the current state
            }
            
            
            if(update_state){
              ////Sample Proportionally
              //Get random uniforms
              arma::Col<double> u_random(p,fill::randu);
              NumericVector u=Rcpp::wrap(u_random);
              //Compute the needed values
              NumericVector choose_min=log(-log(u)) - (output_current_X_bounded);//We can sample from
              //Find the index of the minimum entry
              GetMin min_coord(choose_min);
              parallelReduce(0,dim_size,min_coord);
              //Swap that coordinate
              X(min_coord.min_index,replica)=1-X(min_coord.min_index,replica);
            }//End If for updating state
            
          }else{//For the hotter temperatures we use step by step update
            //// Process to perform step by step instead of rejection free steps     
            // single_step_update(NumericVector currentX, NumericMatrix Q,int p, int bal_func, double current_temp, double theta, double current_log_bound)
            single_step_output=single_step_update(current_X,Q_mat_R,p,bal_func,current_temp,theta,current_log_bound);
            new_samples=1;
            bool accept_jump=single_step_output(0);
            int chosen_coord_single = single_step_output(1);
            if(accept_jump){
              X(chosen_coord_single,replica)=1-X(chosen_coord_single,replica);   
            }
          }
          
          ///// Updating before the next iteration of the loop
          samples_replica+=new_samples; // Update number of samples obtained from the replica 
        }//End loop to update a single replica
      }//End loop to update all replicas
      
      //// Start replica swap process
      swap_count+=1;//Increase the count of swaps
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
        //Define replicas to swap
        auto find_t = std::find(index_process.begin(), index_process.end(), t);
        auto find_t_1 = std::find(index_process.begin(), index_process.end(), t+1);
        size_t index_t=std::distance(index_process.begin(), find_t);
        size_t index_t_1=std::distance(index_process.begin(), find_t_1);
        Xtemp_from=X(_,index_t);
        Xtemp_to=X(_,index_t_1);
        
        //// Computing swap probability
        swap_prob=(temp(t)-temp(t+1))*(loglik_R(Xtemp_to,Q_mat_R) - loglik_R(Xtemp_from,Q_mat_R));
        swap_prob=exp(swap_prob);
        
        // Rcpp::Rcout <<"Swap prob "<< swap_prob << std::endl;
        ppp=randu();
        if(ppp<swap_prob){//In case the swap is accepted
          swap_success(t)+=1;//Increase the number of successful swaps of temp t
          do_swap.elem(find(index_process==t)).ones();
          do_swap.elem(find(index_process==t+1)).ones();
        }
      }
      //Update index process
      resulting_swap=epsilon_indic % prop_swap % do_swap;
      index_process+=resulting_swap;
      ind_pro_hist.row((s*total_swaps)+swap_count)=index_process.t();
      ////End of replica swap process
      
    }// End loop of iterations
    std::clock_t end_time = std::clock(); // Stop timer
    // Calculate the time taken in seconds
    double duration = static_cast<double>(end_time - start) / CLOCKS_PER_SEC;
    time_taken[s] = duration;
    // Store result of the simulation
    // vec temp_rate=swap_success / swap_total;
    vec temp_rate=swap_success / final_swap;
    Rcpp::Rcout << "Swap count: " << swap_count<< std::endl;
    Rcpp::Rcout << "Final swap: " << final_swap<< std::endl;
    swap_rate.row(s)=temp_rate.t();
  }//End loop simulations
  List ret;
//Print found likelihoods
  Rcpp::Rcout <<"Found liks: \n"<<exp(loglikelihood_visited)<< std::endl; 
  ret["ip"]=ind_pro_hist;
  ret["swap_rate"]=swap_rate;
  ret["loglik_visited"]=loglikelihood_visited;
  ret["iter_visit"]=iter_to_visit;
  ret["time_visit"]=time_visited;
  ret["state_visit"]=states_visited;
  
  ret["total_iter"]=total_iterations;
  ret["time_taken"]=time_taken;

  ret["max_bounds"]=max_log_bound_vector;
  ret["final_bounds"]=log_bound_vector;

  

  return ret;
}
