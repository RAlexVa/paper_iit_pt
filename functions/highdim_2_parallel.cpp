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
// [[Rcpp::export]]
arma::Mat<double> initializeRandom(const int num_rows,const int num_cols, const double prob) {
  
  // Initialize a matrix with random values between 0 and 1
  arma::mat A(num_rows, num_cols,fill::randu);
  
  // Threshold the random values to 0 or 1
  A = arma::conv_to<arma::Mat<double>>::from(A > prob);
  
  return(A);
}
// [[Rcpp::export]]
arma::Mat<double> initializeRandom_w_modes(const int num_rows,const int num_cols, const mat mode_matrix) {
  arma::mat A(num_rows, num_cols,fill::zeros);
  int check_nrows=mode_matrix.n_rows;
  if((num_rows!=check_nrows)){
    Rcpp::Rcout << "Error: Number of rows dont match" << std::endl;
    return(A);}
  
  int number_modes=mode_matrix.n_cols;
  
  for(int r=0;r<num_cols;r++){
    // double rand_mode=number_modes*randu();
    // int choose_mode=rand_mode*1;
    int choose_mode = arma::randi<int>(arma::distr_param(0, number_modes-1));
    
    A.col(r)=mode_matrix.col(choose_mode);//State close to that mode
    int N;
    if(num_rows<1000){
      N= arma::randi<int>(arma::distr_param(1,num_rows));
    }else{
      N = arma::randi<int>(arma::distr_param(100,300));//Number of coords to flip
      
    }
    // Rcpp::Rcout << "N:" <<N<< std::endl;
    Rcpp::Rcout << "Chosen mode:" <<choose_mode<<", Change coords: "<<N<< std::endl;
    uvec indices = arma::randperm(num_rows, N);//Choose coords to flip
    
    
    for (arma::uword i = 0; i < indices.n_elem; ++i) {
      A(indices(i),r) = 1 - A(indices(i),r);
    }
  }
  return(A);
}


// [[Rcpp::export]]
arma::Mat<double> create_mode_matrix(const int num_rows,const int num_cols) {
  arma::mat Q(num_rows, num_cols,fill::zeros);
  if((num_cols>=1) && (num_cols<=7)){
for(int c=0;c<num_rows;c++){
      if(c%2==0){Q(c,0)=1;}//First mode index 0
      if(c%2==1){Q(c,1)=1;}//second mode index 1
      
      if(num_cols>2){
        if(c<(num_rows/2)){//First half is 1
          Q(c,2)=1;//Third mode index 2
        }
        if(num_cols>3){
          if(c>=(num_rows/2)){//Second half is 1
            Q(c,3)=1;//Fourth mode index 3
          }
          if(num_cols>4){
            if((c>=(num_rows/4)) && (c<(num_rows*3/4))){//First and last quarter are 0
              Q(c,4)=1;//Fifth mode index 4
            }
            if(num_cols>5){
              if((c<(num_rows/4)) || (c>=(num_rows*3/4))){//First and last quarter are 1
                Q(c,5)=1;//Sixth mode index 5
              }
              if(num_cols>6){
                Q(c,6)=1;  //All 1s
                //Seventh mode index 6
              }
            }
          }
        }
      }
      
    }//End for loop for columns
  }else{//Check number of temperatures
    Rcpp::Rcout << "Error: number of modes is not between 1 and 7" << std::endl;
    Q.fill(-1);
  }
  return(Q);
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
double loglik_exptail(const arma::vec& X,const arma::mat& M,const double& theta){
  // double theta=3;
  
  if(M.n_cols!=2){
    double lik_computed=0;
    //Each column in M is a mode
    for(uword c=0;c<M.n_cols;c++){//For loop for modes
      double dist_mode=sum(abs(X-M.col(c)))*theta;
      if(dist_mode<708){//If to avoid underflowing
        lik_computed+=exp(-dist_mode);//Add exp
      }
    }
    
    return(log1p(lik_computed));//Return 1+log()
  }else{
    double loglik_computed;
    vec mod1=M.col(0);
    vec mod2=M.col(1);
    double dif1=sum(abs(X-mod1));
    double dif2=sum(abs(X-mod2));

    if(dif2<=dif1){
      loglik_computed = -dif2*theta + log1p(exp((dif2-dif1)*theta));
    }
    if(dif2>dif1){
      loglik_computed = -dif1*theta + log1p(exp((dif1-dif2)*theta));
    }
   return(loglik_computed);
  }
  
}
double loglik_R_exptail(NumericVector& X,const NumericMatrix& M, const double& theta){
  
  if(M.ncol()==2){
    double loglik_computed;
    //Define vectors to represent each mode
    NumericVector mod1=M(_,0);
    NumericVector mod2=M(_,1);
    // Compute distances to each mode  
    double dif1=sum(abs(X-mod1));
    double dif2=sum(abs(X-mod2));
    
    if(dif2<=dif1){
      loglik_computed = -dif2*theta + log1p(exp((dif2-dif1)*theta));
    }
    if(dif2>dif1){
      loglik_computed = -dif1*theta + log1p(exp((dif1-dif2)*theta));
    }
    return(loglik_computed);
    //End part for 2 modes
  }else{
    //Each column in M is a mode
    double lik_computed;
    for(int c=0;c<M.ncol();c++){//For loop for modes
      double dist_mode=sum(abs(X-M(_,c)))*theta;
      if(dist_mode<708){//If to avoid underflowing
        lik_computed+=exp(-dist_mode);//Add exp
      }
    }
    
    return(log1p(lik_computed));

  }//End part for more than 2 modes
  
  
  
}
// [[Rcpp::export]]
double loglik(const arma::vec& X,const arma::mat& M,const double& theta){
  // double theta=3;
  double max_dist=400.0;//How long do we allow the tail to go
  bool close_enough=false;  
    double loglik_computed=-0.0;
    double loglik_second_part=-0.0;
    // Rcpp::Rcout <<"Created Lik: "<<lik_computed<< std::endl;
    //Each column in M is a mode
    for(uword c=0;c<M.n_cols;c++){//For loop for modes
      double dist_mode=arma::accu(abs(X-M.col(c)));
      // Rcpp::Rcout <<"mode: "<<c<<" dist: "<<dist_mode<< std::endl;
      //&& (dist_mode*theta)<706
      if(dist_mode<max_dist ){
        // Check the distance, so there's space between modes
        //Check that computing exp of the log-lik doesnt underflow
        loglik_computed-=(dist_mode*theta);
        close_enough=true;
//IMPORTANT: with max_dist<p/2
//at most 1 mode will contribute to the log-likelihood
//That's why we can do this.
      }else{
        loglik_second_part-=(dist_mode*theta);
      }
  
    }//End for loop for modes
    // Rcpp::Rcout <<"New LogLik: "<<loglik_computed<<"New LogLikx100: "<<loglik_computed*100<< std::endl;
    // Rcpp::Rcout <<"Epxm1: "<<expm1(loglik_computed)<<"Log1p "<<log1p(expm1(loglik_computed)+1)<< std::endl;
    // if(close_enough){
    //   return(log1p(expm1(loglik_computed)+1)); 
    // }else{
    //   return(0);
    // }
    if(close_enough){
      return(loglik_computed + log1p(exp(loglik_second_part)));
    }else{
      return(loglik_second_part);
    }

}
double loglik_R(NumericVector& X,const NumericMatrix& M, const double& theta){
  int dim_x = X.length();
  int row_M = M.nrow();
  int col_M = M.ncol();
  
  
  arma::Col<double> VEC(dim_x);
  arma::Mat<double> MAT(row_M,col_M);
  for(int r=0;r<dim_x;r++){//For loop for rows
    VEC(r)=X(r);
    for(int c=0;c<col_M;c++){
      MAT(r,c)=M(r,c);
    }
  }
  
 return(loglik(VEC,MAT,theta));
 
//   double max_dist=235.0;//How long do we allow the tail to go
//   bool close_enough=false;  
//     //Each column in M is a mode
//     double loglik_computed=0;
//     for(int c=0;c<M.ncol();c++){//For loop for modes
//       double dist_mode=sum(abs(X-M(_,c)));
//       if(dist_mode<max_dist && (dist_mode*theta)<706){//If to avoid underflowing
//         loglik_computed-=(dist_mode*theta);
//         close_enough=true;
// //IMPORTANT: with max_dist<p/2
// //at most 1 mode will contribute to the log-likelihood
// //That's why we can do this.
//       }
//     }
//     if(close_enough){
//       return(log1p(expm1(loglik_computed)+1)); 
//     }else{
//       return(0);
//     }
}



///// Full definition of internal functions of workers 
double IIT_visit_neighbors::apply_bal_func_internal(double x,const int chosen){
  return apply_bal_func(x,chosen);
}

double IIT_visit_neighbors::loglik_internal(const arma::Col<double>& X,const arma::Mat<double>& M, const double& theta){
  return loglik(X,M,theta);
}

double IIT_visit_bounded::apply_bal_func_bounded_internal(double x,double log_bound){
  return bound_sq(x,log_bound);
}

double IIT_visit_bounded::loglik_internal(const arma::Col<double>& X,const arma::Mat<double>& M, const double& theta){
  return loglik(X,M,theta);
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
  //Number of states to keep track
  // mat iter_to_visit(num_states_visited,total_sim);
  // mat loglikelihood_visited(num_states_visited,total_sim);
  // loglikelihood_visited.fill(-10000);//Initialize a very negative loglikelihood
  // cube states_visited(p,num_states_visited,total_sim,fill::zeros);
  Rcpp::Rcout <<"p= "<<p<<" num_states_visited= "<<num_states_visited<<" total_sim= "<<total_sim<< std::endl;
  // mat modes_visited(numiter * total_sim,T);//Matrix to store the modes visited and temperature
  const std::size_t dim_size = static_cast <size_t> (p); 
  const std::size_t number_modes = static_cast <size_t> (num_modes); 
  //// Define matrix with modes
  mat Q_matrix(p,num_modes);
  Q_matrix=create_mode_matrix(p,num_modes);
  // for(int i=0;i<p;i++){
  //   if(i%2==0){Q_matrix(i,1)=1;}
  //   if(i%2==1){Q_matrix(i,0)=1;}
  // }
  NumericMatrix Q_mat_R=Rcpp::wrap(Q_matrix);//Numeric Matrix version of the Q_matrix
  // NumericVector mode1=Q_mat_R(_,0);
  // NumericVector mode2=Q_mat_R(_,1);
  NumericMatrix distance_modes(num_modes,T);
  // NumericMatrix distance_mode1(total_sim,T);
  // NumericMatrix distance_mode2(total_sim,T);
  // distance_mode1.fill(100000);
  // distance_mode2.fill(100000);
  distance_modes.fill(100000); //To store distance to modes for each simulation
  
  cube distance_modes_full(num_modes,T,total_sim); //Final report of distance to modes
  // NumericMatrix full_distance_mode1(total_sim,T,NA);
  // NumericMatrix full_distance_mode2(total_sim,T,NA);
  // bool first_find_m1=false;
  // bool first_find_m2=false;
  // NumericMatrix time_find_m1(total_sim,T);
  // NumericMatrix time_find_m2(total_sim,T);
  NumericMatrix time_find_modes(num_modes,T);
  // time_find_m1.fill(-1);
  // time_find_m2.fill(-1);
  time_find_modes.fill(-1);
  cube time_find_modes_full(num_modes,T,total_sim);
  
  uvec check_mode_visit(num_modes);//Vector to break the loop when temp1 visits the mode
  
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
    inter_mat=initializeRandom_w_modes(p,T,Q_matrix);
    X=Rcpp::wrap(inter_mat);
    initialX=clone(X);
    swap_total.zeros();
    swap_success.zeros();
    check_mode_visit.fill(0);
    //// Start loop for burn_in period
    for(int i=0;i<burn_in;i++){
      if (i % 1000 == 1) {Rcpp::Rcout << "PT-IIT Simulation: " << s+startsim << " Burn_in period, iteration: " << i << std::endl;}
      for(int replica=0;replica<T;replica++){//For loop for replica update
        current_temp=temp(index_process(replica));
        
        NumericVector output_current_X(p);
        // NumericVector current_X=Rcpp::wrap(X.col(replica));
        NumericVector current_X=X(_,replica);
        
        IIT_visit_neighbors visit_current_X(current_X,
                                            Q_mat_R,
                                            bal_func,
                                            current_temp,
                                            theta,
                                            output_current_X,
                                            dim_size,
                                            number_modes);
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
          // epsilon_indic.elem(find(index_process==t)).ones();
          // prop_swap.elem(find(index_process==t)).ones(); //we swap temperature t
          // prop_swap.elem(find(index_process==t+1)).ones(); //With t+1
          ////Compute swap probability
          // Xtemp_from=Rcpp::wrap(X.col(t));
          // Xtemp_to=Rcpp::wrap(X.col(t+1));
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
                                      number_modes);
            IIT_visit_neighbors Z_v22(Xtemp_to,
                                      Q_mat_R,
                                      bal_func,
                                      temp(t+1),
                                      theta,
                                      output_Z_v22,
                                      dim_size,
                                      number_modes);
            IIT_visit_neighbors Z_v12(Xtemp_from,
                                      Q_mat_R,
                                      bal_func,
                                      temp(t+1),
                                      theta,
                                      output_Z_v12,
                                      dim_size,
                                      number_modes);
            IIT_visit_neighbors Z_v21(Xtemp_to,
                                      Q_mat_R,
                                      bal_func,
                                      temp(t),
                                      theta,
                                      output_Z_v21,
                                      dim_size,
                                      number_modes);
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
          swap_prob=(temp(t)-temp(t+1))*(loglik_R(Xtemp_to,Q_mat_R,theta) - loglik_R(Xtemp_from,Q_mat_R,theta)); 
          swap_prob=Z_fact_correc*exp(swap_prob);
          // Rcpp::Rcout <<"Swap prob "<< swap_prob << std::endl;
          ppp=Rcpp::runif(1);
          if(ppp(0)<swap_prob){//In case the swap is accepted
            //Swap vectors in the matrix X
            for(int coord=0;coord<p;coord++){
              X(coord,t+1)=Xtemp_from[coord];
              X(coord,t)=Xtemp_to[coord];
            }
            // X.col(t+1)=Xtemp_from;
            // X.col(t)=Xtemp_to;
          }//In case a swap is accepted
        }//Finish loop to swap many replicas
      }//Finish IF for replica swap
    }// Finish burn-in period

    /////////////////////// Finish burn-in period
    swap_count=0; //Reset swap count
    std::clock_t start = std::clock(); // Start timer for simulation s
    // Sleep(5000);
    //// Start the loop for all iterations in simulation s
    for(int i=0;i<numiter;i++){
      if (i % 1000 == 1) {Rcpp::Rcout << "PT-IIT Simulation: " << s+startsim << " Iteration: " << i << std::endl;}
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
                                            number_modes);
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
        
        // Rcpp::Rcout <<"Temp: "<<replica<< " Min Coord: " << min_coord.min_index  << std::endl;
        // // Rcpp::Rcout <<"Temp: "<<replica<< " Current X:\n " << output_current_X << std::endl;
        // Rcpp::Rcout <<"idx 72: "<<output_current_X[71]<< " idx 476: " <<output_current_X[475]<<" idx 1: " <<output_current_X[0]<< std::endl;
        // Rcpp::Rcout <<"Minimum : "<<output_current_X[min_coord.min_index]<< std::endl;
        // Rcpp::Rcout <<"Temp: "<<current_temp<< " Now with the rand: " << std::endl;
        // Rcpp::Rcout <<"idx 72: "<<choose_min[71]<< " idx 476: " <<choose_min[475]<<" idx 1: " <<choose_min[0]<< std::endl;
        // Rcpp::Rcout <<"Minimum : "<<choose_min[min_coord.min_index]<< std::endl;
        
//// Check distance to modes
        for(int mode_counter=0;mode_counter<num_modes;mode_counter++){
          double dist_mode=sum(abs(current_X-Q_mat_R(_,mode_counter)));
          if(distance_modes(mode_counter,temperature_index)>dist_mode){//In case we find a smaller distance to the mode
            // Rcpp::Rcout <<"Mode : "<<mode_counter<<" dist : "<<dist_mode<< std::endl;
            distance_modes(mode_counter,temperature_index)=dist_mode;
            std::clock_t time_find_mode = std::clock();
            double secs_find_mode = static_cast<double>(time_find_mode - start) / CLOCKS_PER_SEC;
            time_find_modes(mode_counter,temperature_index)=secs_find_mode;
            if(dist_mode==0){
              check_mode_visit(mode_counter)=1;//Turn to 1 when visit the mode
            }
          }
        }
 
        //Swap that coordinate
        X(min_coord.min_index,replica)=1-X(min_coord.min_index,replica);
        // X.col(replica)=vec(output(0)); //Update current state of the chain
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
                                      number_modes);
            IIT_visit_neighbors Z_v22(Xtemp_to,
                                      Q_mat_R,
                                      bal_func,
                                      temp(t+1),
                                      theta,
                                      output_Z_v22,
                                      dim_size,
                                      number_modes);
            IIT_visit_neighbors Z_v12(Xtemp_from,
                                      Q_mat_R,
                                      bal_func,
                                      temp(t+1),
                                      theta,
                                      output_Z_v12,
                                      dim_size,
                                      number_modes);
            IIT_visit_neighbors Z_v21(Xtemp_to,
                                      Q_mat_R,
                                      bal_func,
                                      temp(t),
                                      theta,
                                      output_Z_v21,
                                      dim_size,
                                      number_modes);
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
          swap_prob=(temp(t)-temp(t+1))*(loglik_R(Xtemp_to,Q_mat_R,theta) - loglik_R(Xtemp_from,Q_mat_R,theta)); 
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
     if(all(check_mode_visit)){//Condition to break the loop when visit all modes
       Rcpp::Rcout << "Found all modes: " << std::endl;
       Rcpp::Rcout << "PT-IIT Simulation: " << s+startsim << " iteration: " << i << std::endl;
       break;
     } 

    }// End loop of iterations
    std::clock_t end_time = std::clock(); // Stop timer
    // Calculate the time taken in seconds
    double duration = static_cast<double>(end_time - start) / CLOCKS_PER_SEC;
    time_taken[s] = duration;
    // Store result of the simulation
    vec temp_rate=swap_success / swap_total;
    swap_rate.row(s)=temp_rate.t();
    // Rcpp::Rcout << "Index Process:\n " <<index_process << std::endl;
    // Rcpp::Rcout << "State X:\n " <<X << std::endl;
    // Rcpp::Rcout <<"Final state "<< X << std::endl;
    //Store the distance and time of visiting modes in the FULL cube
    distance_modes_full.slice(s)=Rcpp::as<arma::mat>(distance_modes);
    time_find_modes_full.slice(s)=Rcpp::as<arma::mat>(time_find_modes);
  }//End loop simulations
  
  List ret;
  // Rcpp::Rcout <<"Dist1\n "<< distance_mode1 << std::endl;
  // Rcpp::Rcout <<"Dist2\n "<< distance_mode2 << std::endl;
  // Rcpp::Rcout <<"Time1\n "<< time_find_m1 << std::endl;
  // Rcpp::Rcout <<"Time2\n "<< time_find_m2 << std::endl;
  ret["ip"]=ind_pro_hist;
  ret["swap_rate"]=swap_rate;
  // ret["states"]=states_visited;
  // ret["loglik_visited"]=loglikelihood_visited;
  // ret["iter_visit"]=iter_to_visit;
  ret["time_taken"]=time_taken;
  ret["distance_modes"]=distance_modes_full;
  // ret["distance_mode1"]=distance_mode1;
  // ret["distance_mode2"]=distance_mode2;
  // ret["time_mode1"]=time_find_m1;
  // ret["time_mode2"]=time_find_m2;
  ret["time_modes"]=time_find_modes_full;
  ret["initial_X"]=initialX;
  // ret["distance_origin"]=full_distance_origin;
  return ret;
}

// [[Rcpp::export]]
List PT_a_IIT_sim(int p,int startsim,int endsim, int total_swaps,int sample_inter_swap,int burn_in, vec temp, const int bal_func,const std::string& filename,int num_states_visited,const std::vector<int>& starting_coord, double decreasing_constant,std::string reduc_model, double theta, int num_modes){
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
  //Number of states to keep track
  mat iter_to_visit(num_states_visited,total_sim);
  // mat loglikelihood_visited(num_states_visited,total_sim);
  // loglikelihood_visited.fill(-10000);//Initialize a very negative loglikelihood
  // cube states_visited(p,num_states_visited,total_sim,fill::zeros);
  // Rcpp::Rcout <<"p= "<<p<<" num_states_visited= "<<num_states_visited<<" total_sim= "<<total_sim<< std::endl;
  // double temporal_loglik;
  // uword found_min; // to find the minimum
  //// Define modes
  mat Q_matrix(p,num_modes);
  Q_matrix=create_mode_matrix(p,num_modes);
  // for(int i=0;i<p;i++){
  //   if(i%2==0){Q_matrix(i,1)=1;}
  //   if(i%2==1){Q_matrix(i,0)=1;}
  // }
  NumericMatrix Q_mat_R=Rcpp::wrap(Q_matrix);//Numeric Matrix version of the Q_matrix
  // NumericVector mode1=Q_mat_R(_,0);
  // NumericVector mode2=Q_mat_R(_,1);
  NumericMatrix distance_modes(num_modes,T);
  // NumericMatrix distance_mode1(total_sim,T);
  // NumericMatrix distance_mode2(total_sim,T);
  // distance_mode1.fill(100000);
  // distance_mode2.fill(100000);
  distance_modes.fill(100000); //To store distance to modes for each simulation
  cube distance_modes_full(num_modes,T,total_sim); //Final report of distance to modes
  
  // NumericMatrix time_find_m1(total_sim,T);
  // NumericMatrix time_find_m2(total_sim,T);
  // time_find_m1.fill(-1);
  // time_find_m2.fill(-1);
  NumericMatrix time_find_modes(num_modes,T);
  time_find_modes.fill(-1);
  cube time_find_modes_full(num_modes,T,total_sim);
  const std::size_t dim_size = static_cast <size_t> (p);
  const std::size_t number_modes = static_cast <size_t> (num_modes); 


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

  uvec check_mode_visit(num_modes);//Vector to break the loop when temp1 visits the mode
  
  
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
    inter_mat=initializeRandom_w_modes(p,T,Q_matrix);
    X=Rcpp::wrap(inter_mat);

    check_mode_visit.fill(0);
    log_bound_vector.zeros();//Reset log-bounds, all log-bounds start at 0
    swap_success.zeros();
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
                                              number_modes);
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

          NumericVector bounded_vector=output_current_X - current_log_bound;
          // Rcpp::Rcout <<"Declaring getmax"<< std::endl;
          SumExp get_sum(bounded_vector);
          // Rcpp::Rcout <<"Before parallelReduce SumExp"<< std::endl;
          parallelReduce(0,dim_size,get_sum);
          // Rcpp::Rcout <<"After parallelReduce SumExp"<< std::endl;
          Z=get_sum.Z/p;//Divide over the number of neighbors
          new_samples=1+R::rgeom(Z);//Get multiplicity list
          if(new_samples<1){
            Rcpp::Rcout <<"Error: geometric in "<< "simulation: " << s+startsim << " Burn-in period after " << track_burn_in <<"simulations,  temp:"<<current_temp<< std::endl;
            Rcpp::Rcout <<"new_samples= "<<new_samples<< ", Z=" << Z << " log-bound= " << current_log_bound << std::endl;
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

        }//End while loop to update repicas
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
        swap_prob=(temp(t)-temp(t+1))*(loglik_R(Xtemp_to,Q_mat_R,theta) - loglik_R(Xtemp_from,Q_mat_R,theta));
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
// Rcpp::Rcout <<"Dist1\n "<< distance_mode1 << std::endl;
// Rcpp::Rcout <<"Dist2\n "<< distance_mode2 << std::endl;
// Rcpp::Rcout <<"Time1\n "<< time_find_m1 << std::endl;
// Rcpp::Rcout <<"Time2\n "<< time_find_m2 << std::endl;
    swap_count=0; //Reset swap count
    std::clock_t start = std::clock(); // Start timer for simulation s
    
    // Sleep(5000);    
    //// Start the loop for all iterations in simulation s
    for(int i=0;i<total_swaps;i++){
    if (i % 100 == 1) {Rcpp::Rcout << "PT A-IIT Simulation: " << s+startsim << " Swap: " << i<<" Prob_decrease_bound: " << prob_to_dec << std::endl;}
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

          NumericVector output_current_X_bounded(p);
          IIT_visit_bounded visit_current_X_bounded(current_X,
                                            Q_mat_R,
                                            current_temp,
                                            theta,
                                            output_current_X_bounded,
                                            dim_size,
                                            number_modes,
                                            current_log_bound);
          
          parallelFor(0,dim_size,visit_current_X_bounded);//Apply ParallelFor

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
            update_state=false;
          }

          ///// Check distance to modes
          for(int mode_counter=0;mode_counter<num_modes;mode_counter++){
            double dist_mode=sum(abs(current_X-Q_mat_R(_,mode_counter)));
            if(distance_modes(mode_counter,temperature_index)>dist_mode){//In case we find a smaller distance to the mode
              distance_modes(mode_counter,temperature_index)=dist_mode;
              std::clock_t time_find_mode = std::clock();
              double secs_find_mode = static_cast<double>(time_find_mode - start) / CLOCKS_PER_SEC;
              time_find_modes(mode_counter,temperature_index)=secs_find_mode;
              if(dist_mode==0){
                check_mode_visit(mode_counter)=1;//Turn to 1 when visit the mode
              }
            }
          }
          


          ///// Updating before the next iteration of the loop
          samples_replica+=new_samples; // Update number of samples obtained from the replica
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
        swap_prob=(temp(t)-temp(t+1))*(loglik_R(Xtemp_to,Q_mat_R,theta) - loglik_R(Xtemp_from,Q_mat_R,theta));
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
      if(all(check_mode_visit)){//Condition to break the loop when visit all modes
        Rcpp::Rcout << "Found all modes: " << std::endl;
        Rcpp::Rcout << "PT A-IIT Simulation: " << s+startsim << " swap: " << i << std::endl;
        break;
      } 
      
    }// End loop of iterations
    std::clock_t end_time = std::clock(); // Stop timer
    // Calculate the time taken in seconds
    double duration = static_cast<double>(end_time - start) / CLOCKS_PER_SEC;
    time_taken[s] = duration;
    // Store result of the simulation
    vec temp_rate=swap_success / swap_total;
    swap_rate.row(s)=temp_rate.t();
    // Rcpp::Rcout <<"Final state "<< X << std::endl;
    distance_modes_full.slice(s)=Rcpp::as<arma::mat>(distance_modes);
    time_find_modes_full.slice(s)=Rcpp::as<arma::mat>(time_find_modes);
  }//End loop simulations
  List ret;
  // Rcpp::Rcout <<"Dist1\n "<< distance_mode1 << std::endl;
  // Rcpp::Rcout <<"Dist2\n "<< distance_mode2 << std::endl;
  // Rcpp::Rcout <<"Time1\n "<< time_find_m1 << std::endl;
  // Rcpp::Rcout <<"Time2\n "<< time_find_m2 << std::endl;
  ret["ip"]=ind_pro_hist;
  ret["swap_rate"]=swap_rate;
  // ret["states"]=states_visited;
  // ret["loglik_visited"]=loglikelihood_visited;
  // ret["iter_visit"]=iter_to_visit;
  ret["total_iter"]=total_iterations;
  ret["time_taken"]=time_taken;
  // ret["modes_visit"]=full_iter_visit_modes;
  // ret["distance_mode1"]=distance_mode1;
  // ret["distance_mode2"]=distance_mode2;
  ret["distance_modes"]=distance_modes_full;
  // ret["distance_origin"]=full_distance_origin;
  ret["max_bounds"]=max_log_bound_vector;
  ret["final_bounds"]=log_bound_vector;
  // ret["time_mode1"]=time_find_m1;
  // ret["time_mode2"]=time_find_m2;
  ret["time_modes"]=time_find_modes_full;
  return ret;
}
