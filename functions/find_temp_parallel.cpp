//#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;
using namespace arma;

// [[Rcpp::plugins("cpp17")]] 


///// List of balancing functions that apply to log-likelihoods
double ret_min(double a,double b,double c){
  arma::Col<double> temp={a,b,c};
  return(min(temp));
}
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
double bound_sq(double x, double log_gamma){
  double temp1=x/2-log_gamma;
  return ret_min(temp1,x,0);
}
//This function we don't export to R because it generates errors
double invoke(double x, double (*func)(double)) {
  return func(x);
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


int L1_distance(const arma::Col<double>& v1, const arma::Col<double>& v2) {
  // Check if the vectors have the same size
  if (v1.n_elem != v2.n_elem) {
    Rcpp::Rcout <<"Vectors are not the same size" << std::endl;
    return -1;
  }else{
    int dist=sum(abs(v1-v2));
    return dist;
  }
}


// [[Rcpp::export]]
double loglik(const arma::Col<double>& X,const arma::Mat<double>& M, const double& theta){
  double loglik_computed;
  if(M.n_cols==2){
    //Define vectors to represent each mode
    arma::Col<double> mod1=M.col(0);
    arma::Col<double> mod2=M.col(1);
    // Compute distances to each mode  
    double dif1=L1_distance(X,mod1);
    double dif2=L1_distance(X,mod2);
    
    if(dif2<=dif1){
      loglik_computed = -dif2*theta + log1p(exp((dif2-dif1)*theta));
    }
    if(dif2>dif1){
      loglik_computed = -dif1*theta + log1p(exp((dif1-dif2)*theta));
    }
    return(loglik_computed);
    
  }else{
    Rcpp::Rcout << "Error matrix number of columns isn't 2" << std::endl;
    return(-10000);
  }

  
  
}

arma::Mat<double> initializeRandom(const int num_rows,const int num_cols, const double prob) {

  // Initialize a matrix with random values between 0 and 1
  arma::mat A(num_rows, num_cols,fill::randu);

  // Threshold the random values to 0 or 1
  A = arma::conv_to<arma::Mat<double>>::from(A > prob);

  return(A);
}


//// Definition of worker
struct Parallel_replica_update : public Worker
{
  //// Data
  // RMatrix<double> X_n_matrix;//State matrix
  arma::Mat<double> X_n_matrix;//State matrix
  // const RMatrix<double> Q_matrix; //Matrix with modes
  const arma::Mat<double> Q_matrix; //Matrix with modes
  const int bal_func; //Specified balancing function
  // const RVector<double> temperatures;//Vector of temperatures
  const arma::Col<double> temperatures;//Vector of temperatures
  const double theta; //Parameter for 
  // RMatrix<double> output;
  arma::Mat<double> output;
  const std::size_t dim_p; // dimension of the problem (and length of X columns)
  const std::size_t num_temp; // dimension of the problem (and length of X columns)
  
  
  //// Functions to use from outside
  arma::Col<double> IIT_update_p(arma::Col<double> X,
                    const arma::Mat<double>& M,
                    const int chosen_bf,
                    const double& temperature,
                    const double& theta);
  
  ////  Functions to transform data
  // arma::Mat<double> convert_X(){
  //   RMatrix<double> tmp_matrix = X_n_matrix;
  //   arma::Mat<double> MAT(tmp_matrix.begin(), dim_p, num_temp, false,true);
  //   return MAT;
  // }
  // 
  // arma::Mat<double> convert_Q(){
  //   RMatrix<double> tmp_matrix = Q_matrix;
  //   arma::Mat<double> MAT(tmp_matrix.begin(), dim_p, num_temp, false,true);
  //   return MAT;
  // }
  // 
  // arma::Col<double> convert_temps(){
  //   RVector<double> tmp_vec = temperatures;
  //   arma::Col<double> vec(tmp_vec.begin(),num_temp,false,true);
  //   return vec;
  // }
  
  
  
  //// Main constructor
  Parallel_replica_update(
    // NumericMatrix X_n_matrix_in,
    arma::Mat<double> X_n_matrix_in,
    // const NumericMatrix Q_matrix_in,
    const  arma::Mat<double> Q_matrix_in,
    const int bal_func_in,
    // const NumericVector temperature_in,
    const arma::Col<double> temperature_in,
    const double theta_in,
    // NumericMatrix output_in,
    arma::Mat<double> output_in)://This binds the class members (defined above) to the constructor arguments
    X_n_matrix(X_n_matrix_in),
    Q_matrix(Q_matrix_in),
    bal_func(bal_func_in),
    temperatures(temperature_in),
    theta(theta_in),
    output(output_in),
    dim_p(X_n_matrix.n_rows),
    num_temp(X_n_matrix.n_cols)
    // dim_p(X_n_matrix.n_rows), //I'm not sure how to get the dim
    // num_tmep(x_n_matrix.n_cols) //I'm not sure how to get the dim
  {}
  
  
  
  //// Operator
  void operator()(std::size_t begin, std::size_t end){
    //First transform all the inputs
    // arma::Mat<double> X = convert_X();
    // arma::Mat<double> Q_matrix = convert_Q();
    // arma::Col<double> temperatures = convert_temps();
    
    
    for(std::size_t r=begin; r<end;r++){ // For for each replica
      double current_temp = temperatures(r);
      arma::Col temporal_X=IIT_update_p(X_n_matrix,
                                        Q_matrix,
                                        bal_func,
                                        current_temp,
                                        theta);
        output.col(r)=temporal_X;
    }
    
  }//End operator
};// End of worker

// Full definition of update function used inside the worker
arma::Col<double> Parallel_replica_update::IIT_update_p(arma::Col<double> X,
                                                        const arma::Mat<double>& M,
                                                        const int chosen_bf,
                                                        const double& temperature,
                                                        const double& theta){
  int total_neighbors = X.n_rows; // total number of neighbors is p spacial
  vec probs(total_neighbors, fill::zeros); //vector of probabilities
  double logpi_current = loglik(X,M,theta);  // Compute likelihood of current state
  
  ////// Compute weight for all neighbors
  double temporal=0;
  vec newX;
  for(int j=0; j<total_neighbors;j++){
    newX = X;
    newX.row(j) = 1-X.row(j);
    temporal=loglik(newX,M,theta)-logpi_current;
    probs(j)=apply_bal_func(temporal*temperature, chosen_bf);//Apply balancing function to log probability times temperature
  }
  //////Choose the next neighbor
  vec u(total_neighbors,fill::randu);
  vec probs_choose = log(-log(u)) - probs;
  int neigh_pos = (std::min_element(probs_choose.begin(), probs_choose.end()))-probs_choose.begin();
  X.row(neigh_pos) = 1-X.row(neigh_pos); //modify the coordinate of the chosen neighbor
  return X;
}




//Function to run simulations
// [[Rcpp::export]]
arma::Mat<double> PT_IIT_parallel_sim(int p, int num_iter, arma::Col<double> temperatures, const std::string bal_func, double theta){
  int T=temperatures.n_rows;
  arma::Mat<double> output(p,T);
  // Change String of balancing function into int
  int chosen_bal_func;
  if(bal_func=="sq"){
    chosen_bal_func=2;
  }else if(bal_func=="min"){
    chosen_bal_func=1;
  }else{
    Rcpp::Rcout << "Incorrect definition of Balancing function.\n Defaulting to sq " << std::endl;
    chosen_bal_func=2;
  }
  
  //// Define two modes
  arma::Mat<double> Q_matrix(p,2);
  for(int i=0;i<p;i++){
    if(i%2==0){Q_matrix(i,1)=1;}
    if(i%2==1){Q_matrix(i,0)=1;}
  }
  
  double ppp=randu();
  arma::Mat<double> X(p,T);
  X=initializeRandom(p,T,ppp);//Randomly initialize the state of each replica.
  
  
  
  for(int i=0;i<num_iter;i++){
    arma::Mat<double> output(p,T);
    
    Parallel_replica_update par_rep_update(
        X,Q_matrix,
        chosen_bal_func,
        temperatures,
        theta,
        output);// Initialize the worker
    
    const std::size_t begin = static_cast <size_t> (0);
    const std::size_t end = static_cast <size_t> (T); 
    
    parallelFor(begin,end,par_rep_update);//Perform the for
      X=output;  
  }
  
  return X;
}



