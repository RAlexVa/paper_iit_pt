//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <thread>
#include <iostream>
#include <random>
#include <execution>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace RcppParallel;



//////// Other functions //////////

// [[Rcpp::export]]
vec num_to_vec(int n, int d){
  vec X(d);
  X.zeros();
  int temp;
  if(n<0 || n>=std::pow(2,d)){
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

// Return maximum of 3 numbers
// [[Rcpp::export]]
cube sum_cube(){
  cube test_cube(2,3,4,fill::randu);
  // cube test_cube(2,3,4,fill::ones);
  mat slice3=test_cube.slice(2);
  Rcpp::Rcout <<"Sum: "<< sum(slice3.col(1)) << std::endl;
  // Rcpp::Rcout <<"Slice 3: "<< slice3 << std::endl;
  return(test_cube);
}


// [[Rcpp::export]]
void entries_vec(uword& replica, vec& vector){
  uword replica_left;
  uword replica_right;
  if(replica==0){
    replica_left=vector.size()-1;
    replica_right=1;
  }else{
    if(replica==vector.size()-1){
      replica_left=replica-1;
      replica_right=0;
    }else{
      replica_left=replica-1;
      replica_right=replica+1;
    }
  }
  Rcpp::Rcout <<"Replica left: \n"<< vector(replica_left) << std::endl;
  Rcpp::Rcout <<"Replica right: \n"<< vector(replica_right) << std::endl;
}



// [[Rcpp::export]]
void find_min(vec& X, double& new_value){
  uword mm=X.index_min();
  X(mm)=new_value;
  
}


// [[Rcpp::export]]
void print_geom(double& Z){
  Rcpp::Rcout << "Geom value: " << R::rgeom(Z) << std::endl;
}

////////// testing functions //////////

// [[Rcpp::export]]
double test_update(double log_bound,bool update, double prob_to_dec, double temperature, double decreasing_constant){
  if(update){//If it's defined to keep updating the bounding constant
    //Then try Arithmetic reduction of the bounding constant
    if(prob_to_dec>0){//If we consider a probability to decrease the constant
      double test_prob=0;
      if(prob_to_dec<1){
        vec ppp=Rcpp::runif(1);//Get a random number
        test_prob=ppp(0);
      }
      if(test_prob<prob_to_dec){//In case the update is accepted
        double temporal_log_b=log_bound/temperature;
        // Rcpp::Rcout <<"log_b= "<<log_bound<<" temp: "<<temperature<<" temporal_log_b: "<< temporal_log_b<<std::endl;
        double delta_bound = decreasing_constant/exp(temporal_log_b);
        // Rcpp::Rcout <<"delta bound: "<< delta_bound<<std::endl;
        //log(a-b)=log(a)+log(1-b/a) ≈ log(a) - b/a if b/a is very small
        //a=exp(log_bound), b=delta_bound
        if(delta_bound<.06){
          log_bound=temperature*(temporal_log_b - delta_bound);
        }else{
          log_bound = temperature*log(exp(temporal_log_b)-(decreasing_constant)); //Reduce constant
        }
        
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
  
  return log_bound;
}

// [[Rcpp::export]]
void print_log_bound(int iterations, double initial_bound, double prob_to_dec, double temperature, double decreasing_constant){
  Rcpp::Rcout <<"Initial bound: "<< initial_bound<<std::endl;
  double current_bound=initial_bound;
  for(int i=0; i<iterations;i++){
    current_bound=test_update(current_bound,true,prob_to_dec, temperature, decreasing_constant);
    Rcpp::Rcout <<(i+1)<<" New bound: "<< current_bound<<std::endl;
  }
}

// [[Rcpp::export]]
void check_threads(){
  unsigned num_threads = std::thread::hardware_concurrency();
  Rcpp::Rcout <<"Number of threads in local: " <<num_threads<< std::endl;
}

// [[Rcpp::export]]
void check_rng(){
  std::mt19937 gen(1);
  
  // Seed the engine with an unsigned int
  // gen.seed(1);
  std::cout << "1. after seed by 1: " << gen() << '\n';
  std::cout << "2. after seed by 1: " << gen() << '\n';
  std::cout << "3. after seed by 1: " << gen() << '\n';
  std::cout << "4. after seed by 1: " << gen() << '\n';
  
  // Seed the engine with two unsigned ints
  std::seed_seq sseq{1, 2};
  gen.seed(sseq);
  std::cout << "after seed by {1,2}: " << gen() << '\n';
}

// [[Rcpp::export]]
NumericMatrix checking_num(int p){
  //// Define two modes
  NumericMatrix Q_matrix(p,2);
  for(int i=0;i<p;i++){
    if(i%2==0){Q_matrix(i,1)=1;}
    if(i%2==1){Q_matrix(i,0)=1;}
  }
  return Q_matrix;
}

arma::Mat<double> initializeRandom(const int num_rows,const int num_cols, const double prob) {
  
  // Initialize a matrix with random values between 0 and 1
  arma::mat A(num_rows, num_cols,fill::randu);
  
  // Threshold the random values to 0 or 1
  A = arma::conv_to<arma::Mat<double>>::from(A > prob);
  
  return(A);
}


// [[Rcpp::export]]
NumericMatrix check_x(int p, int T){
  
  double ppp=randu();
  arma::Mat<double> X(p,T);
  X=initializeRandom(p,T,ppp);//Randomly initialize the state of each replica.
  return Rcpp::wrap(X);
  
  
  
}



/////////////////////////////////////////////////
//// benchmark performance of updates
/////////////////////////////////////////////////
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



/////////////////////////////////////////////////////////////////////////////////
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
    Rcpp::Rcout <<"M: " <<M<< std::endl;
    return(-10000);
  }
  
  
  
}



///// Here include the for in the parallelization AND make it reproducible /////
//// Definition of worker
struct Parallel_replica_update_rep : public Worker
{
  //// Data
  RMatrix<double> X_n_matrix;//State matrix
  const RMatrix<double> Q_matrix; //Matrix with modes
  const int bal_func; //Specified balancing function
  const RVector<double> temperatures;//Vector of temperatures
  const double theta; //Parameter for
  RMatrix<double> output;
  const std::size_t dim_p; // dimension of the problem (and length of X columns)
  const std::size_t num_temp; // dimension of the problem (and length of X columns)
  const int num_iter;//Number of iterations in the for loop
  std::vector<std::mt19937_64>& rngs;//Vector of RNGs
  //// Functions to use from outside
  arma::Col<double> IIT_update_p(arma::Col<double> X,
                                 const arma::Mat<double>& M,
                                 const int chosen_bf,
                                 const double& temperature,
                                 const double& theta,
                                 std::mt19937_64& rng);
  
  ////  Functions to transform data
  arma::Mat<double> convert_X(){
    RMatrix<double> tmp_matrix = X_n_matrix;
    arma::Mat<double> MAT(tmp_matrix.begin(), dim_p, num_temp, false,true);
    return MAT;
  }
  arma::Mat<double> convert_Q(){
    RMatrix<double> tmp_matrix = Q_matrix;
    arma::Mat<double> MAT(tmp_matrix.begin(), dim_p, 2, false,true);
    return MAT;
  }
  arma::Col<double> convert_temps(){
    RVector<double> tmp_vec = temperatures;
    arma::Col<double> vec(tmp_vec.begin(),num_temp,false,true);
    return vec;
  }
  
  //// Main constructor
  Parallel_replica_update_rep(
    NumericMatrix X_n_matrix_in,
    const NumericMatrix Q_matrix_in,
    const int bal_func_in,
    const NumericVector temperature_in,
    const double theta_in,
    NumericMatrix output_in,
    const std::size_t dim_p_in,
    const std::size_t num_temp_in,
    const int num_iter_in,
    std::vector<std::mt19937_64>& rngs_in)://This binds the class members (defined above) to the constructor arguments
    X_n_matrix(X_n_matrix_in),
    Q_matrix(Q_matrix_in),
    bal_func(bal_func_in),
    temperatures(temperature_in),
    theta(theta_in),
    output(output_in),
    dim_p(dim_p_in),
    num_temp(num_temp_in),
    num_iter(num_iter_in),
    rngs(rngs_in)
  {}
  
  //// Operator
  void operator()(std::size_t begin, std::size_t end){
    //First transform all the inputs
    arma::Mat<double> X = convert_X();
    arma::Mat<double> Q_matrix = convert_Q();
    arma::Col<double> temperatures = convert_temps();
    int dim_p_int = static_cast<int>(dim_p);
    
    for(std::size_t r=begin; r<end;r++){ // For for each replica
      double current_temp = temperatures(r);
      arma::Col<double> selected_X=X.col(r); // Get the status of the replica
      arma::Col<double> temporal_X(dim_p);
      std::mt19937_64 selected_rng=rngs[r];
      for(int k=0;k<num_iter;k++){//Repeat for as many iterations
        temporal_X=IIT_update_p(selected_X,
                                Q_matrix,
                                bal_func,
                                current_temp,
                                theta,
                                selected_rng);
        selected_X=temporal_X;//Assign new state for the loop
      }
      // After num_iter updates, export the resulting file
      for(int i=0;i<dim_p_int;i++){
        output(i,r)=temporal_X[i];
      }//End loop for the export
    }// End loop for each replica
    
  }//End operator
};// End of worker

// Full definition of update function used inside the worker
arma::Col<double> Parallel_replica_update_rep::IIT_update_p(arma::Col<double> X,
                                                            const arma::Mat<double>& M,
                                                            const int chosen_bf,
                                                            const double& temperature,
                                                            const double& theta,
                                                            std::mt19937_64& rng){
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
  // vec u(total_neighbors,fill::randu);
  vec u(total_neighbors);
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  for(int j=0; j<total_neighbors; j++) {
    u(j) = dist(rng);
  }
  vec probs_choose = log(-log(u)) - probs;
  int neigh_pos = (std::min_element(probs_choose.begin(), probs_choose.end()))-probs_choose.begin();
  X.row(neigh_pos) = 1-X.row(neigh_pos); //modify the coordinate of the chosen neighbor
  return X;
}




// [[Rcpp::export]]
void check_min_find(vec X, int value){
  auto find_value = std::find(X.begin(), X.end(), value);
  Rcpp::Rcout << "Find value " <<find_value<< std::endl;
  size_t index_value=std::distance(X.begin(), find_value);
  Rcpp::Rcout << "Index: " <<index_value<< std::endl;
  Rcpp::Rcout << "Value: " <<X[index_value]<< std::endl;
}

// [[Rcpp::export]]
void check_bo0l_vec(int d){
  uvec vec_bool(d);
  vec_bool.fill(0);
  for(int i=0;i<100;i++){
    vec_bool(i)=1;
    
    if(all(vec_bool)){
      Rcpp::Rcout << "Finish loop " << std::endl;
      break;
    }else
    {
      Rcpp::Rcout << "Not finish" << std::endl;
    }
    
    
  }
}








