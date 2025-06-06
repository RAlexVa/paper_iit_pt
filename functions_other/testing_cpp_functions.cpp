//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <thread>
#include <iostream>
#include <random>
#include <execution>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

////////// Other functions //////////

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

