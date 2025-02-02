//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

////////// Other functions //////////

// [[Rcpp::export]]
vec num_to_vec(int n, int d){
  vec X(d);
  X.zeros();
  int temp;
  if(n<0 | n>=std::pow(2,d)){
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

////////// testing functions //////////
