//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

////////// Other functions //////////


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



////////// testing functions //////////
// [[Rcpp::export]]
void processStrings(const std::vector<std::string>& inputVector) {
  std::cout << bal_func(144,inputVector[0]) << std::endl;
  std::cout << inputVector[0] << std::endl;
  for (const auto& str : inputVector) {
    std::cout << bal_func(54,str) << std::endl;
  }
}

// [[Rcpp::export]]
mat testassignment(int p, vec temp,const std::vector<std::string>& bal_function){
  int T=temp.n_rows; // Count number of temperatures
  mat X(p,T);
  X.zeros();
  double current_temp;
  vec index_process(T); 
  List output;
  for(int i=0;i<T;i++){ // Reset index process vector at the start of each simulation
    // Rcpp::Rcout <<"Fills index process "<< i+1 << std::endl;
    index_process.row(i)=i;
  }
  
  for(int replica=0;replica<T;replica++){//For loop for replicas
    current_temp=temp(index_process(replica));// Extract temperature of the replica
    // Rcpp::Rcout <<"Inside replica loop "<< replica << std::endl;
    //Depending on the chosen method
    //// Update each replica independently
    output=IIT_update_w(X.col(replica),bal_function[replica],current_temp);
    X.col(replica)=vec(output(0));
  }
  return X;
}