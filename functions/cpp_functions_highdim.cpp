//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <fstream>
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



////////// Creating sparse matrix //////////

// void readSparseMatrix(const std::string& filename)
// [[Rcpp::export]]
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
  
  // Read the subsequent lines for the non-zero entries
  size_t row, col;
  double value;
  for (size_t i = 0; i < n_nonzeros; ++i) {
    file >> row >> col >> value;
    // Adjust for 1-based indexing in the file
    matrix(row - 1, col - 1) = value;
  }
  file.close();
  // Rcpp::Rcout <<"Matrix\n " <<matrix<< std::endl;
  return matrix;
}


// double loglik(vec X, mat M){
//   return(arma::as_scalar(X.t() * M * X));
// }


// [[Rcpp::export]]
void testing_loglik(const std::string& filename, vec X){
  sp_mat M=readSparseMatrix(filename);
  double temp;
  int n=X.n_rows;
  vec tempX(n);
  for(int i=0;i<n;i++){
    tempX=X;
    tempX(i)=1-tempX(i);
    temp=arma::as_scalar(X.t() * M * X);
    Rcpp::Rcout <<"Vector\n " <<tempX.t()<< std::endl;
    Rcpp::Rcout <<"loglik:  " <<temp<< std::endl;
  }
}