#ifndef WORKERS_H
#define WORKERS_H

#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;
using namespace arma;


// Inside here there are 2 functions
//loglik_internal
//apply_bal_func_internal
//These need to be defined after
struct IIT_visit_neighbors : public Worker
{
  //// Data
  RVector<double> X_n;//State matrix
  const RMatrix<double> Q_matrix; //Matrix with modes
  const int bal_func; //Specified balancing function
  const double temperature;//Chosen temperature
  const double theta; //Parameter for likelihood
  RVector<double> output;
  const std::size_t dim_p; // dimension of the problem (and length of X columns)
  
  //// Functions to use from outside
  double loglik_internal(const arma::Col<double>& X,const arma::Mat<double>& M, const double& theta);
  
  double apply_bal_func_internal(double x,const int chosen);
  ////  Functions to transform data
  arma::Col<double> convert_X(){
    RVector<double> tmp_X = X_n;
    arma::Col<double> VEC(tmp_X.begin(), dim_p, false,true);
    return VEC;
  }
  arma::Mat<double> convert_Q(){
    RMatrix<double> tmp_matrix = Q_matrix;
    arma::Mat<double> MAT(tmp_matrix.begin(), dim_p, 2, false,true);
    return MAT;
  }
  
  //// Main constructor
  IIT_visit_neighbors(
    NumericVector X_n_in,
    const NumericMatrix Q_matrix_in,
    const int bal_func_in,
    const double temperature_in,
    const double theta_in,
    NumericVector output_in,
    const size_t dim_p_in)://This binds the class members (defined above) to the constructor arguments
    X_n(X_n_in),
    Q_matrix(Q_matrix_in),
    bal_func(bal_func_in),
    temperature(temperature_in),
    theta(theta_in),
    output(output_in),
    dim_p(dim_p_in)
  {}
  
  
  
  //// Operator
  void operator()(std::size_t begin, std::size_t end){
    //First transform all the inputs
    arma::Col<double> X = convert_X();
    arma::Mat<double> M = convert_Q();
    // int dim_p_int = X.n_rows;
    double logpi_current=loglik_internal(X,M,theta);
    
    for(std::size_t n = begin; n < end;n++){ // For for each neighbor
      arma::Col<double> temp_X=X;
      temp_X(n)=1-temp_X(n);//Switch the nth coordinate
      double mid_step=loglik_internal(temp_X,M,theta)-logpi_current;
      output[n]=apply_bal_func_internal(mid_step*temperature,bal_func);
      
    }// End for loop
    
  }//End operator
};// End of worker to visit neighbors

struct SumExp : public Worker
{   
  
  const RVector<double> X;// source vector
  double Z; // accumulated value
  
  // constructors
  SumExp(const NumericVector X) : X(X), Z(0) {}
  SumExp(const SumExp& sum, Split) : X(sum.X), Z(0) {}
  
  // accumulate just the element of the range I've been asked to
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      Z += std::exp(X[i]);
    }
  }// End of operator
  
  // join my value with that of another Sum
  void join(const SumExp& rhs) { 
    Z += rhs.Z; 
  }
};//End of worker to compute the Z factor

struct GetMin : public Worker
{   
  
  const RVector<double> X;// input vector
  double min_value; // MIN value found
  std::size_t min_index; // Index of MIN value
  
  // constructors
  GetMin(const NumericVector X) : X(X), min_value(INFINITY),min_index(0) {}
  GetMin(const GetMin& min_cons, Split) : X(min_cons.X), min_value(INFINITY),min_index(0) {}
  
  // accumulate just the element of the range I've been asked to
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      if(X[i]<min_value){
        min_value=X[i];
        min_index=i;
      }
    }
  }// End of operator
  
  // join my value with that of another Sum
  void join(const GetMin& other) {
    if(other.min_value<min_value){
      min_value = other.min_value;
      min_index = other.min_index;
    }
  }
};//End of worker to compute the Z factor


#endif
