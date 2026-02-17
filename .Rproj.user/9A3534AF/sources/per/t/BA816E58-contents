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
  const std::size_t num_modes; //Number of modes, # of columns of Q_matrix
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
    arma::Mat<double> MAT(tmp_matrix.begin(), dim_p, num_modes, false,true);
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
    const size_t dim_p_in,
    const size_t num_modes_in)://This binds the class members (defined above) to the constructor arguments
    X_n(X_n_in),
    Q_matrix(Q_matrix_in),
    bal_func(bal_func_in),
    temperature(temperature_in),
    theta(theta_in),
    output(output_in),
    dim_p(dim_p_in),
    num_modes(num_modes_in){}
  
  
  
  //// Operator
  void operator()(std::size_t begin, std::size_t end){
    //First transform all the inputs
    arma::Col<double> X = convert_X();
    arma::Mat<double> M = convert_Q();
    // int dim_p_int = X.n_rows;
    double logpi_current=loglik_internal(X,M,theta);
    // Rcpp::Rcout <<"Logpicurrent in parallel: "<<logpi_current << std::endl;
    for(std::size_t n = begin; n < end;n++){ // For for each neighbor
      arma::Col<double> temp_X=X;
      temp_X(n)=1-temp_X(n);//Switch the nth coordinate
      double mid_step=loglik_internal(temp_X,M,theta)-logpi_current;
      output[n]=apply_bal_func_internal(mid_step*temperature,bal_func);
      
    }// End for loop
    
  }//End operator
};// End of worker to visit neighbors

struct IIT_visit_bounded : public Worker
{
  //// Data
  RVector<double> X_n;//State matrix
  const RMatrix<double> Q_matrix; //Matrix with modes
  const double temperature;//Chosen temperature
  const double theta; //Parameter for likelihood
  RVector<double> output;
  const std::size_t dim_p; // dimension of the problem (and length of X columns)
  const std::size_t num_modes;
  double log_bound;
  int bal_func;
  
  //// Functions to use from outside
  double loglik_internal(const arma::Col<double>& X,const arma::Mat<double>& M, const double& theta);
  
  double apply_bal_func_bounded_internal(double x, double log_bound, int bal_func);
  ////  Functions to transform data
  arma::Col<double> convert_X_bounded(){
    RVector<double> tmp_X = X_n;
    arma::Col<double> VEC(tmp_X.begin(), dim_p, false,true);
    return VEC;
  }
  arma::Mat<double> convert_Q_bounded(){
    RMatrix<double> tmp_matrix = Q_matrix;
    arma::Mat<double> MAT(tmp_matrix.begin(), dim_p, num_modes, false,true);
    return MAT;
  }
  
  //// Main constructor
  IIT_visit_bounded(
    NumericVector X_n_in,
    const NumericMatrix Q_matrix_in,
    const double temperature_in,
    const double theta_in,
    NumericVector output_in,
    const size_t dim_p_in,
    const size_t num_modes_in,
    double log_bound_in,
    int bal_func_in)://This binds the class members (defined above) to the constructor arguments
    X_n(X_n_in),
    Q_matrix(Q_matrix_in),
    temperature(temperature_in),
    theta(theta_in),
    output(output_in),
    dim_p(dim_p_in),
    num_modes(num_modes_in),
    log_bound(log_bound_in),
    bal_func(bal_func_in){}
  
  
  
  //// Operator
  void operator()(std::size_t begin, std::size_t end){
    //First transform all the inputs
    arma::Col<double> X = convert_X_bounded();
    arma::Mat<double> M = convert_Q_bounded();
    double logpi_current=loglik_internal(X,M,theta);
    
    for(std::size_t n = begin; n < end;n++){ // For for each neighbor
      arma::Col<double> temp_X=X;
      temp_X(n)=1-temp_X(n);//Switch the nth coordinate
      double mid_step=loglik_internal(temp_X,M,theta)-logpi_current;
      output[n]=apply_bal_func_bounded_internal(mid_step*temperature,log_bound,bal_func);
      
    }// End for loop
  }//End operator
};// End of worker to visit neighbors

struct Ratio_probabilities : public Worker
{
  //// Data
  RVector<double> X_n;//State matrix
  const RMatrix<double> Q_matrix; //Matrix with modes
  const double temperature;//Chosen temperature
  const double theta; //Parameter for likelihood
  RVector<double> output;
  const std::size_t dim_p; // dimension of the problem (and length of X columns)
  const std::size_t num_modes; //Number of modes, # of columns of Q_matrix
  //// Functions to use from outside
  double loglik_internal(const arma::Col<double>& X,const arma::Mat<double>& M, const double& theta);
  
  ////  Functions to transform data
  arma::Col<double> convert_X(){
    RVector<double> tmp_X = X_n;
    arma::Col<double> VEC(tmp_X.begin(), dim_p, false,true);
    return VEC;
  }
  arma::Mat<double> convert_Q(){
    RMatrix<double> tmp_matrix = Q_matrix;
    arma::Mat<double> MAT(tmp_matrix.begin(), dim_p, num_modes, false,true);
    return MAT;
  }
  
  //// Main constructor
  Ratio_probabilities(
    NumericVector X_n_in,
    const NumericMatrix Q_matrix_in,
    const double temperature_in,
    const double theta_in,
    NumericVector output_in,
    const size_t dim_p_in,
    const size_t num_modes_in)://This binds the class members (defined above) to the constructor arguments
    X_n(X_n_in),
    Q_matrix(Q_matrix_in),
    temperature(temperature_in),
    theta(theta_in),
    output(output_in),
    dim_p(dim_p_in),
    num_modes(num_modes_in){}
  
  
  
  //// Operator
  void operator()(std::size_t begin, std::size_t end){
    //First transform all the inputs
    arma::Col<double> X = convert_X();
    arma::Mat<double> M = convert_Q();
    // int dim_p_int = X.n_rows;
    double logpi_current=loglik_internal(X,M,theta);
    // Rcpp::Rcout <<"Logpicurrent in parallel: "<<logpi_current << std::endl;
    for(std::size_t n = begin; n < end;n++){ // For for each neighbor
      arma::Col<double> temp_X=X;
      temp_X(n)=1-temp_X(n);//Switch the nth coordinate
      output[n]=temperature*(loglik_internal(temp_X,M,theta)-logpi_current);
    }// End for loop
  }//End operator
};// End of worker to visit neighbors

struct Apply_bal_func_parallel : public Worker
{
  //// Data
  RVector<double> X_n;//State matrix
  RVector<double> output;
  const std::size_t dim_p; // dimension of the problem (and length of X columns)
  double log_bound;
  int bal_func;
  
  //// Functions to use from outside
  double apply_bal_func_bounded_internal(double x, double log_bound, int bal_func);
  ////  Functions to transform data
  arma::Col<double> convert_X_bounded(){
    RVector<double> tmp_X = X_n;
    arma::Col<double> VEC(tmp_X.begin(), dim_p, false,true);
    return VEC;
  }
  
  //// Main constructor
  Apply_bal_func_parallel(
    NumericVector X_n_in,
    NumericVector output_in,
    const size_t dim_p_in,
    double log_bound_in,
    int bal_func_in)://This binds the class members (defined above) to the constructor arguments
    X_n(X_n_in),
    output(output_in),
    dim_p(dim_p_in),
    log_bound(log_bound_in),
    bal_func(bal_func_in){}
  
  
  //// Operator
  void operator()(std::size_t begin, std::size_t end){
    //First transform all the inputs
    arma::Col<double> X = convert_X_bounded();
    
    for(std::size_t n = begin; n < end;n++){ // For for each neighbor
      arma::Col<double> temp_X=X;//Create a copy of the input vector
      //The output is just applying the balancing function to each entry
      output[n]=apply_bal_func_bounded_internal(temp_X(n),log_bound,bal_func);
      
    }// End for loop
  }//End operator
};// End of worker to visit neighbors

struct IIT_visit_new_bound : public Worker//Not used, dont work with 2 outputs
{
  //// Data
  RVector<double> X_n;//State matrix
  const RMatrix<double> Q_matrix; //Matrix with modes
  const double temperature;//Chosen temperature
  const double theta; //Parameter for likelihood
  RVector<double> output;
  RVector<double> output_bound;
  const std::size_t dim_p; // dimension of the problem (and length of X columns)
  const std::size_t num_modes;
  double log_bound;
  int bal_func;

  //// Functions to use from outside
  double loglik_internal(const arma::Col<double>& X,const arma::Mat<double>& M, const double& theta);

  double apply_bal_func_bounded_internal(double x, double log_bound, int bal_func);
  ////  Functions to transform data
  arma::Col<double> convert_X_bounded(){
    RVector<double> tmp_X = X_n;
    arma::Col<double> VEC(tmp_X.begin(), dim_p, false,true);
    return VEC;
  }
  arma::Mat<double> convert_Q_bounded(){
    RMatrix<double> tmp_matrix = Q_matrix;
    arma::Mat<double> MAT(tmp_matrix.begin(), dim_p, num_modes, false,true);
    return MAT;
  }

  //// Main constructor
  IIT_visit_new_bound(
    NumericVector X_n_in,
    const NumericMatrix Q_matrix_in,
    const double temperature_in,
    const double theta_in,
    NumericVector output_in,
    NumericVector output_bound_in,
    const size_t dim_p_in,
    const size_t num_modes_in,
    double log_bound_in,
    int bal_func_in)://This binds the class members (defined above) to the constructor arguments
    X_n(X_n_in),
    Q_matrix(Q_matrix_in),
    temperature(temperature_in),
    theta(theta_in),
    output(output_in),
    output_bound(output_bound_in),
    dim_p(dim_p_in),
    num_modes(num_modes_in),
    log_bound(log_bound_in),
    bal_func(bal_func_in){}



  //// Operator
  void operator()(std::size_t begin, std::size_t end){
    //First transform all the inputs
    arma::Col<double> X = convert_X_bounded();
    arma::Mat<double> M = convert_Q_bounded();
    double logpi_current=loglik_internal(X,M,theta);

    for(std::size_t n = begin; n < end;n++){ // For for each neighbor
      arma::Col<double> temp_X=X;
      temp_X(n)=1-temp_X(n);//Switch the nth coordinate
      double mid_step=loglik_internal(temp_X,M,theta)-logpi_current;
      output_bound[n]=mid_step*temperature;
      output[n]=apply_bal_func_bounded_internal(mid_step*temperature,log_bound,bal_func);

    }// End for loop
  }//End operator
};// End of worker to visit neighbors

struct GibbsSampler : public Worker
{
  //// Data
  RVector<double> X_n;//State matrix
  RVector<double> X_mode;//State matrix
  const double temperature;//Chosen temperature
  const double theta; //Parameter for likelihood
  RVector<double> output;
  const std::size_t dim_p; // dimension of the problem (and length of X columns)
  std::vector<std::mt19937_64>& rngs;//Vector of RNGs
  
  //// Functions to use from outside
  std::size_t flip_coord(std::size_t coord,bool approaching,double theta,std::mt19937_64& rng);
  
  // arma::Col<double> convert_X(){
  //   RVector<double> tmp_X = X_n;
  //   arma::Col<double> VEC(tmp_X.begin(), dim_p, false,true);
  //   return VEC;
  // }
  
  
  //// Main constructor
  GibbsSampler(
    NumericVector X_n_in,
    NumericVector X_mode_in,
    const double temperature_in,
    const double theta_in,
    NumericVector output_in,
    const std::size_t dim_p_in,
    std::vector<std::mt19937_64>& rngs_in)://Vector of RNGs
    X_n(X_n_in),
    X_mode(X_mode_in),
    temperature(temperature_in),
    theta(theta_in),
    output(output_in),
    dim_p(dim_p_in),
    rngs(rngs_in)
  {}
  
  
  
  //// Operator
  void operator()(std::size_t begin, std::size_t end){
    
    for(std::size_t n = begin ; n < end ; n++){ // For for each neighbor
      bool check_coord=X_n[n]!=X_mode[n];//Different coordinate means the swap will approach the mode
      output[n]=flip_coord(X_n[n],check_coord,theta*temperature,rngs[n]);
    }// End for loop
    
  }//End operator
};// End of worker to visit neighbors



struct SumExp_bf : public Worker//Not used, we do a double pass on the vector
{

  const RVector<double> X;// source vector
  const int bal_func; //Balancing function
  double Z; // accumulated value

  //Functions to use from outside
  double apply_bal_func_internal(double x,const int chosen);
  // constructors
  SumExp_bf(const NumericVector X,const int bal_func_in) : X(X),bal_func(bal_func_in), Z(0) {}
  SumExp_bf(const SumExp_bf& sum, Split) : X(sum.X),bal_func(bal_func), Z(0) {}

  // accumulate just the element of the range I've been asked to
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      Z += std::exp(apply_bal_func_internal(X[i],bal_func));
    }
  }// End of operator

  // join my value with that of another Sum
  void join(const SumExp_bf& rhs) {
    Z += rhs.Z;
  }
};//End of worker to compute the Z factor by adding log-probs

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
};//End of worker to compute the Z factor by adding log-probs

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
};//End of worker to find minimum of a vector

struct GetMax : public Worker
{   
  const RVector<double> X;// input vector
  double max_value; // MIN value found
  std::size_t max_index; // Index of MIN value
  // constructors
  GetMax(const NumericVector X) : X(X), max_value(-9999999),max_index(0) {}
  GetMax(const GetMax& max_cons, Split) : X(max_cons.X), max_value(-9999999),max_index(0) {}
  
  // accumulate just the element of the range I've been asked to
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      if(X[i]>max_value){
        max_value=X[i];
        max_index=i;
      }
    }
  }// End of operator
  
  // join my value with that of another Sum
  void join(const GetMax& other) {
    if(other.max_value>max_value){
      max_value = other.max_value;
      max_index = other.max_index;
    }
  }
};//End of worker to find maximum of a vector

struct GetMaxAbs : public Worker
{   
  const RVector<double> X;// input vector
  double max_value; // MIN value found
  std::size_t max_index; // Index of MIN value
  // constructors
  GetMaxAbs(const NumericVector X) : X(X), max_value(-9999999),max_index(0) {}
  GetMaxAbs(const GetMaxAbs& max_cons, Split) : X(max_cons.X), max_value(-9999999),max_index(0) {}
  
  // accumulate just the element of the range I've been asked to
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      double temp=abs(X[i]);
      if(temp>max_value){
        max_value=temp;
        max_index=i;
      }
    }
  }// End of operator
  
  // join my value with that of another Sum
  void join(const GetMaxAbs& other) {
    if(other.max_value>max_value){
      max_value = other.max_value;
      max_index = other.max_index;
    }
  }
};//End of worker to find maximum of a vector


#endif
