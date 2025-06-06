#include <RcppArmadillo.h>
#include <execution>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
////////////To compare parellelization ////////////
///// List of balancing functions that apply to log-likelihoods
double ret_min(double a,double b,double c){
  vec temp={a,b,c};
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

double loglik(const arma::vec& X,const arma::mat& M, const double& theta){
  double loglik_computed;
  
  
  if(M.n_cols==2){
    //Define vectors to represent each mode
    vec mod1=M.col(0);
    vec mod2=M.col(1);
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
  
  // if(M.n_cols==3){
  //   //Define vectors to represent each mode
  //   vec mod1=M.col(0);
  //   vec mod2=M.col(1);
  //   vec mod3=M.col(2);
  //   // Compute distances to each mode  
  //   double dif1=L1_distance(X,mod1);
  //   double dif2=L1_distance(X,mod2);
  //   double dif3=L1_distance(X,mod3);
  // }  
  
  
}

mat initializeRandom(const int& num_rows,const int& num_cols, const double& prob) {
  
  // Initialize a matrix with random values between 0 and 1
  arma::mat A = arma::randu<arma::mat>(num_rows, num_cols);
  
  // Threshold the random values to 0 or 1
  A = arma::conv_to<arma::mat>::from(A > prob);
  
  return(A);
}


List IIT_update_w(vec X,const arma::mat& M,String chosen_bf,
                  double temperature, double theta){
  int total_neighbors = X.n_rows; // total number of neighbors is p spacial
  vec probs(total_neighbors, fill::zeros); //probabilities
  
  ////// Compute likelihood of current state
  double logpi_current=0;
  // uvec current_coord = find(X==1);
  logpi_current = loglik(X,M,theta);
  ////// Finish computing likelihood of current state
  
  // Rcpp::Rcout << "current loglik: "<< logpi_current << std::endl;
  ////// Compute weight for all neighbors
  double temporal=0;
  vec newX;
  for(int j=0; j<total_neighbors;j++){
    // Rcpp::Rcout << "Starts checking neighbors  "<< j<<std::endl; 
    newX = X;
    newX.row(j) = 1-X.row(j);
    //Rcpp::Rcout << newX << std::endl;
    // uvec coord = find(newX==1);
    //Rcpp::Rcout << coord << std::endl;
    temporal=loglik(newX,M,theta)-logpi_current;
    //Apply balancing function to log probability times temperature ////
    probs(j)=bal_func(temporal*temperature, chosen_bf);
  }
  // Thread-safe random numbers for selection
  vec u(total_neighbors,fill::randu);
  
  vec probs_choose = log(-log(u)) - probs;
  
  //Find the index of the minimum element. 
  //This corresponds to choosing that neighbor
  int neigh_pos = (std::min_element(probs_choose.begin(), probs_choose.end()))-probs_choose.begin();
  // Rcpp::Rcout <<"probs vector: "<< probs << std::endl;
  // Rcpp::Rcout <<"chosen neighbor: "<< neigh_pos << std::endl;
  
  X.row(neigh_pos) = 1-X.row(neigh_pos); //modify the coordinate of the chosen neighbor
  List result;
  result["X"]=X; // Return new state
  result["Z"]=sum(exp(probs))/total_neighbors; //Compute Z factor
  return result;
}

// [[Rcpp::export]]
mat run_sim(int p, const int iter,const vec temp,const std::string bal_function,
             double theta=3,int seed=5){
  int T=temp.n_rows;
  double current_temp;
  vec current_col;
  List output;
  //Define seed
  arma_rng::set_seed(seed);
  
  //// Define two modes
  mat Q_matrix(p,T);
  for(int i=0;i<p;i++){
    if(i%2==0){Q_matrix(i,1)=1;}
    if(i%2==1){Q_matrix(i,0)=1;}
  }
  
  double ppp=randu();
  Rcpp::Rcout <<"random number: " <<ppp<< std::endl;
  mat X(p,T);
  X=initializeRandom(p,T,ppp);//Randomly initialize the state of each replica.
  vec col_indices(T);
  std::iota(col_indices.begin(), col_indices.end(), 0);
  
  // For loop to update
  
  for(int i=0;i<iter;i++){
    
    // std::for_each(std::execution::par,
    std::for_each(
                  col_indices.begin(),
                  col_indices.end(),
                  [&](int j){

      current_temp=temp(j);
      current_col=X.col(j);
      output=IIT_update_w(current_col,Q_matrix,bal_function,current_temp,theta);
      X.col(j) = as<arma::vec>(output["X"]);
                  });//End of for each
    
    // for(int replica=0;replica<T;replica++){//For loop for replicas
    //   current_temp=temp(replica);// Extract temperature of the replica
    //   output=IIT_update_w(X.col(replica),Q_matrix,bal_function,current_temp, theta);
    //   X.col(replica)=vec(output(0)); //Update current state of the chain
    // }//End loop to update replicas
    
  }// End loop of interswap
  

  
  
  
  return X;
}//End of function





