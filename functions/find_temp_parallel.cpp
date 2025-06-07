// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;
using namespace arma;

#include <ctime> 
#include "workers.h"

/////##### FUNCTIONS #####/////

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
double round_to(double value, int decimals){
  double mult=pow(10.0,decimals);
  value = round(value*mult)/mult;
  return value;
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
    Rcpp::Rcout <<"M: " <<M<< std::endl;
    return(-10000);
  }
  
  
  
}

double loglik_R(NumericVector& X,const NumericMatrix& M, const double& theta){
  double loglik_computed;
  if(M.ncol()==2){
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
    
  }else{
    Rcpp::Rcout << "Error matrix number of columns isn't 2" << std::endl;
    Rcpp::Rcout <<"M: " <<M<< std::endl;
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


///// Full definition of internal functions of workers 
double IIT_visit_neighbors::apply_bal_func_internal(double x,const int chosen){
  return apply_bal_func(x,chosen);
}

double IIT_visit_neighbors::loglik_internal(const arma::Col<double>& X,const arma::Mat<double>& M, const double& theta){
  return loglik(X,M,theta);
}

// Function to run simulation

// [[Rcpp::export]]
List test_1(int p, double temperature, int bal_func, double theta){
  int T=1;
  
  const std::size_t begin = static_cast <size_t> (0);
  const std::size_t end = static_cast <size_t> (p); 
  
  //// Define two modes
  NumericMatrix Q_matrix(p,2);
  for(int i=0;i<p;i++){
    if(i%2==0){Q_matrix(i,1)=1;}
    if(i%2==1){Q_matrix(i,0)=1;}
  }
  
  double ppp=randu();
  arma::Mat<double> inter_mat(p,T);
  inter_mat=initializeRandom(p,T,ppp);//Randomly initialize the state of each replica.
  arma::Col<double> original_X = inter_mat.col(0);
  //Convert 
  NumericVector X=Rcpp::wrap(original_X);
  NumericVector output (p);
  // Declare constructor to visit all neighbors
  IIT_visit_neighbors iit_neighbors(X,
                                    Q_matrix,
                                    bal_func,
                                    temperature,
                                    theta,
                                    output,
                                    end);
  //Apply for to visit all neighbors
  parallelFor(begin,end,iit_neighbors);
  
  //Declare constructor to add log-probabilities
  SumExp sum_exp(output);
  //Get the sum of probabilities
  parallelReduce(0,end,sum_exp);
  
  //Substract minimum
  arma::Col<double> u_random(p,fill::randu);
  NumericVector u=Rcpp::wrap(u_random);
  NumericVector choose_min=log(-log(u)) - (output);
  //Declare constructor to find minimum
  GetMin get_min(choose_min);
  parallelReduce(0,end,get_min);
  
  output.push_back(sum_exp.Z);
  choose_min.push_back(get_min.min_value);
  choose_min.push_back(get_min.min_index);
  
  List ret;
  ret["output"] = output;
  ret["mins"]=choose_min;
  return ret;
  // return output;
  
}

// [[Rcpp::export]]
List temperature_PT_IIT(int p,int interswap, double temp_ini, int bal_func, const double& theta){
  // Inputs are:
  // p:dimension
  //interswap: number of iterations to try between replica swaps
  // t1: initial temperature 1 (not changes)
  // t2: initial temperature 2 (changes) t1/(1+exp(rho_n)) = t1/exp1m(rho_n)
  //threshold: Para detener la busqueda 
  // Numero de temperaturas (o temperatura mínima a alcanzar)
  //Theta: parametro a usar para la log-likelihood
  
  //// Initialize variables to use in the code
  int T=2; //The total number of temperatures will be 2
  NumericMatrix X(p,T); // To store the current state of the joint chain, as many rows as neighbors, as many columns as temperatures
  // List output; //To store output of the update function
  
  //// Parameters for the temperature finding   
  double rho=-3.5;//1.3;//-2.6; // To define initial second temperature
  double threshold=0.001;//0.001;//.1;//0.001;//Stop when the updates differ less than this much
  double target_swap_rate=0.2345;//Target swap rate
  
  int precision=3;//To round values
  
  bool stay_in_loop=true;
  double swap_prob;
  NumericVector Xtemp_from(p);
  NumericVector Xtemp_to(p);
  //// Initialize temperatures
  vec temp={temp_ini,round_to(temp_ini/(1+exp(rho)),precision)};
  double avg_swap_prob=0;
  int count_convergence=0;
  
  double current_temp;
  int swap_count=0; //To count how many swaps were performed
  
  //// Define two modes
  NumericMatrix Q_matrix(p,2);
  for(int i=0;i<p;i++){
    if(i%2==0){Q_matrix(i,1)=1;}
    if(i%2==1){Q_matrix(i,0)=1;}
  }
  
  //// Initialize states X_0
  double ppp=randu();
  arma::Mat<double> inter_mat(p,T);
  inter_mat=initializeRandom(p,T,ppp);//Randomly initialize the state of each replica.
  X=Rcpp::wrap(inter_mat);
  
  //// Paramemters for parallelization
  // const std::size_t begin = static_cast <size_t> (0);
  const std::size_t end = static_cast <size_t> (p); 
  
  std::clock_t start_time = std::clock(); /// Start timer

  while(stay_in_loop){
    //// Start the loop for all iterations
    for(int i=0;i<interswap;i++){
      for(int replica=0;replica<T;replica++){//For loop for replicas
        current_temp=temp(replica);// Extract temperature of the replica
//// Visit all neighbors
        NumericVector output_bf(p); //Declare vector to store info of visiting neighbors
//// Declare constructor to visit all neighbors
        IIT_visit_neighbors iit_neighbors(X(_,replica),
                                          Q_matrix,
                                          bal_func,
                                          current_temp,
                                          theta,
                                          output_bf,
                                          end);

        parallelFor(0,end,iit_neighbors);//Apply ParallelFor
        
        //Declare constructor to add log-probabilities
        // SumExp sum_exp(output_bf);
        //Get the sum of probabilities
        // parallelReduce(0,end,sum_exp);
        
//// Choose the next neighbor state

        arma::Col<double> u_random(p,fill::randu);
        NumericVector u=Rcpp::wrap(u_random);
        NumericVector choose_min=log(-log(u)) - (output_bf);
        //Declare constructor to find minimum
        GetMin get_min(choose_min);
        parallelReduce(0,end,get_min);//Find the Minimum
//// Update the state of the replica
double temp_coord=X(get_min.min_index,replica);
        X(get_min.min_index,replica)=1-temp_coord;
      }//End loop to update replicas
    }// End loop of interswap
    
    //// Start replica swap process
    
    swap_count+=1;//Increase the count of swaps
    // n_iter+=1; //Increase denominator
    Xtemp_from=X(_,0);
    Xtemp_to=X(_,1);
    
    
    //// Computation of Z factors to correct bias
    double Z_fact_correc;//to temporarily store Z_factor correction
    // double Z_temp11;
    // double Z_temp12;
    // double Z_temp21;
    // double Z_temp22;
    
    //Declare vector to store info of visiting neighbors
    NumericVector output_Z_v11(p); 
    NumericVector output_Z_v22(p); 
    NumericVector output_Z_v12(p); 
    NumericVector output_Z_v21(p); 
//// Declare constructor to visit all neighbors
    IIT_visit_neighbors Z_v11(Xtemp_from,
                              Q_matrix,
                              bal_func,
                              temp(0),
                              theta,
                              output_Z_v11,
                              end);
    IIT_visit_neighbors Z_v22(Xtemp_to,
                              Q_matrix,
                              bal_func,
                              temp(1),
                              theta,
                              output_Z_v22,
                              end);
    IIT_visit_neighbors Z_v12(Xtemp_from,
                              Q_matrix,
                              bal_func,
                              temp(1),
                              theta,
                              output_Z_v12,
                              end);
    IIT_visit_neighbors Z_v21(Xtemp_to,
                              Q_matrix,
                              bal_func,
                              temp(0),
                              theta,
                              output_Z_v21,
                              end);
//// Apply ParallelFor
    parallelFor(0,end,Z_v11);
    parallelFor(0,end,Z_v22);
    parallelFor(0,end,Z_v12);
    parallelFor(0,end,Z_v21);
    
//// Declare constructor to add log-probabilities
    SumExp sum_Z11(output_Z_v11);
    SumExp sum_Z22(output_Z_v22);
    SumExp sum_Z12(output_Z_v12);
    SumExp sum_Z21(output_Z_v21);
//// Get the sum of probabilities
    parallelReduce(0,end,sum_Z11);
    parallelReduce(0,end,sum_Z22);
    parallelReduce(0,end,sum_Z12);
    parallelReduce(0,end,sum_Z21);
    
    Z_fact_correc=(sum_Z12.Z*sum_Z21.Z)/(sum_Z11.Z*sum_Z22.Z);
    // Z_fact_correc=Z_temp12*Z_temp21/(Z_temp11*Z_temp22);
    
    //// Computing swap probability
    swap_prob=(temp(0)-temp(1))*(loglik_R(Xtemp_to,Q_matrix,theta) - loglik_R(Xtemp_from,Q_matrix,theta)); 
    swap_prob=ret_min(Z_fact_correc*exp(swap_prob),1,1);
    
    
    if(swap_count==1){
      avg_swap_prob=swap_prob;
    }else{
      avg_swap_prob=(avg_swap_prob*(swap_count-1)+swap_prob)/swap_count;
    }
    
    //// Update temperature
    rho=rho + (swap_prob-target_swap_rate)/swap_count;
    
    // if(rho<1e-5){
    //   temp(1)=round_to(temp_ini/(2+expm1(rho)),precision);
    // }else{
    //   temp(1)=round_to(temp_ini/(1+exp(rho)),precision);
    // }
    if(rho<1e-5){
      temp(1)=temp_ini/(2+expm1(rho));
    }else{
      temp(1)=temp_ini/(1+exp(rho));
    }
    
    //Check current threshold
    // Rcpp::Rcout <<"Avg. swap_rate: "<<avg_swap_prob << std::endl; 
    if((target_swap_rate-avg_swap_prob)<threshold && (target_swap_rate-avg_swap_prob)>-threshold){
      count_convergence+=1;
      if(count_convergence>=3){
        stay_in_loop=false;
      }
      
    }else{
      count_convergence=0;
    }
    
    
    if(swap_count % 100 == 0){
      Rcpp::Rcout <<"Swap: "<<swap_count<<" avg. swap prob: "<<avg_swap_prob <<" new temperature: "<< temp(1) << std::endl; 
    }
    
    if(swap_count == 500000){// Force finishing of algorithm
      stay_in_loop=false;
    } 
    
  }
  std::clock_t end_time = std::clock(); // Stop timer
  // Calculate the time taken in seconds
  double duration = static_cast<double>(end_time - start_time) / CLOCKS_PER_SEC;
  Rcpp::Rcout <<"FINAL RESULTS:\nSwap: "<<swap_count<<" avg. swap prob: "<<avg_swap_prob <<" new temperature: "<< temp(1) << std::endl; 
  
  // return round_to(temp(1),3);
  List ret;
  ret["temp"]=round_to(temp(1),precision);
  ret["swap"]=swap_count;
  ret["swap_rate"]=round_to(avg_swap_prob,precision*2);
  ret["seconds"]=duration;
  return ret;
}





