//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <fstream>
#include <ctime> 
#include <random>
#include <execution>
#include <mutex>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


//// To implement Random Number Generator with specific seed
// Thread-safe RNG wrapper
class ReproducibleRNG {
private:
  std::vector<std::mt19937> rngs; //Vector to store random number generators
  size_t seed; //Seed that user defines
  
public:
  // Function to create Random Number Generators, 1 for each thread
  ReproducibleRNG(size_t base_seed, size_t num_threads) : seed(base_seed) {
    std::mt19937 seeder(base_seed);
    for (size_t i = 0; i < num_threads; ++i) {
      rngs.emplace_back(seeder());
    }
  }
  
  std::mt19937& get_rng(size_t thread_id) {
    return rngs[thread_id % rngs.size()];
  }
};


////////// Other functions //////////
double ret_min(double a,double b,double c){
  vec temp={a,b,c};
  return(min(temp));
}


mat initializeRandom(const int& num_rows,const int& num_cols, const double& prob) {
  
  // Initialize a matrix with random values between 0 and 1
  arma::mat A = arma::randu<arma::mat>(num_rows, num_cols);
  
  // Threshold the random values to 0 or 1
  A = arma::conv_to<arma::mat>::from(A > prob);
  
  return(A);
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

double round_to(double value, int decimals){
  double mult=pow(10.0,decimals);
  value = round(value*mult)/mult;
  return value;
}


////////// Balancing functions //////////

///// List of balancing functions that apply to log-likelihoods
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

//// Bounded balancing function based on sqrt 
double bound_sq(double x, double log_gamma){
  double temp1=x/2-log_gamma;
  return ret_min(temp1,x,0);
}

///// Functions to call balancing functions inside other functions

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



// [[Rcpp::export]]
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

////////// Updating functions //////////

// Modified IIT_update_w for thread-safe operation
List IIT_update_w_parallel(vec X, const arma::mat& M, String chosen_bf, 
                           double temperature, double theta, 
                           std::mt19937& rng) {
  int total_neighbors = X.n_rows;
  vec probs(total_neighbors, fill::zeros);
  double logpi_current = loglik(X, M, theta);
  
  // Compute weights for all neighbors
  for(int j = 0; j < total_neighbors; j++) {
    vec newX = X;
    newX.row(j) = 1 - X.row(j);
    double temporal = loglik(newX, M, theta) - logpi_current;
    probs(j) = bal_func(temporal * temperature, chosen_bf);
  }
  
  // Thread-safe random number generation
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  vec u(total_neighbors);
  for(auto& val : u) val = dist(rng);
  
  vec probs_choose = log(-log(u)) - probs;
  int neigh_pos = (std::min_element(probs_choose.begin(), probs_choose.end())) - probs_choose.begin();
  
  X.row(neigh_pos) = 1 - X.row(neigh_pos);
  List result;
  result["X"] = X;
  result["Z"] = sum(exp(probs)) / total_neighbors;
  return result;
}

////////// Code to find temperatures for Parallel Tempering //////////

// [[Rcpp::export]]
List temperature_PT_IIT(int p,int interswap, double temp_ini,
                        const std::string bal_function,
                        const double& theta, int seed = 42, unsigned n_thread){
  // Inputs are:
  // p:dimension
  //interswap: number of iterations to try between replica swaps
  // t1: initial temperature 1 (not changes)
  // t2: initial temperature 2 (changes) t1/(1+exp(rho_n)) = t1/exp1m(rho_n)
  //threshold: Para detener la busqueda 
  // Numero de temperaturas (o temperatura mínima a alcanzar)
  //Theta: parametro a usar para la log-likelihood
  
  // Initialize RNG system
  if(n_thread==0){
    unsigned num_threads = std::thread::hardware_concurrency(); //Identify number of threads
  }else{
    unsigned num_threads = n_thread;
  }
    

  ReproducibleRNG rng_system(seed, num_threads);
  
  
  //// Initialize variables to use in the code
  int T=2; //The total number of temperatures will be 2
  mat X(p,T); // To store the current state of the joint chain, as many rows as neighbors, as many columns as temperatures
  List output; //To store output of the update function
  
  //// Parameters for the temperature finding   
  double rho=-3;//1.3;//-2.6; // To define initial second temperature
  double threshold=0.001;//.1;//0.001;//Stop when the updates differ less than this much
  double target_swap_rate=0.2345;//Target swap rate
  
  int precision=3;//To round values
  
  bool stay_in_loop=true;
  double swap_prob;
  vec Xtemp_from(p);
  vec Xtemp_to(p);
  //// Initialize temperatures
  vec temp={temp_ini,round_to(temp_ini/(1+exp(rho)),precision)};
  double avg_swap_prob=0;
  int count_convergence=0;
  
  double current_temp;
  int swap_count=0; //To count how many swaps were performed
  // int n_iter=0;
  //// Define two modes
  mat Q_matrix(p,2);
  for(int i=0;i<p;i++){
    if(i%2==0){Q_matrix(i,1)=1;}
    if(i%2==1){Q_matrix(i,0)=1;}
  }
  
  
  vec ppp=Rcpp::runif(1);
  X=initializeRandom(p,2,ppp(0));//Randomly initialize the state of each replica.
  std::clock_t start = std::clock(); /// Start timer
  while(stay_in_loop){
    //// Start the loop for all iterations
    for(int i=0;i<interswap;i++){
      
      for(int replica=0;replica<T;replica++){//For loop for replicas
        current_temp=temp(replica);// Extract temperature of the replica
        
        //// Update each replica independently
        // IIT_update_w(vec X,const arma::mat& M,String chosen_bf, double temperature, double theta)
        output=IIT_update_w(X.col(replica),Q_matrix,bal_function,current_temp, theta);
        
        X.col(replica)=vec(output(0)); //Update current state of the chain
      }//End loop to update replicas
    }// End loop of interswap
    
    //// Start replica swap process
    
    swap_count+=1;//Increase the count of swaps
    // n_iter+=1; //Increase denominator
    Xtemp_from=X.col(0);
    Xtemp_to=X.col(1);
    
    
    //// Computation of Z factors to correct bias
    double Z_fact_correc;//to temporarily store Z_factor correction
    double Z_temp11;
    double Z_temp12;
    double Z_temp21;
    double Z_temp22;
    
    output=IIT_update_w(Xtemp_from,Q_matrix,bal_function,temp(0),theta);
    Z_temp11=output(1);
    output=IIT_update_w(Xtemp_to,Q_matrix,bal_function,temp(1),theta);
    Z_temp22=output(1);
    output=IIT_update_w(Xtemp_from,Q_matrix,bal_function,temp(1),theta);
    Z_temp12=output(1);
    output=IIT_update_w(Xtemp_to,Q_matrix,bal_function,temp(0),theta);
    Z_temp21=output(1);
    
    Z_fact_correc=Z_temp12*Z_temp21/(Z_temp11*Z_temp22);
    
    //// Computing swap probability
    swap_prob=(temp(0)-temp(1))*(loglik(Xtemp_to,Q_matrix,theta) - loglik(Xtemp_from,Q_matrix,theta)); 
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
  std::clock_t end = std::clock(); // Stop timer
  // Calculate the time taken in seconds
  double duration = static_cast<double>(end - start) / CLOCKS_PER_SEC;
  Rcpp::Rcout <<"FINAL RESULTS:\nSwap: "<<swap_count<<" avg. swap prob: "<<avg_swap_prob <<" new temperature: "<< temp(1) << std::endl; 
  
  // return round_to(temp(1),3);
  List ret;
  ret["temp"]=round_to(temp(1),precision);
  ret["swap"]=swap_count;
  ret["swap_rate"]=round_to(avg_swap_prob,precision*2);
  ret["seconds"]=duration;
  return ret;
}



