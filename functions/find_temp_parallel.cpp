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

std::vector<std::mt19937_64> initialize_rngs(int n, int base_seed) {
  std::vector<std::mt19937_64> rngs(n);
  for(int i = 0; i < n; i++) {
    rngs[i].seed(base_seed + i);  // Unique seed for each RNG
  }
  return rngs;
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

// Full definition to update coord
std::size_t GibbsSampler::flip_coord(std::size_t coord,bool approaching,double theta,std::mt19937_64& rng){
  std::size_t new_coord=coord;//To store new coordinate
  double u;//To store random number
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  u = dist(rng);
  if(approaching){//In case the proposed coordinate is approaching the mode
    if(u<(1/(1+exp(-theta)))){new_coord=(1-coord);}
    // Rcpp::Rcout <<"Accept approach"<< std::endl; 
  }else{//In case the proposed coordinate is moving away from the mode
    if(u<(1/(1+exp(theta)))){new_coord=(1-coord);}
    // Rcpp::Rcout <<"Accept move far"<< std::endl; 
  }
  return new_coord;
}

// Function to run simulation

// [[Rcpp::export]]
List find_temp_gibbs_A_IIT(int p,int interswap, int burn_in,double temp_ini, int bal_func, const double& theta, int base_seed, int direction){
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
  double threshold=0.0005;//.1;//0.001;//Stop when the updates differ less than this much
  double target_swap_rate=0.2345;//Target swap rate
  
  int precision=3;//To round values
  
  bool stay_in_loop=true;
  double swap_prob;
  NumericVector Xtemp_from(p);
  NumericVector Xtemp_to(p);
  //// Initialize temperatures
  double temp_ini_2;
  if(direction>0){
    temp_ini_2=round_to(temp_ini*(1+exp(rho)),precision)
  }else{
    temp_ini_2=round_to(temp_ini/(1+exp(rho)),precision)
  }
  vec temp={temp_ini,temp_ini_2};
  double avg_swap_prob=0;
  int count_convergence=0;
  
  double current_temp;
  int swap_count=0; //To count how many swaps were performed
  
  //// Define two modes
  NumericMatrix Q_matrix(p,2);
  for(int i=0;i<p;i++){
    if(i%2==0){Q_matrix(i,1)=1;
      X(i,1)=1;
    }
    if(i%2==1){Q_matrix(i,0)=1;
      X(i,0)=1;
    }
  }
  
  //// Initialize states X_0
  // double ppp=randu();
  // arma::Mat<double> inter_mat(p,T);
  // inter_mat=initializeRandom(p,T,ppp);//Randomly initialize the state of each replica.
  // X=Rcpp::wrap(inter_mat);
  
  //// Paramemters for parallelization
  // const std::size_t begin = static_cast <size_t> (0);
  const std::size_t end = static_cast <size_t> (p); 
  NumericMatrix chain_m1=clone(X); // To keep track of how the chain for mode 1 evolves
  NumericMatrix chain_m2=clone(X);// To keep track of how the chain for mode 2 evolves
  
  NumericVector X_mode1=Q_matrix(_,0);
  NumericVector X_mode2=Q_matrix(_,1);
  
  std::vector<std::mt19937_64> rngs = initialize_rngs(p,base_seed);//Initialize RNGs
  
  
  //Burn-in period
  for(int i=0;i<burn_in;i++){
    for(int replica=0;replica<T;replica++){//For loop for replicas
      current_temp=temp(replica);// Extract temperature of the replica
      
      
      NumericVector output_m1(p); //Declare vector to store info of first chain
      //// Declare constructor to update chain 1
      GibbsSampler g_sample_m1(chain_m1(_,replica),
                               X_mode1,
                               current_temp,
                               theta,
                               output_m1,
                               end,
                               rngs);
      parallelFor(0,end,g_sample_m1);//Apply ParallelFor
      
      NumericVector output_m2(p); //Declare vector to store info of 2nd chain
      //// Declare constructor to update chain 2
      GibbsSampler g_sample_m2(chain_m2(_,replica),
                               X_mode2,
                               current_temp,
                               theta,
                               output_m2,
                               end,
                               rngs);
      parallelFor(0,end,g_sample_m2);//Apply ParallelFor
      
      //Update chains
      chain_m1(_,replica)=output_m1;
      chain_m2(_,replica)=output_m2;
      //Update state of the chain based on the mixture with the same weight
    }//End loop to update replicas
  }// End loop of burn-in

  double ppp1=randu();
  if(ppp1<0.5){X=clone(chain_m1);}else{X=clone(chain_m2);}  
  
  std::clock_t start_time = std::clock(); /// Start timer  
  //Start finding temperatures
  while(stay_in_loop){
    //// Start the loop for all iterations
    for(int i=0;i<interswap;i++){
      for(int replica=0;replica<T;replica++){//For loop for replicas
        current_temp=temp(replica);// Extract temperature of the replica
        
        
        NumericVector output_m1(p); //Declare vector to store info of first chain
        //// Declare constructor to update chain 1
        GibbsSampler g_sample_m1(chain_m1(_,replica),
                                 X_mode1,
                                 current_temp,
                                 theta,
                                 output_m1,
                                 end,
                                 rngs);
        parallelFor(0,end,g_sample_m1);//Apply ParallelFor
        
        NumericVector output_m2(p); //Declare vector to store info of 2nd chain
        //// Declare constructor to update chain 2
        GibbsSampler g_sample_m2(chain_m2(_,replica),
                                 X_mode2,
                                 current_temp,
                                 theta,
                                 output_m2,
                                 end,
                                 rngs);
        parallelFor(0,end,g_sample_m2);//Apply ParallelFor
        
        //Update chains
        chain_m1(_,replica)=output_m1;
        chain_m2(_,replica)=output_m2;
        //Update state of the chain based on the mixture with the same weight
        double ppp1=randu();
        if(ppp1<0.5){X(_,replica)=output_m1;}else{X(_,replica)=output_m2;}
        
      }//End loop to update replicas
    }// End loop of interswap
    
    //// Start replica swap process
    
    swap_count+=1;//Increase the count of swaps
    // n_iter+=1; //Increase denominator
    Xtemp_from=X(_,0);
    Xtemp_to=X(_,1);
    
    
    //// Computing swap probability
    swap_prob=(temp(0)-temp(1))*(loglik_R(Xtemp_to,Q_matrix,theta) - loglik_R(Xtemp_from,Q_matrix,theta)); 
    swap_prob=ret_min(exp(swap_prob),1,1);
    
    
    if(swap_count==1){
      avg_swap_prob=swap_prob;
    }else{
      avg_swap_prob=(avg_swap_prob*(swap_count-1)+swap_prob)/swap_count;
    }
    
    //// Update temperature
    rho=rho + (swap_prob-target_swap_rate)/swap_count;
    
    double temp_update;
    if(direction>0){
      if(rho<1e-5){
        temp_update=temp_ini*(2+expm1(rho));
      }else{
        temp_update=temp_ini*(1+exp(rho));
      }
    }else{
      if(rho<1e-5){
        temp_update=temp_ini/(2+expm1(rho));
      }else{
        temp_update=temp_ini/(1+exp(rho));
      }
    }
    
    temp(1)=temp_update;
    
    
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
    
    
    if(swap_count % 500 == 0){
      Rcpp::Rcout <<"A-IIT ("<<base_seed<<") p: "<<p<<" Swap: "<<swap_count<<" avg. swap prob: "<<avg_swap_prob <<" new temperature: "<< temp(1) << std::endl; 
      
    }
    
    if(swap_count == 1000000){// Force finishing of algorithm
      stay_in_loop=false;
    } 
    
  }
  std::clock_t end_time = std::clock(); // Stop timer
  // Calculate the time taken in seconds
  double duration = static_cast<double>(end_time - start_time) / CLOCKS_PER_SEC;
  Rcpp::Rcout <<"FINAL RESULTS:\nSwap: "<<swap_count<<" avg. swap prob: "<<avg_swap_prob <<" new temperature: "<< temp(1) << std::endl; 
  
  // return round_to(temp(1),3);
  List ret;
  ret["temp"]=round_to(temp(1),precision*2);
  ret["swap"]=swap_count;
  ret["swap_rate"]=round_to(avg_swap_prob,precision*2);
  ret["seconds"]=duration;
  return ret;
}

// [[Rcpp::export]]
List find_temp_gibbs_PT_IIT(int p, int burn_in,double temp_ini, int bal_func, const double& theta, int gibbs_steps, int direction){
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
  double threshold=0.0005;//0.001;//.1;//0.001;//Stop when the updates differ less than this much
  double target_swap_rate=0.2345;//Target swap rate
  
  int precision=3;//To round values
  
  bool stay_in_loop=true;
  double swap_prob;
  NumericVector Xtemp_from(p);
  NumericVector Xtemp_to(p);
  //// Initialize temperatures
  double temp_ini_2;
  if(direction>0){
    temp_ini_2=round_to(temp_ini*(1+exp(rho)),precision)
  }else{
    temp_ini_2=round_to(temp_ini/(1+exp(rho)),precision)
  }
  vec temp={temp_ini,temp_ini_2};
  double avg_swap_prob=0;
  int count_convergence=0;
  
  double current_temp;
  int swap_count=0; //To count how many swaps were performed
  
  //// Define two modes and Initialize states X_0
  //We initialize X_0 at the modes
  NumericMatrix Q_matrix(p,2);
  for(int i=0;i<p;i++){
    if(i%2==0){Q_matrix(i,1)=1;
      X(i,1)=1;
      }
    if(i%2==1){Q_matrix(i,0)=1;
      X(i,0)=1;
      }
  }
  //// Initialize states X_0
  // double ppp=randu();
  // arma::Mat<double> inter_mat(p,T);
  // inter_mat=initializeRandom(p,T,ppp);//Randomly initialize the state of each replica.
  // X=Rcpp::wrap(inter_mat);
  

  
  //// Paramemters for parallelization
  const std::size_t end = static_cast <size_t> (p); 
  NumericVector X_mode1=Q_matrix(_,0);
  NumericVector X_mode2=Q_matrix(_,1);
  
   if(gibbs_steps>p|gibbs_steps<=0){gibbs_steps=p/2;}//In case the Gibbs step size is too big
  //Burn-in period
  // Rcpp::Rcout <<"Starting state: \n"<<X<< std::endl; 
   
  for(int i=0;i<burn_in;i++){
      Rcpp::Rcout <<"Burn_in, iteration: "<<i<< std::endl; 
    for(int replica=0;replica<T;replica++){//For loop for replicas
      current_temp=temp(replica);// Extract temperature of the replica
      // First gather information of current states of both chains
      NumericVector current_X=X(_,replica);//Copy of current state of chain 1
      
      NumericVector output_X(p);//Vector to store output
      IIT_visit_neighbors visit_current_X(current_X,
                                             Q_matrix,
                                             bal_func,
                                             current_temp,
                                             theta,
                                             output_X,
                                             end);
      
      parallelFor(0,end,visit_current_X);//Apply parallelFor
      SumExp info_X(output_X);// Declare constructor to add log-probabilities
      parallelReduce(0,end,info_X);// Get the sum of probabilities
      
      double dist_m1=sum(abs(current_X - X_mode1));
      double dist_m2=sum(abs(current_X - X_mode2));
      double density_current_X = (info_X.Z*exp(current_temp*loglik_R(current_X,Q_matrix,theta)));
      
      // Start for loop for coordinates        
      for(int coord=0;coord<p;coord++){
        /// Updating chain
        NumericVector output_X_neighbor(p); 
        NumericVector X_neighbor=clone(current_X);//Copy the current state of chain 1
        X_neighbor[coord]=1-X_neighbor[coord];//Flip coordinate
        //// Declare constructor to visit all neighbors
        IIT_visit_neighbors visit_X_neighbors(X_neighbor,
                                                    Q_matrix,
                                                    bal_func,
                                                    current_temp,
                                                    theta,
                                                    output_X_neighbor,
                                                    end);
        
        parallelFor(0,end,visit_X_neighbors);//Apply parallelFor
        SumExp info_X_neighbor(output_X_neighbor);/// Declare constructor to add log-probabilities
        parallelReduce(0,end,info_X_neighbor);//// Get the sum of probabilities
        
        double prob_flip;
        
// Check if we flip coordinate
//We check if the coordinate is different from mode 1
// bool check_coord_m1 = (current_X[coord]!=X_mode1[coord]);
// Compute the density of the proposed swap (proposed state)
//With this probability we flip the coordinate
          double density_X_neighbor = (info_X_neighbor.Z*exp(current_temp*loglik_R(X_neighbor,Q_matrix,theta)));
          prob_flip=(density_X_neighbor)/(density_X_neighbor + density_current_X);
        double ppm1=randu();
        if(ppm1<prob_flip){X(coord,replica) = 1-X(coord,replica);}
        
      }//End loop for coordinate 
    }//End loop to update replicas
  }// End loop of interswap
  Rcpp::Rcout <<"Finish Burn_in: "<< std::endl;   
  // Rcpp::Rcout <<"State after burn-in: \n"<<X<< std::endl;
  std::clock_t start_time = std::clock(); /// Start timer  
  //Start finding temperatures
  while(stay_in_loop){
    for(int replica=0;replica<T;replica++){//For loop for replicas
      current_temp=temp(replica);// Extract temperature of the replica
      // First gather information of current states of both chains
      NumericVector current_X=X(_,replica);//Copy of current state of chain 1
      
      NumericVector output_X(p);//Vector to store output
      IIT_visit_neighbors visit_current_X(current_X,
                                          Q_matrix,
                                          bal_func,
                                          current_temp,
                                          theta,
                                          output_X,
                                          end);
      
      parallelFor(0,end,visit_current_X);//Apply parallelFor
      SumExp info_X(output_X);// Declare constructor to add log-probabilities
      parallelReduce(0,end,info_X);// Get the sum of probabilities
      
      // double dist_m1=sum(abs(current_X - X_mode1));
      // double dist_m2=sum(abs(current_X - X_mode2));
      double density_current_X = (info_X.Z*exp(current_temp*loglik_R(current_X,Q_matrix,theta)));// Start for loop for coordinates        
      for(int coord=0;coord<p;coord++){
        /// Updating chain 1
        NumericVector output_X_neighbor(p); 
        NumericVector X_neighbor=clone(current_X);//Copy the current state of chain 1
        X_neighbor[coord]=1-X_neighbor[coord];//Flip coordinate
        //// Declare constructor to visit all neighbors
        IIT_visit_neighbors visit_X_neighbors(X_neighbor,
                                              Q_matrix,
                                              bal_func,
                                              current_temp,
                                              theta,
                                              output_X_neighbor,
                                              end);
        
        parallelFor(0,end,visit_X_neighbors);//Apply parallelFor
        SumExp info_X_neighbor(output_X_neighbor);/// Declare constructor to add log-probabilities
        parallelReduce(0,end,info_X_neighbor);//// Get the sum of probabilities
        
        double prob_flip;
          // Check if we flip coordinate
          //We check if the coordinate is different from mode 1
          // bool check_coord_m1 = (current_X[coord]!=X_mode1[coord]);

            double density_X_neighbor = (info_X_neighbor.Z*exp(current_temp*loglik_R(X_neighbor,Q_matrix,theta)));
            prob_flip=(density_X_neighbor)/(density_X_neighbor + density_current_X);

          double ppm1=randu();
          if(ppm1<prob_flip){X(coord,replica) = 1-X(coord,replica);}

//// Try a swap after specific number of steps
        if((swap_count+coord) % gibbs_steps == 0){

          //// Start replica swap process
          swap_count+=1;//Increase the count of swaps
          // n_iter+=1; //Increase denominator
          Xtemp_from=X(_,0);
          Xtemp_to=X(_,1);
          
          // Rcpp::Rcout <<"Trying swap: "<<swap_count<< std::endl;   
          //// Computation of Z factors to correct bias
          double Z_fact_correc=1;//to temporarily store Z_factor correction
          
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
          // Rcpp::Rcout <<"Finish parallel steps swap: "<<swap_count<< std::endl; 
          Z_fact_correc=(sum_Z12.Z*sum_Z21.Z)/(sum_Z11.Z*sum_Z22.Z);
          //// Computing swap probability
          swap_prob=(temp(0)-temp(1))*(loglik_R(Xtemp_to,Q_matrix,theta) - loglik_R(Xtemp_from,Q_matrix,theta)); 
          swap_prob=ret_min(Z_fact_correc*exp(swap_prob),1,1);
          
          if(swap_count==1){avg_swap_prob=swap_prob;
          }else{avg_swap_prob=(avg_swap_prob*(swap_count-1)+swap_prob)/swap_count;}    
          //// Update temperature
          rho=rho + (swap_prob-target_swap_rate)/swap_count;
          double temp_update;
          if(direction>0){
            if(rho<1e-5){
              temp_update=temp_ini*(2+expm1(rho));
            }else{
              temp_update=temp_ini*(1+exp(rho));
            }
          }else{
            if(rho<1e-5){
              temp_update=temp_ini/(2+expm1(rho));
            }else{
              temp_update=temp_ini/(1+exp(rho));
            }
          }
          temp(1)=temp_update;
          //Check current threshold 
          if((target_swap_rate-avg_swap_prob)<threshold && (target_swap_rate-avg_swap_prob)>-threshold){
            count_convergence+=1;
            if(count_convergence>=3){
              stay_in_loop=false;}
          }else{count_convergence=0;}
          if(swap_count % 200 == 0){
            Rcpp::Rcout <<"PT-IIT p: "<<p<<" Swap: "<<swap_count<<" avg. swap prob: "<<avg_swap_prob <<" new temperature: "<< temp(1) << std::endl; 
            // Rcpp::Rcout <<"Current state: \n"<<X<< std::endl; 
          }
          if(swap_count == 1000000){// Force finishing of algorithm
            stay_in_loop=false;} 
          /// Finish replica swap process   
          
        }//End IF for swap prob checking
      }//End loop for coordinates
      // Rcpp::Rcout <<"Finish 1 replica, swap: "<<swap_count<< std::endl; 
    }//End loop to update replicas
    // Rcpp::Rcout <<"Finish both replicas, swap: "<<swap_count<< std::endl;       
  }// End while loop
  std::clock_t end_time = std::clock(); // Stop timer
  // Calculate the time taken in seconds
  double duration = static_cast<double>(end_time - start_time) / CLOCKS_PER_SEC;
  Rcpp::Rcout <<"FINAL RESULTS:\nSwap: "<<swap_count<<" avg. swap prob: "<<avg_swap_prob <<" new temperature: "<< temp(1) << std::endl; 
  
  // return round_to(temp(1),3);
  List ret;
  ret["temp"]=round_to(temp(1),precision*2);
  ret["swap"]=swap_count;
  ret["swap_rate"]=round_to(avg_swap_prob,precision*2);
  ret["seconds"]=duration;
  return ret;
}



