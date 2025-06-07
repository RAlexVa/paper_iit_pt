//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <fstream>
#include <ctime> 
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


////////// Other functions //////////

double ret_min(double a,double b,double c){
  vec temp={a,b,c};
  return(min(temp));
}


// [[Rcpp::export]]
mat initializeRandom(const int& num_rows,const int& num_cols, const double& prob) {
  
  // Initialize a matrix with random values between 0 and 1
  arma::mat A = arma::randu<arma::mat>(num_rows, num_cols);
  
  // Threshold the random values to 0 or 1
  A = arma::conv_to<arma::mat>::from(A > prob);
  
  return(A);
}


// [[Rcpp::export]]
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
// [[Rcpp::export]]
double bf_sq(double x){
  return x/2;
}

// [[Rcpp::export]]
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
// [[Rcpp::export]]
double bound_sq(double x, double log_gamma){
  double temp1=x/2-log_gamma;
  return ret_min(temp1,x,0);
}

///// Functions to call balancing functions inside other functions

//This function we don't export to R because it generates errors
double invoke(double x, double (*func)(double)) {
  return func(x);
}

// [[Rcpp::export]]
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

// [[Rcpp::export]]
List IIT_update_w(vec X,const arma::mat& M,String chosen_bf, double temperature, double theta){
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
  //////Choose the next neighbor
  vec u = Rcpp::runif(total_neighbors);
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
List a_IIT_update(vec X,const arma::mat& M, String chosen_bf, const double& temperature, double log_bound, const bool& update, double prob_to_dec, const double& decreasing_constant, double max_logbound_found, double theta){
  int total_neighbors = X.n_rows; // total number of neighbors is p spacial
  double logpi_current= loglik(X,M,theta);// likelihood of current state
  
  const double threshold = 1e-5;//Threshold for updating maximum bound. It has to change at least this
  vec logprobs(total_neighbors, fill::zeros); //probabilities
  vec max_logprobs(total_neighbors,fill::zeros);//vector to store max-log-probabilities
  double temporal=0;
  vec newX;
  

  if(update){//If it's defined to keep updating the bounding constant
    ////// Compute weight for all neighbors
    for(int j=0; j<total_neighbors;j++){
      newX = X;
      newX.row(j) = 1-X.row(j);//Change coordinate of the state to define new neighbor
      temporal= temperature*(loglik(newX,M,theta)-logpi_current);
      logprobs(j)=temporal;//Store raw log_probability
      max_logprobs(j)=abs(temporal); // Store the max log-probability, either pi_y/pi_x or pi_x/pi_y
    }// End of loop to compute raw log-probability of neighbors
    
    double checking_max_logprob=bal_func(max(max_logprobs),"sq");
    if(max_logbound_found<checking_max_logprob && (checking_max_logprob-max_logbound_found)>=threshold){
      // Rcpp::Rcout <<"Diff en log bound: "<< checking_max_logprob-max_logbound_found<<std::endl;
      // Rcpp::Rcout <<"Previous max log bound: "<< max_logbound_found*1000000<<std::endl;
      max_logbound_found=checking_max_logprob;
      // Rcpp::Rcout <<"New max log bound: "<< max_logbound_found*1000000<<std::endl;
      log_bound=max_logbound_found;//First update according to the maximum 
    }
    //Then try Arithmetic reduction of the bounding constant
    if(prob_to_dec>0){//If we consider a probability to decrease the constant
      double test_prob=0;
      if(prob_to_dec<1){
        vec ppp=Rcpp::runif(1);//Get a random number
        test_prob=ppp(0);
      }
      if(test_prob<prob_to_dec){//In case the update is accepted
        double temporal_log_b=log_bound/temperature;
        double delta_bound = decreasing_constant/exp(temporal_log_b);
        //log(a-b)=log(a)+log(1-b/a) ≈ log(a) - b/a if b/a is very small
        //a=exp(temporal_log_b), b=decreasing_constant
        if(delta_bound>=1){//If the decreasing constant is too big
          //log(1-b/a) would be undefined
          log_bound=0;//Go the minimum 
        }else{
          log_bound = temperature*(temporal_log_b + log1p(-delta_bound));
        }
        
        // Rcpp::Rcout <<"New log_b: "<< log_bound<<std::endl;
        // if(temperature==0.05){
        //   Rcpp::Rcout <<"Decreasing delta: "<< delta_bound<<" Current bound: "<<exp(log_bound)<<" new bound: "<<log_bound<<std::endl;
        //   Rcpp::Rcout <<"Decreasing log-bound to "<< log_bound<<std::endl;
        // }
        if(log_bound<0){log_bound=0;} //Minimum bound is 1, minimum log_bound is 0
        
      }else{
        // Rcpp::Rcout <<"Rejected a bounding constant decrease"<< std::endl;
      }
    }
    //// Apply bounding B.F. with updated constant    
    for(int j=0;j<total_neighbors;j++){
      logprobs(j)=bound_sq(logprobs(j),log_bound);
    }
  }else{//In case we don't update the bouding constant
    
    //// Apply bounded B.F. without modifying the constant
    for(int j=0;j<total_neighbors;j++){
      newX = X;
      newX.row(j) = 1-X.row(j); //Change coordinate of the state to define new neighbor
      temporal = temperature*(loglik(newX,M,theta)-logpi_current);
      logprobs(j)=bound_sq(temporal,log_bound);
    }
  }
  
  /////////////
  ////IMPORTANT: we're only using bound sqrt root for the adaptive IIT
  /////////////
  
  //////Choose the next neighbor
  vec u = Rcpp::runif(total_neighbors);
  vec probs_choose = log(-log(u)) - logprobs;
  
  //Find the index of the minimum element. 
  //This corresponds to choosing that neighbor
  int neigh_pos = (std::min_element(probs_choose.begin(), probs_choose.end()))-probs_choose.begin();
  // Rcpp::Rcout <<"probs vector: "<< probs << std::endl;
  // Rcpp::Rcout <<"chosen neighbor: "<< neigh_pos << std::endl;
  
  X.row(neigh_pos) = 1-X.row(neigh_pos); //modify the coordinate of the chosen neighbor
  List ret;
  ret["X"]=X;
  ret["Z"]=sum(exp(logprobs))/total_neighbors; // Compute Z factor with uniform proposal distribution
  ret["logbound"]=log_bound;
  ret["max_logbound"]=max_logbound_found;
  return ret;
}

////////// Code to find temperatures for Parallel Tempering //////////

// [[Rcpp::export]]
List temperature_PT_IIT(int p,int interswap, double temp_ini, const std::string bal_function, const double& theta){
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


// [[Rcpp::export]]
List temperature_PT_a_IIT(int p,int interswap, double temp_ini, const std::string bal_function, const double& theta){
  //(int p,int interswap, double temp_ini, const std::string bal_function, const double& theta)
  //sample_inter_swap ----------- interswap
  //// Initialize variables to use in the code
  int T=2; //The total number of temperatures will be 2
  mat X(p,T); // To store the current state of the joint chain, as many rows as neighbors, as many columns as temperatures
  List output; //To store output of the update function
  int new_samples; //To store multiplicity list
  vec log_bound_vector(T,fill::zeros);
  vec max_log_bound_vector(T,fill::zeros);
  double current_log_bound;
  double Z;
  bool update_state=true;
  vec count_iterations(T,fill::zeros);
  int temp_count_iter=0;
  
  //// Parameters for the temperature finding   
  double rho=-3;//1.3;//-2.6; // To define initial second temperature
  double threshold=0.001;//.0003;//.01;////0.001;//Stop when the updates differ less than this much
  double target_swap_rate=0.2345;//Target swap rate
  int count_convergence=0;
  int precision=3;//To round values
  
  bool stay_in_loop=true;
  double swap_prob;
  vec Xtemp_from(p);
  vec Xtemp_to(p);
  //// Initialize temperatures
  vec temp={temp_ini,round_to(temp_ini/(1+exp(rho)),precision)};
  double avg_swap_prob=0;
  
  
  
  
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
    log_bound_vector.zeros();
    max_log_bound_vector.zeros();
      for(int replica=0;replica<T;replica++){//For loop for replicas
        // Rcpp::Rcout << "Sim: " << s+startsim << " Swap: " << i <<"replica: "<<replica<< std::endl;
        int samples_replica=0;
        temp_count_iter=0;
        while(samples_replica<interswap){//Loop to create samples for each replica until we reach the defined threshold
          temp_count_iter+=1;
          current_temp=temp(replica);// Extract temperature of the replica
          current_log_bound=log_bound_vector(replica);// Extract log-bound of the corresponding temperature
  
          // We always update the bounding constant, it can increase and never decrease
          output=a_IIT_update(X.col(replica),Q_matrix,bal_function,current_temp,current_log_bound,true,0,0,max_log_bound_vector(replica),theta);

          update_state=true;
          //// Compute weight
          Z = output(1); //Extract the Z-factor
          new_samples=1+R::rgeom(Z);
          if(new_samples<1){
            Rcpp::Rcout <<"Error: geometric in " << " Swap: " << swap_count <<" temperature:"<<current_temp<< std::endl;
            Rcpp::Rcout <<"new_samples= "<<new_samples<< ", Z=" << Z << " log-bound= " << current_log_bound << std::endl;
          }
          if((samples_replica+new_samples)>interswap){//If we're going to surpass the required number of samples
            new_samples = interswap-samples_replica;//We force to stop at sample_inter_swap
            update_state=false;
            if(swap_count==0){
              count_iterations(replica)=temp_count_iter;
            }else{
              count_iterations(replica)=(count_iterations(replica)*swap_count + temp_count_iter)/(swap_count+1);
            }
            
          }

         
          ///// Updating before the next iteration of the loop
          samples_replica+=new_samples; // Update number of samples obtained from the replica
          if(update_state){
            X.col(replica)=vec(output(0)); //Update current state of the chain
          }
          log_bound_vector(replica)=output(2); //Update log-bound
          max_log_bound_vector(replica)=output(3); //Update maximum log-bound found
        }//End loop to update a single replica
      }//End loop to update all replicas


      swap_count+=1;//Increase the count of swaps


        Xtemp_from=X.col(0);
        Xtemp_to=X.col(1);

        //// Computing swap probability
        swap_prob=(temp(0)-temp(1))*(loglik(Xtemp_to,Q_matrix,theta) - loglik(Xtemp_from,Q_matrix,theta));
        swap_prob=ret_min(exp(swap_prob),1,1);
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
        // Rcpp::Rcout <<"Target swap_rate: "<<target_swap_rate << std::endl; 
        // Rcpp::Rcout <<"Threshold: "<<threshold << std::endl; 
        // Rcpp::Rcout <<"Diff swap_rate: "<<target_swap_rate-avg_swap_prob << std::endl; 
        // Rcpp::Rcout <<"Abs.Diff swap_rate: "<<abs(target_swap_rate-avg_swap_prob) << std::endl; 
        
        if((target_swap_rate-avg_swap_prob)<threshold && (target_swap_rate-avg_swap_prob)>-threshold){
          count_convergence+=1;
          if(count_convergence>=10){
            stay_in_loop=false; 
          }

        }else{
          count_convergence=0;
        }
        
        
        if(swap_count % 500 == 0){
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
  List ret;
  
  ret["temp"]=round_to(temp(1),precision);
  ret["iter"]=count_iterations;
  ret["swap"]=swap_count;
  ret["swap_rate"]=round_to(avg_swap_prob,precision*2);
  ret["seconds"]=duration;
  return ret;
}



