#include <RcppArmadillo.h>
#include <RcppParallel.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
using namespace Rcpp;
using namespace RcppParallel;
using namespace arma;



//This function we don't export to R because it generates errors
double invoke(double x, double (*func)(double)) {
  return func(x);
}
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


double loglik_R(IntegerVector& X,const IntegerMatrix& M, const double& theta){
  double loglik_computed;
  if(M.ncol()==2){
    //Define vectors to represent each mode
    IntegerVector mod1=M(_,0);
    IntegerVector mod2=M(_,1);
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

double loglik_R_par(IntegerVector& X,const IntegerMatrix& M, const double& theta){
  double loglik_computed;
  if(M.ncol()==2){
    //Define vectors to represent each mode
    IntegerVector mod1=M(_,0);
    IntegerVector mod2=M(_,1);
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

// mat initializeRandom(const int& num_rows,const int& num_cols, const double& prob) {
//   
//   // Initialize a matrix with random values between 0 and 1
//   arma::mat A = arma::randu<arma::mat>(num_rows, num_cols);
//   
//   // Threshold the random values to 0 or 1
//   A = arma::conv_to<arma::mat>::from(A > prob);
//   
//   return(A);
// }

[[Rcpp::export]]
IntegerMatrix initializeRandom(const int num_rows,const int num_cols, const double prob) {
  // Initialize a matrix with random values between 0 and 1
  arma::mat A(num_rows, num_cols,fill::randu);
  // Threshold the random values to 0 or 1
  A = arma::conv_to<arma::Mat<double>>::from(A > prob);

  arma::Mat<int> B = arma::conv_to<arma::Mat<int>>::from(A);

  return(Rcpp::wrap(B));
}


/////////////// Parallel workers
struct GibbsSampler : public Worker
{
  //// Data
  RVector<int> X_n;//State matrix
  RVector<int> X_mode;//State matrix
  const double temperature;//Chosen temperature
  const double theta; //Parameter for likelihood
  RVector<int> output;
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
    IntegerVector X_n_in,
    IntegerVector X_mode_in,
    const double temperature_in,
    const double theta_in,
    IntegerVector output_in,
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

// Inside here there are 2 functions
//loglik_internal
//apply_bal_func_internal
//These need to be defined after
struct IIT_visit_neighbors : public Worker
{
  //// Data
  RVector<int> X_n;//State matrix
  const RMatrix<int> Q_matrix; //Matrix with modes
  const int bal_func; //Specified balancing function
  const double temperature;//Chosen temperature
  const double theta; //Parameter for likelihood
  RVector<int> output;
  const std::size_t dim_p; // dimension of the problem (and length of X columns)
  
  //// Functions to use from outside
  double loglik_internal(IntegerVector X,const IntegerMatrix M, const double theta);
  
  double apply_bal_func_internal(double x,const int chosen);
  ////  Functions to transform data
  IntegerVector convert_X(){
    RVector<int> tmp_X = X_n;
    Rcpp::IntegerVector VEC(tmp_X.begin(),tmp_X.end());
    // Rcpp::IntegerVector VEC;
    // for(size_t i=0;i<dim_p;i++){
    //   VEC(i)=tmp_X(i);
    // }
    return VEC;
  }
  IntegerMatrix convert_Q(){
    RMatrix<int> tmp_matrix = Q_matrix;
    Rcpp::IntegerMatrix MAT(tmp_matrix.nrow(),tmp_matrix.ncol());
    std::copy(tmp_matrix.begin(),tmp_matrix.end(),MAT.begin());
    // for(size_t i=0;i<dim_p;i++){
    //   for(size_t j=0;j<2;j++){
    //     MAT(i,j)=tmtmp_matrixp_X(i,j);
    //   }
    // }
    return MAT;
  }
  
  //// Main constructor
  IIT_visit_neighbors(
    IntegerVector X_n_in,
    const IntegerMatrix Q_matrix_in,
    const int bal_func_in,
    const double temperature_in,
    const double theta_in,
    IntegerVector output_in,
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
    IntegerVector X = convert_X();
    IntegerMatrix M = convert_Q();
    // int dim_p_int = X.n_rows;
    double logpi_current=loglik_internal(X,M,theta);
    for(std::size_t n = begin; n < end;n++){ // For for each neighbor
      IntegerVector temp_X=clone(X);//Important to use clone
      temp_X[n]=1-temp_X[n];//Switch the nth coordinate
      double mid_step=loglik_internal(temp_X,M,theta)-logpi_current;
      output[n]=apply_bal_func_internal(mid_step*temperature,bal_func);
      // output[n]=1;
    }// End for loop
    
  }//End operator
};// End of worker to visit neighbors

struct SumExp : public Worker
{   
  
  const RVector<int> X;// source vector
  double Z; // accumulated value
  
  // constructors
  SumExp(const IntegerVector X) : X(X), Z(0) {}
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


///// Full definition of internal functions of workers 
double IIT_visit_neighbors::apply_bal_func_internal(double x,const int chosen){
  return apply_bal_func(x,chosen);
}

// double IIT_visit_neighbors::loglik_internal(const arma::Col<double>& X,const arma::Mat<double>& M, const double& theta){
//   return loglik(X,M,theta);
// }

double IIT_visit_neighbors::loglik_internal(IntegerVector X,const IntegerMatrix M, const double theta){
  return loglik_R(X,M,theta);
}

// Function to initialize the rng according to the seed
std::vector<std::mt19937_64> initialize_rngs(int n, int base_seed) {
  std::vector<std::mt19937_64> rngs(n);
  for(int i = 0; i < n; i++) {
    rngs[i].seed(base_seed + i);  // Unique seed for each RNG
  }
  return rngs;
}


double round_to(double value, int decimals){
  double mult=pow(10.0,decimals);
  value = round(value*mult)/mult;
  return value;
}




// [[Rcpp::export]]
List find_temp_gibbs(int p, int num_iter, int burn_in, double temp_ini,int bal_func, double theta, int base_seed, bool adapting_factors){
  int T=2; //Number of temperatures (for finding the correct temps)
  
  //Initializations
  bool stay_in_loop=true;
  int precision=3;
  double pu=randu();
  IntegerMatrix X = initializeRandom(p,T,pu);//Initialize states
  std::vector<std::mt19937_64> rngs = initialize_rngs(p,base_seed);//Initialize RNGs
  const std::size_t end = static_cast <size_t> (p); 
  
  IntegerVector Xtemp_from(p);
  IntegerVector Xtemp_to(p);

  
  //// Parameters for the temperature finding   
  double rho=-3.5;//1.3;//-2.6; // To define initial second temperature
  double threshold=0.001;//0.001;//.1;//0.001;//Stop when the updates differ less than this much
  double target_swap_rate=0.2345;//Target swap rate
  vec temp_vector={temp_ini,round_to(temp_ini/(1+exp(rho)),precision)};//Initialize temperatures
  double avg_swap_prob=0;
  int count_convergence=0;
  int swap_count=0; //To count how many swaps were performed
  double swap_prob;
  
  //// Define two modes
  IntegerMatrix Q_matrix(p,2);
  // mat Q_matrix(p,2);
  for(int i=0;i<p;i++){
    if(i%2==0){Q_matrix(i,1)=1;}
    if(i%2==1){Q_matrix(i,0)=1;}
  }
  
  
  IntegerMatrix chain_m1=clone(X); // To keep track of how the chain for mode 1 evolves
  IntegerMatrix chain_m2=clone(X);// To keep track of how the chain for mode 2 evolves
  
  //Create mode vectors
  IntegerVector X_mode1(p);
  IntegerVector X_mode2(p);
  for(int i=0;i<p;i++){
    if(i%2==0){X_mode1[i]=1;}
    if(i%2==1){X_mode2[i]=1;}
  }

  
///// Burn-in period  
  
  for(int b=0; b<burn_in; b++){//Burn-in iterations
    for(int replica=0; replica<T;replica++){
      // Rcpp::Rcout <<"Chain 1\n"<<chain_m1<<"\n Chain 2 \n"<<chain_m2 << std::endl;
      IntegerVector output_m1(p);
      GibbsSampler g_sample_m1(chain_m1(_,replica),
                            X_mode1,
                            temp_vector(replica),
                            theta,
                            output_m1,
                            end,
                            rngs);
      parallelFor(0,end,g_sample_m1);//Apply for
      // chain_m1(_,replica)=output_m1; //Update chain
      
      IntegerVector output_m2(p);
      GibbsSampler g_sample_m2(chain_m2(_,replica),
                            X_mode2,
                            temp_vector(replica),
                            theta,
                            output_m2,
                            end,
                            rngs);
      parallelFor(0,end,g_sample_m2);//Apply for
//Update chain

      chain_m1(_,replica)=output_m1;
      chain_m2(_,replica)=output_m2;
      
    }
  }
  
// Interswap iterations  
// while(stay_in_loop){
for(int temp_counter=0;temp_counter<3;temp_counter++){
  for(int i=0; i<num_iter;i++){
    for(int replica=0; replica<T;replica++){
      // Rcpp::Rcout <<"Chain 1\n"<<chain_m1<<"\n Chain 2 \n"<<chain_m2 << std::endl;
      IntegerVector output_m1(p);
      IntegerVector output_m2(p);
      GibbsSampler g_sample_m1(chain_m1(_,replica),
                               X_mode1,
                               temp_vector(replica),
                               theta,
                               output_m1,
                               end,
                               rngs);
      parallelFor(0,end,g_sample_m1);//Apply for
      // chain_m1(_,replica)=output_m1; //Update chain
      
      GibbsSampler g_sample_m2(chain_m2(_,replica),
                               X_mode2,
                               temp_vector(replica),
                               theta,
                               output_m2,
                               end,
                               rngs);
      parallelFor(0,end,g_sample_m2);//Apply for
      // chain_m2(_,replica)=output_m2; //Update chain

      // for(int coord=0;coord<p;coord++){
      //   chain_m1(coord,replica)=output_m1(coord);
      //   chain_m2(coord,replica)=output_m2(coord);
      // }
      chain_m1(_,replica)=output_m1;
      chain_m2(_,replica)=output_m2;
      
      //Update state of the chain based on the mixture with the same weight
      double ppp1=randu();
      if(ppp1<0.5){X(_,replica)=output_m1;}else{X(_,replica)=output_m2;}
    }

  }
  
  // Rcpp::Rcout <<"Chain 1\n"<<chain_m1<<"\n Chain 2 \n"<<chain_m2 << std::endl;
  // Rcpp::Rcout <<"State X:\n"<<X << std::endl;
//Compute replica swap probabilities

swap_count+=1;//Increase the count of swaps
  Xtemp_from=X(_,0);
  Xtemp_to=X(_,1);
  // arma::Col<double> Xtemp_from(X.column(0).begin(),X.column(0).end());
  // arma::Col<double> Xtemp_to(X.column(1).begin(),X.column(1).end());
  double Z_fact_correc=1;//to temporarily store Z_factor correction
  //Declare vector to store info of visiting neighbors
  IntegerVector output_Z_v11(p); 
  IntegerVector output_Z_v22(p); 
  IntegerVector output_Z_v12(p); 
  IntegerVector output_Z_v21(p); 
  if(adapting_factors){//For PT-IIT
  
  //// Computation of Z factors to correct bias
  
  
  //// Declare constructor to visit all neighbors
  IIT_visit_neighbors Z_v11(Xtemp_from,
                            Q_matrix,
                            bal_func,
                            temp_vector(0),
                            theta,
                            output_Z_v11,
                            end);
  IIT_visit_neighbors Z_v22(Xtemp_to,
                            Q_matrix,
                            bal_func,
                            temp_vector(1),
                            theta,
                            output_Z_v22,
                            end);
  IIT_visit_neighbors Z_v12(Xtemp_from,
                            Q_matrix,
                            bal_func,
                            temp_vector(1),
                            theta,
                            output_Z_v12,
                            end);
  IIT_visit_neighbors Z_v21(Xtemp_to,
                            Q_matrix,
                            bal_func,
                            temp_vector(0),
                            theta,
                            output_Z_v21,
                            end);
  // Apply ParallelFor
  parallelFor(0,end,Z_v11);
  parallelFor(0,end,Z_v22);
  parallelFor(0,end,Z_v12);
  parallelFor(0,end,Z_v21);

  //// Declare constructor to add log-probabilities
  SumExp sum_Z11(output_Z_v11);
  SumExp sum_Z22(output_Z_v22);
  SumExp sum_Z12(output_Z_v12);
  SumExp sum_Z21(output_Z_v21);
  // //// Get the sum of probabilities
  parallelReduce(0,end,sum_Z11);
  parallelReduce(0,end,sum_Z22);
  parallelReduce(0,end,sum_Z12);
  parallelReduce(0,end,sum_Z21);

  // Z_fact_correc=(sum_Z12.Z*sum_Z21.Z)/(sum_Z11.Z*sum_Z22.Z);
  
}else{//For A-IIT
  Z_fact_correc=1;
}
  
//// Computing swap probability
swap_prob=(temp_vector(0)-temp_vector(1))*(loglik_R(Xtemp_to,Q_matrix,theta) - loglik_R(Xtemp_from,Q_matrix,theta)); 
swap_prob=ret_min(Z_fact_correc*exp(swap_prob),1,1);

// Rcpp::Rcout <<"Prob TO: "<<loglik_R(Xtemp_to,Q_matrix,theta)<<"\n Prob FROM "<<loglik_R(Xtemp_from,Q_matrix,theta) << std::endl;


if(swap_count==1){
  avg_swap_prob=swap_prob;
}else{
  avg_swap_prob=(avg_swap_prob*(swap_count-1)+swap_prob)/swap_count;
}

//// Update temperature
rho=rho + (swap_prob-target_swap_rate)/swap_count;


if(rho<1e-5){
  temp_vector(1)=temp_ini/(2+expm1(rho));
}else{
  temp_vector(1)=temp_ini/(1+exp(rho));
}

//Check current threshold
if((target_swap_rate-avg_swap_prob)<threshold && (target_swap_rate-avg_swap_prob)>-threshold){
  count_convergence+=1;
  if(count_convergence>=3){
    stay_in_loop=false;
  }
  
}else{
  count_convergence=0;
}


if(swap_count % 100 == 0){
  Rcpp::Rcout <<"Swap: "<<swap_count<<" avg. swap prob: "<<avg_swap_prob <<" new temperature: "<< temp_vector(1) << std::endl; 
}

if(swap_count == 700000){// Force finishing of algorithm
  stay_in_loop=false;
} 

}
Rcpp::Rcout <<"FINAL RESULTS:\nSwap: "<<swap_count<<" avg. swap prob: "<<avg_swap_prob <<" new temperature: "<< temp_vector(1) << std::endl; 

List ret;
ret["temp"]=round_to(temp_vector(1),precision);
ret["swap"]=swap_count;
ret["swap_rate"]=round_to(avg_swap_prob,precision*2);

return ret;  

}


// [[Rcpp::export]]
std::size_t test_flip_coord(std::size_t coord,bool approaching,double theta){
  std::size_t new_coord=coord;//To store new coordinate
  double u=randu();
  if(approaching){//In case the proposed coordinate is approaching the mode
    if(u<(1/(1+exp(-theta)))){new_coord=(1-coord);}
    // Rcpp::Rcout <<"Accept approach"<< std::endl;
  }else{//In case the proposed coordinate is moving away from the mode
    if(u<(1/(1+exp(theta)))){new_coord=(1-coord);}
    // Rcpp::Rcout <<"Accept move far"<< std::endl;
  }
  return new_coord;
}

