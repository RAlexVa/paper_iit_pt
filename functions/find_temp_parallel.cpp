//#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;
using namespace arma;

// [[Rcpp::plugins("cpp17")]] 


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

arma::Mat<double> initializeRandom(const int num_rows,const int num_cols, const double prob) {

  // Initialize a matrix with random values between 0 and 1
  arma::mat A(num_rows, num_cols,fill::randu);

  // Threshold the random values to 0 or 1
  A = arma::conv_to<arma::Mat<double>>::from(A > prob);

  return(A);
}


///// Here just parallelize the update for each replica /////

//// Definition of worker
struct Parallel_replica_update : public Worker
{
  //// Data
  RMatrix<double> X_n_matrix;//State matrix
  // arma::Mat<double> X_n_matrix;//State matrix
  const RMatrix<double> Q_matrix; //Matrix with modes
  // const arma::Mat<double> Q_matrix; //Matrix with modes
  const int bal_func; //Specified balancing function
  const RVector<double> temperatures;//Vector of temperatures
  // const arma::Col<double> temperatures;//Vector of temperatures
  const double theta; //Parameter for 
  RMatrix<double> output;
  // arma::Mat<double> output;
  const std::size_t dim_p; // dimension of the problem (and length of X columns)
  const std::size_t num_temp; // dimension of the problem (and length of X columns)
  
  
  //// Functions to use from outside
  arma::Col<double> IIT_update_p(arma::Col<double> X,
                                 const arma::Mat<double>& M,
                                 const int chosen_bf,
                                 const double& temperature,
                                 const double& theta);
  
  ////  Functions to transform data
  arma::Mat<double> convert_X(){
    RMatrix<double> tmp_matrix = X_n_matrix;
    arma::Mat<double> MAT(tmp_matrix.begin(), dim_p, num_temp, false,true);
    return MAT;
  }
  
  arma::Mat<double> convert_Q(){
    RMatrix<double> tmp_matrix = Q_matrix;
    arma::Mat<double> MAT(tmp_matrix.begin(), dim_p, 2, false,true);
    return MAT;
  }
  
  arma::Col<double> convert_temps(){
    RVector<double> tmp_vec = temperatures;
    arma::Col<double> vec(tmp_vec.begin(),num_temp,false,true);
    return vec;
  }
  
  
  
  //// Main constructor
  Parallel_replica_update(
    NumericMatrix X_n_matrix_in,
    // arma::Mat<double> X_n_matrix_in,
    const NumericMatrix Q_matrix_in,
    // const  arma::Mat<double> Q_matrix_in,
    const int bal_func_in,
    const NumericVector temperature_in,
    // const arma::Col<double> temperature_in,
    const double theta_in,
    NumericMatrix output_in,
    // arma::Mat<double> output_in,
    const std::size_t dim_p_in,
    const std::size_t num_temp_in)://This binds the class members (defined above) to the constructor arguments
    X_n_matrix(X_n_matrix_in),
    Q_matrix(Q_matrix_in),
    bal_func(bal_func_in),
    temperatures(temperature_in),
    theta(theta_in),
    output(output_in),
    dim_p(dim_p_in),
    num_temp(num_temp_in)
    //// dim_p(X_n_matrix.n_rows), //I'm not sure how to get the dim
    //// num_temp(x_n_matrix.n_cols) //I'm not sure how to get the dim
  {}
  
  
  
  //// Operator
  void operator()(std::size_t begin, std::size_t end){
    //First transform all the inputs
    arma::Mat<double> X = convert_X();
    arma::Mat<double> Q_matrix = convert_Q();
    arma::Col<double> temperatures = convert_temps();
    int dim_p_int = static_cast<int>(dim_p);
    for(std::size_t r=begin; r<end;r++){ // For for each replica
      double current_temp = temperatures(r);
      
      arma::Col<double> temporal_X=IIT_update_p(X.col(r),
                                                Q_matrix,
                                                bal_func,
                                                current_temp,
                                                theta);
      // arma::Col<double> temporal_X(dim_p);
      // temporal_X.zeros();
      for(int i=0;i<dim_p_int;i++){
        output(i,r)=temporal_X[i];
      }
      
    }
    
  }//End operator
};// End of worker

// Full definition of update function used inside the worker
arma::Col<double> Parallel_replica_update::IIT_update_p(arma::Col<double> X,
                                                        const arma::Mat<double>& M,
                                                        const int chosen_bf,
                                                        const double& temperature,
                                                        const double& theta){
  int total_neighbors = X.n_rows; // total number of neighbors is p spacial
  vec probs(total_neighbors, fill::zeros); //vector of probabilities
  double logpi_current = loglik(X,M,theta);  // Compute likelihood of current state
  
  ////// Compute weight for all neighbors
  double temporal=0;
  vec newX;
  for(int j=0; j<total_neighbors;j++){
    newX = X;
    newX.row(j) = 1-X.row(j);
    temporal=loglik(newX,M,theta)-logpi_current;
    probs(j)=apply_bal_func(temporal*temperature, chosen_bf);//Apply balancing function to log probability times temperature
  }
  //////Choose the next neighbor
  vec u(total_neighbors,fill::randu);
  // vec u;
  // u.fill(0.5);
  vec probs_choose = log(-log(u)) - probs;
  int neigh_pos = (std::min_element(probs_choose.begin(), probs_choose.end()))-probs_choose.begin();
  X.row(neigh_pos) = 1-X.row(neigh_pos); //modify the coordinate of the chosen neighbor
  return X;
}




//Function to run simulations
// [[Rcpp::export]]
NumericMatrix PT_IIT_parallel_sim(int p, int num_iter, arma::Col<double> original_temperatures, const std::string bal_func, double theta){
  int T=original_temperatures.n_rows;
  
  const std::size_t dim_p = static_cast<size_t>(p);
  const std::size_t num_temps = static_cast<size_t>(T);
  
  const std::size_t begin = static_cast <size_t> (0);
  const std::size_t end = static_cast <size_t> (T); 
  
  // Change String of balancing function into int
  int chosen_bal_func;
  if(bal_func=="sq"){
    chosen_bal_func=2;
  }else if(bal_func=="min"){
    chosen_bal_func=1;
  }else{
    Rcpp::Rcout << "Incorrect definition of Balancing function.\n Defaulting to sq " << std::endl;
    chosen_bal_func=2;
  }
  
  //// Define two modes
  NumericMatrix Q_matrix(p,2);
  for(int i=0;i<p;i++){
    if(i%2==0){Q_matrix(i,1)=1;}
    if(i%2==1){Q_matrix(i,0)=1;}
  }
  
  double ppp=randu();
  arma::Mat<double> original_X(p,T);
  original_X=initializeRandom(p,T,ppp);//Randomly initialize the state of each replica.
  
  
  //After initialization of vectors and matrices we transform them to Nuneric variables
  NumericMatrix X=Rcpp::wrap(original_X);
  NumericVector temperatures=Rcpp::wrap(original_temperatures);
  // Rcpp::Rcout <<"Before parallel X:\n "<<X<< std::endl;
  // Now we start the for loop
  for(int i=0;i<num_iter;i++){
    NumericMatrix output(p,T);
    
    Parallel_replica_update par_rep_update(
        X,Q_matrix,
        chosen_bal_func,
        temperatures,
        theta,
        output,
        dim_p,
        num_temps);// Initialize the worker
    
    
    
    parallelFor(begin,end,par_rep_update);//Perform the for
    X=output;  
  }
  
  return X;
}




/////////////////////////////////////////////////////////////////////////////////
///// Here include the for in the parallelization /////
//// Definition of worker
struct Parallel_replica_update_wfor : public Worker
{
  //// Data
  RMatrix<double> X_n_matrix;//State matrix
  // arma::Mat<double> X_n_matrix;//State matrix
  const RMatrix<double> Q_matrix; //Matrix with modes
  // const arma::Mat<double> Q_matrix; //Matrix with modes
  const int bal_func; //Specified balancing function
  const RVector<double> temperatures;//Vector of temperatures
  // const arma::Col<double> temperatures;//Vector of temperatures
  const double theta; //Parameter for 
  RMatrix<double> output;
  // arma::Mat<double> output;
  const std::size_t dim_p; // dimension of the problem (and length of X columns)
  const std::size_t num_temp; // dimension of the problem (and length of X columns)
  const int num_iter;//Number of iterations in the for loop
  
  
  //// Functions to use from outside
  arma::Col<double> IIT_update_p(arma::Col<double> X,
                                 const arma::Mat<double>& M,
                                 const int chosen_bf,
                                 const double& temperature,
                                 const double& theta);
  
  ////  Functions to transform data
  arma::Mat<double> convert_X(){
    RMatrix<double> tmp_matrix = X_n_matrix;
    arma::Mat<double> MAT(tmp_matrix.begin(), dim_p, num_temp, false,true);
    return MAT;
  }
  
  arma::Mat<double> convert_Q(){
    RMatrix<double> tmp_matrix = Q_matrix;
    arma::Mat<double> MAT(tmp_matrix.begin(), dim_p, 2, false,true);
    return MAT;
  }
  
  arma::Col<double> convert_temps(){
    RVector<double> tmp_vec = temperatures;
    arma::Col<double> vec(tmp_vec.begin(),num_temp,false,true);
    return vec;
  }
  
  
  
  //// Main constructor
  Parallel_replica_update_wfor(
    NumericMatrix X_n_matrix_in,
    // arma::Mat<double> X_n_matrix_in,
    const NumericMatrix Q_matrix_in,
    // const  arma::Mat<double> Q_matrix_in,
    const int bal_func_in,
    const NumericVector temperature_in,
    // const arma::Col<double> temperature_in,
    const double theta_in,
    NumericMatrix output_in,
    // arma::Mat<double> output_in,
    const std::size_t dim_p_in,
    const std::size_t num_temp_in,
    const int num_iter_in)://This binds the class members (defined above) to the constructor arguments
    X_n_matrix(X_n_matrix_in),
    Q_matrix(Q_matrix_in),
    bal_func(bal_func_in),
    temperatures(temperature_in),
    theta(theta_in),
    output(output_in),
    dim_p(dim_p_in),
    num_temp(num_temp_in),
    num_iter(num_iter_in)
    //// dim_p(X_n_matrix.n_rows), //I'm not sure how to get the dim
    //// num_temp(x_n_matrix.n_cols) //I'm not sure how to get the dim
  {}
  
  
  
  //// Operator
  void operator()(std::size_t begin, std::size_t end){
    //First transform all the inputs
    arma::Mat<double> X = convert_X();
    arma::Mat<double> Q_matrix = convert_Q();
    arma::Col<double> temperatures = convert_temps();
    int dim_p_int = static_cast<int>(dim_p);
    
    for(std::size_t r=begin; r<end;r++){ // For for each replica
      double current_temp = temperatures(r);
      arma::Col<double> selected_X=X.col(r); // Get the status of the replica
      arma::Col<double> temporal_X(dim_p);
      for(int k=0;k<num_iter;k++){//Repeat for as many iterations
        temporal_X=IIT_update_p(selected_X,
                                Q_matrix,
                                bal_func,
                                current_temp,
                                theta);
        selected_X=temporal_X;//Assign new state for the loop
      }
      // After num_iter updates, export the resulting file
      for(int i=0;i<dim_p_int;i++){
        output(i,r)=temporal_X[i];
      }
    }
    
  }//End operator
};// End of worker

// Full definition of update function used inside the worker
arma::Col<double> Parallel_replica_update_wfor::IIT_update_p(arma::Col<double> X,
                                                             const arma::Mat<double>& M,
                                                             const int chosen_bf,
                                                             const double& temperature,
                                                             const double& theta){
  int total_neighbors = X.n_rows; // total number of neighbors is p spacial
  vec probs(total_neighbors, fill::zeros); //vector of probabilities
  double logpi_current = loglik(X,M,theta);  // Compute likelihood of current state
  
  ////// Compute weight for all neighbors
  double temporal=0;
  vec newX;
  for(int j=0; j<total_neighbors;j++){
    newX = X;
    newX.row(j) = 1-X.row(j);
    temporal=loglik(newX,M,theta)-logpi_current;
    probs(j)=apply_bal_func(temporal*temperature, chosen_bf);//Apply balancing function to log probability times temperature
  }
  //////Choose the next neighbor
  vec u(total_neighbors,fill::randu);
  // vec u;
  // u.fill(0.5);
  vec probs_choose = log(-log(u)) - probs;
  int neigh_pos = (std::min_element(probs_choose.begin(), probs_choose.end()))-probs_choose.begin();
  X.row(neigh_pos) = 1-X.row(neigh_pos); //modify the coordinate of the chosen neighbor
  return X;
}




//Function to run simulations
// [[Rcpp::export]]
NumericMatrix PT_IIT_parallel_sim_wfor(int p, int num_iter, arma::Col<double> original_temperatures, const std::string bal_func, double theta){
  int T=original_temperatures.n_rows;
  
  const std::size_t dim_p = static_cast<size_t>(p);
  const std::size_t num_temps = static_cast<size_t>(T);
  
  const std::size_t begin = static_cast <size_t> (0);
  const std::size_t end = static_cast <size_t> (T); 
  
  // Change String of balancing function into int
  int chosen_bal_func;
  if(bal_func=="sq"){
    chosen_bal_func=2;
  }else if(bal_func=="min"){
    chosen_bal_func=1;
  }else{
    Rcpp::Rcout << "Incorrect definition of Balancing function.\n Defaulting to sq " << std::endl;
    chosen_bal_func=2;
  }
  
  //// Define two modes
  NumericMatrix Q_matrix(p,2);
  for(int i=0;i<p;i++){
    if(i%2==0){Q_matrix(i,1)=1;}
    if(i%2==1){Q_matrix(i,0)=1;}
  }
  
  double ppp=randu();
  arma::Mat<double> original_X(p,T);
  original_X=initializeRandom(p,T,ppp);//Randomly initialize the state of each replica.
  // original_X.zeros();
  
  //After initialization of vectors and matrices we transform them to Nuneric variables
  NumericMatrix X=Rcpp::wrap(original_X);
  NumericVector temperatures=Rcpp::wrap(original_temperatures);
  // Rcpp::Rcout <<"Before parallel X:\n "<<X<< std::endl;
  // Now we apply the parallelized loop
  
  NumericMatrix output(p,T);//Declare output
  
  Parallel_replica_update_wfor par_rep_update_wfor(
      X,Q_matrix,
      chosen_bal_func,
      temperatures,
      theta,
      output,
      dim_p,
      num_temps,
      num_iter);// Initialize the worker
  
  parallelFor(begin,end,par_rep_update_wfor);//Perform the for
  
  
  return output;
}



/////////////////////////////////////////////////////////////////////////////////
///// Here include the for in the parallelization AND try to make it replicable /////
//// Definition of worker
struct Parallel_replica_update_rep : public Worker
{
  //// Data
  RMatrix<double> X_n_matrix;//State matrix
  const RMatrix<double> Q_matrix; //Matrix with modes
  const int bal_func; //Specified balancing function
  const RVector<double> temperatures;//Vector of temperatures
  const double theta; //Parameter for 
  RMatrix<double> output;
  const std::size_t dim_p; // dimension of the problem (and length of X columns)
  const std::size_t num_temp; // dimension of the problem (and length of X columns)
  const int num_iter;//Number of iterations in the for loop
  std::vector<std::mt19937_64>& rngs;//Vector of RNGs
  //// Functions to use from outside
  arma::Col<double> IIT_update_p(arma::Col<double> X,
                                 const arma::Mat<double>& M,
                                 const int chosen_bf,
                                 const double& temperature,
                                 const double& theta,
                                 std::mt19937_64& rng);
  
  ////  Functions to transform data
  arma::Mat<double> convert_X(){
    RMatrix<double> tmp_matrix = X_n_matrix;
    arma::Mat<double> MAT(tmp_matrix.begin(), dim_p, num_temp, false,true);
    return MAT;
  }
  arma::Mat<double> convert_Q(){
    RMatrix<double> tmp_matrix = Q_matrix;
    arma::Mat<double> MAT(tmp_matrix.begin(), dim_p, 2, false,true);
    return MAT;
  }
  arma::Col<double> convert_temps(){
    RVector<double> tmp_vec = temperatures;
    arma::Col<double> vec(tmp_vec.begin(),num_temp,false,true);
    return vec;
  }
  
  //// Main constructor
  Parallel_replica_update_rep(
    NumericMatrix X_n_matrix_in,
    const NumericMatrix Q_matrix_in,
    const int bal_func_in,
    const NumericVector temperature_in,
    const double theta_in,
    NumericMatrix output_in,
    const std::size_t dim_p_in,
    const std::size_t num_temp_in,
    const int num_iter_in,
    std::vector<std::mt19937_64>& rngs_in)://This binds the class members (defined above) to the constructor arguments
    X_n_matrix(X_n_matrix_in),
    Q_matrix(Q_matrix_in),
    bal_func(bal_func_in),
    temperatures(temperature_in),
    theta(theta_in),
    output(output_in),
    dim_p(dim_p_in),
    num_temp(num_temp_in),
    num_iter(num_iter_in),
    rngs(rngs_in)
  {}
  
  //// Operator
  void operator()(std::size_t begin, std::size_t end){
    //First transform all the inputs
    arma::Mat<double> X = convert_X();
    arma::Mat<double> Q_matrix = convert_Q();
    arma::Col<double> temperatures = convert_temps();
    int dim_p_int = static_cast<int>(dim_p);
    
    for(std::size_t r=begin; r<end;r++){ // For for each replica
      double current_temp = temperatures(r);
      arma::Col<double> selected_X=X.col(r); // Get the status of the replica
      arma::Col<double> temporal_X(dim_p);
      std::mt19937_64 selected_rng=rngs[r];
      for(int k=0;k<num_iter;k++){//Repeat for as many iterations
        temporal_X=IIT_update_p(selected_X,
                                Q_matrix,
                                bal_func,
                                current_temp,
                                theta,
                                selected_rng);
        selected_X=temporal_X;//Assign new state for the loop
      }
      // After num_iter updates, export the resulting file
      for(int i=0;i<dim_p_int;i++){
        output(i,r)=temporal_X[i];
      }//End loop for the export 
    }// End loop for each replica
    
  }//End operator
};// End of worker

// Full definition of update function used inside the worker
arma::Col<double> Parallel_replica_update_rep::IIT_update_p(arma::Col<double> X,
                                                            const arma::Mat<double>& M,
                                                            const int chosen_bf,
                                                            const double& temperature,
                                                            const double& theta,
                                                            std::mt19937_64& rng){
  int total_neighbors = X.n_rows; // total number of neighbors is p spacial
  vec probs(total_neighbors, fill::zeros); //vector of probabilities
  double logpi_current = loglik(X,M,theta);  // Compute likelihood of current state
  
  ////// Compute weight for all neighbors
  double temporal=0;
  vec newX;
  for(int j=0; j<total_neighbors;j++){
    newX = X;
    newX.row(j) = 1-X.row(j);
    temporal=loglik(newX,M,theta)-logpi_current;
    probs(j)=apply_bal_func(temporal*temperature, chosen_bf);//Apply balancing function to log probability times temperature
  }
  //////Choose the next neighbor
  // vec u(total_neighbors,fill::randu);
  vec u(total_neighbors);
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  for(int j=0; j<total_neighbors; j++) {
    u(j) = dist(rng);
  }
  vec probs_choose = log(-log(u)) - probs;
  int neigh_pos = (std::min_element(probs_choose.begin(), probs_choose.end()))-probs_choose.begin();
  X.row(neigh_pos) = 1-X.row(neigh_pos); //modify the coordinate of the chosen neighbor
  return X;
}



// Function to initialize the rng according to the seed
std::vector<std::mt19937_64> initialize_rngs(int n, int base_seed) {
  std::vector<std::mt19937_64> rngs(n);
  for(int i = 0; i < n; i++) {
    rngs[i].seed(base_seed + i);  // Unique seed for each RNG
  }
  return rngs;
}
//Function to run simulations
// [[Rcpp::export]]
NumericMatrix PT_IIT_parallel_sim_rep(int p, int num_iter, arma::Col<double> original_temperatures, const std::string bal_func, double theta, int base_seed){
  int T=original_temperatures.n_rows;
  const std::size_t dim_p = static_cast<size_t>(p);
  const std::size_t num_temps = static_cast<size_t>(T);
  const std::size_t begin = static_cast <size_t> (0);
  const std::size_t end = static_cast <size_t> (T); 
  // Change String of balancing function into int
  int chosen_bal_func;
  if(bal_func=="sq"){
    chosen_bal_func=2;
  }else if(bal_func=="min"){
    chosen_bal_func=1;
  }else{
    Rcpp::Rcout << "Incorrect definition of Balancing function.\n Defaulting to sq " << std::endl;
    chosen_bal_func=2;
  }
  //// Define two modes
  NumericMatrix Q_matrix(p,2);
  for(int i=0;i<p;i++){
    if(i%2==0){Q_matrix(i,1)=1;}
    if(i%2==1){Q_matrix(i,0)=1;}
  }
  double ppp=randu();
  arma::Mat<double> original_X(p,T);
  original_X=initializeRandom(p,T,ppp);//Randomly initialize the state of each replica.
  std::vector<std::mt19937_64> defined_rngs=initialize_rngs(T,base_seed); // Initialize RNGs
  //After initialization of vectors and matrices we transform them to Nuneric variables
  NumericMatrix X=Rcpp::wrap(original_X);
  NumericVector temperatures=Rcpp::wrap(original_temperatures);
  // Now we apply the parallelized loop
  
  NumericMatrix output(p,T);//Declare output
  Parallel_replica_update_rep par_rep_update_rep(//Declare Worker
      X,Q_matrix,
      chosen_bal_func,
      temperatures,
      theta,
      output,
      dim_p,
      num_temps,
      num_iter,
      defined_rngs);// Initialize the worker
  
  parallelFor(begin,end,par_rep_update_rep);//Perform the for loop
  
  return output;
}












/////////////////////////////////////////////////////////////////////////////////
///// Here Parallelize IIT_update /////

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
    double logpi_current=loglik(X,M,theta);
      
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


// Full definition of internal functions
double IIT_visit_neighbors::apply_bal_func_internal(double x,const int chosen){
  return apply_bal_func(x,chosen);
}

double IIT_visit_neighbors::loglik_internal(const arma::Col<double>& X,const arma::Mat<double>& M, const double& theta){
  return loglik(X,M,theta);
}


// [[Rcpp::export]]
List test_1(int p, double temperature, int bal_func, double theta, int base_seed){
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


