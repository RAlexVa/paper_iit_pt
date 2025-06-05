#include <Rcpp.h>
// [[Rcpp::depends(dqrng, BH)]]
#include <dqrng_generator.h>
#include <dqrng_sample.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

struct RandomFill : public RcppParallel::Worker {
  RcppParallel::RMatrix<int> output;
  uint64_t seed;
  
  RandomFill(Rcpp::IntegerMatrix output, const uint64_t seed) : output(output), seed(seed) {};
  
  void operator()(std::size_t begin, std::size_t end) {
    auto rng = dqrng::generator<pcg64>(seed, end);
    for (std::size_t col = begin; col < end; ++col) {
      auto sampled = dqrng::sample::sample<std::vector<int>, uint32_t>(*rng, 100000, output.nrow(), true);
      RcppParallel::RMatrix<int>::Column column = output.column(col);
      std::copy(sampled.begin(), sampled.end(), column.begin());
    }
  }
};

// [[Rcpp::export]]
Rcpp::IntegerMatrix parallel_random_matrix(const int n, const int m, const int ncores, int seed_def) {
  Rcpp::IntegerMatrix res(n, m);
  RandomFill randomFill(res, seed_def);
  RcppParallel::parallelFor(0, m, randomFill, m/ncores + 1);
  return res;
}