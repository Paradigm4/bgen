
#include <array>
#include <algorithm>
#include <stdexcept>

#include "utils.h"


namespace bgen {

// Returns value of Binomial Coefficient C(n, k)
std::uint32_t n_choose_k(int n, int k) {
  std::uint32_t res = 1;

  // Since C(n, k) = C(n, n-k)
  if ( k > n - k ) {
    k = n - k;
  }

  // Calculate value of
  // [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
  for (std::uint32_t i = 0; i < (std::uint32_t) k; ++i) {
    res *= (n - i);
    res /= (i + 1);
  }

  return res;
}

/// check if the minor allele is certain (to 99.9999999999999& confidence)
///
///  Take the frequency, and number of individuals checked so far, and see if the
///  99.99..(fifteen nines) confidence interval overlaps 0.5. If not, then we can
///  be sure we've identified the minor allele, even without checking the full
///  population.
///
///  @param freq estimated minor allele frequency
///  @param n_checked number of individsuals checked so far
///  @param z standard normal deviate (eg 1.96 for 95% CI, here we use 10.0 for
///    stronger confidence, and the fact the normal approximation for confidence
///    intervals isn't perfect)
///  @return True/False for whether to halt the permuations
bool minor_certain(double freq, int n_checked, double z) {
    double delta = (z * std::sqrt((freq * (1 - freq)) / n_checked));
    
    // check if the confidence interval overlaps 0.5
    return !((freq - delta < 0.5) & (freq + delta > 0.5));
}

// sum array of 8-bit uints with vectorised operations
//
// Summing the ploidy array via numpy is more expensive than it should be, this
// is just a quick vectorized sum to speed that up.
//
// We need to load the data into vector registers, but we need to convert the
// 8-bit vectors to 32-bit to avoid overflow. 32-bit uints should be sufficient,
// since this sums ploidy states, which can be at most 255 per person, so this
// allows at least 269 million individuals (((2 ** 32) * 16) / 255).
//
/// @param x array of floats
/// @param size size of array
/// @returns sum of array
std::uint64_t fast_ploidy_sum(std::uint8_t * x, std::uint32_t & size) {
  std::uint32_t i = 0;
  std::uint64_t total = 0;
  
  if (size > 269000000) {
    // raise error if this gets too many samples. the fix would be to refactor
    // this function to use packed 64-bit ints e.g. _mm512_add_epi64
    throw std::invalid_argument("too many samples for valid summing");
  }

  // include the remainder not used during vectorised operations
  for ( ; i < size; i++) {
    total += x[i];
  }
  return total;
}

// get min and max of the ploidy values in one pass
Range fast_range(std::uint8_t * x, std::uint32_t & size) {
  std::uint8_t min_val = 255;
  std::uint8_t max_val = 0;
  size_t i = 0;



  // include the remainder not used during vectorised operations
  for ( ; i < size; i++) {
    min_val = std::min(min_val, x[i]);
    max_val = std::max(max_val, x[i]);
  }
  return {min_val, max_val};
}

} // namespace bgen
