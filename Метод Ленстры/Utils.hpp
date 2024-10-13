#include <gmpxx.h>
#include <random>
#include <tuple>

std::tuple<mpz_class, mpz_class, mpz_class> extended_euclid(mpz_class& a, mpz_class& b);

void eratosthenes(std::vector<int>& primes, int n);

bool is_prime(std::mt19937& mt, mpz_class& number, int steps);
