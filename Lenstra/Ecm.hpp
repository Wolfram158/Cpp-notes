#include <Utils.cpp>
#include <atomic>
#include <mutex>
#include <gmpxx.h>

class Lenstra_ECM {
public:
    Lenstra_ECM();
    ~Lenstra_ECM();
    void factor(mpz_t& n, int B, mpz_t& C);
    mpz_t* get_result();
    void set_uncompletness();
private:
    std::vector<int> primes;
    mpz_t n;
    mpz_t result;
    int B;
    std::mutex mtx;
    std::atomic<bool> complete;

    bool is_not_satisfied(gmp_randstate_t& state, mpz_t& m);
    bool add_points(
        mpz_t &n,
        mpz_t &a,
        mpz_t &x1,
        mpz_t &y1,
        bool &o1,
        mpz_t &x2,
        mpz_t &y2,
        bool &o2,
        mpz_t& inverse,
        mpz_t& gcd,
        mpz_t& inv,
        mpz_t& input,
        mpz_t& alpha,
        mpz_t& beta,
        mpz_t& x11,
        mpz_t& y11,
        mpz_t& g,
        mpz_t& r,
        mpz_t& s,
        mpz_t& t,
        mpz_t& q,
        mpz_t& u,
        mpz_t& v,
        mpz_t& w
    );
    bool multiply_point(
        mpz_t& n, 
        mpz_t& a, 
        long long num, 
        mpz_t& x,
        mpz_t& y,
        mpz_t& pow_x,
        mpz_t& pow_y,
        mpz_t& inverse,
        mpz_t& gcd,
        mpz_t& inv,
        mpz_t& input,
        mpz_t& alpha,
        mpz_t& beta,
        mpz_t& x11,
        mpz_t& y11,
        mpz_t& g,
        mpz_t& r,
        mpz_t& s,
        mpz_t& t,
        mpz_t& q,
        mpz_t& u,
        mpz_t& v,
        mpz_t& w,
        int steps
    );
    bool try_ecm(
        mpz_t& n, 
        mpz_t& a, 
        int B, 
        mpz_t& C, 
        mpz_t& x,
        mpz_t& y,
        mpz_t& pow_x,
        mpz_t& pow_y,
        mpz_t& inverse,
        mpz_t& gcd,
        mpz_t& inv,
        mpz_t& input,
        mpz_t& alpha,
        mpz_t& beta,
        mpz_t& x11,
        mpz_t& y11,
        mpz_t& g,
        mpz_t& r,
        mpz_t& s,
        mpz_t& t,
        mpz_t& q,
        mpz_t& u,
        mpz_t& v,
        mpz_t& w
    );
};