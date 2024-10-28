#include <gmp.h>
#include <random>

void extended_euclid(mpz_t& a, mpz_t& b, mpz_t& gcd, mpz_t& inv, mpz_t& x, mpz_t& y, mpz_t& g, mpz_t& r, mpz_t& s, mpz_t& t, mpz_t& q, mpz_t& u, mpz_t& v, mpz_t& w) {
    // mpz_t x;
    // mpz_t y;
    // mpz_t g;
    // mpz_t r;
    // mpz_t s;
    // mpz_t t;

    // mpz_init_set_ui(x, 1);
    // mpz_init_set_ui(y, 0);
    // mpz_init_set(g, a);
    // mpz_init_set_ui(r, 0);
    // mpz_init_set_ui(s, 1);
    // mpz_init_set(t, b);

    mpz_set_si(x, 1);
    mpz_set_si(y, 0);
    mpz_set(g, a);
    mpz_set_si(r, 0);
    mpz_set_si(s, 1);
    mpz_set(t, b);

    // mpz_t q;
    // mpz_t u;
    // mpz_t v;
    // mpz_t w;

    // mpz_init(q);
    // mpz_init(u);
    // mpz_init(v);
    // mpz_init(w);

    while (mpz_cmp_si(t, 0) > 0) {
        mpz_fdiv_q(q, g, t);
        mpz_set(u, x);
        mpz_submul(u, q, r);

        mpz_set(v, y);
        mpz_submul(v, q, s);

        mpz_set(w, g);
        mpz_submul(w, q, t);

        mpz_set(x, r);
        mpz_set(y, s);
        mpz_set(g, t);
        mpz_set(r, u);
        mpz_set(s, v);
        mpz_set(t, w);
    }

    mpz_set(gcd, g);
    mpz_set(inv, x);

    // mpz_clear(x);
    // mpz_clear(y);
    // mpz_clear(g);
    // mpz_clear(r);
    // mpz_clear(s);
    // mpz_clear(t);
    // mpz_clear(q);
    // mpz_clear(u);
    // mpz_clear(v);
    // mpz_clear(w);
}

void eratosthenes(std::vector<int>& primes, int n) {
    std::vector<int> m(n + 1, 0);
    for (int i = 2; i < n; i++) {
        if (m[i] == 0) {
            m[i] = i;
            primes.push_back(i);
        }
        for (std::vector<int>::size_type j = 0; j < primes.size() && primes[j] <= m[i] && 
            i * primes[j] <= n; j++) {
            m[i * primes[j]] = primes[j];
        }
    }
}

void fast_power_mod(mpz_t& a, mpz_t& u, mpz_t& mod, mpz_t& result) {
    if (mpz_cmp_ui(u, 0) == 0) {
        mpz_set_ui(result, 1);
        return;
    }
    mpz_t half;
    mpz_init(half);
    mpz_fdiv_q_ui(half, u, 2);
    mpz_t sqrt;
    mpz_init(sqrt);
    fast_power_mod(a, half, mod, sqrt);
    mpz_t umod;
    mpz_init(umod);
    mpz_mod_ui(umod, u, 2);
    mpz_mul(result, sqrt, sqrt);
    if (mpz_cmp_ui(umod, 0) == 0) {
        mpz_mod(result, result, mod);
        return;
    }
    mpz_mul(result, a, result);
    mpz_mod(result, result, mod);
}

struct prev_cur {
    mpz_t prev;
    mpz_t cur;
};

bool is_prime(gmp_randstate_t& state, mpz_t& n, int steps) {
    mpz_t n1;
    mpz_init(n1);
    mpz_sub_ui(n1, n, 1);
    while (steps > 0) {
        mpz_t a;
        mpz_init(a);
        mpz_urandomm(a, state, n);
        if (mpz_cmp_ui(a, 2) < 0) {
            continue;
        }
        steps -= 1;
        int t = 0;
        mpz_t u;
        mpz_init(u);
        mpz_sub_ui(u, n, 1);
        mpz_t mod;
        mpz_init(mod);
        mpz_mod_ui(mod, u, 2);
        while (mpz_cmp(mod, u) == 0) {
            mpz_fdiv_q_ui(u, u, 2);
            mpz_mod_ui(mod, u, 2);
            t += 1;
        }
        mpz_t fst;
        mpz_init(fst);
        fast_power_mod(a, u, n, fst);
        prev_cur x;
        mpz_init(x.prev);
        mpz_set(x.prev, fst);
        mpz_init(x.cur);
        mpz_set_ui(x.cur, 0);
        for (int i = 1; i < t + 1; i++) {
            mpz_mul(x.cur, x.prev, x.prev);
            mpz_mod(x.cur, x.cur, n);
            if (mpz_cmp_ui(x.cur, 1) == 0 && mpz_cmp_ui(x.prev, 1) != 0 && mpz_cmp(x.prev, n1) != 0) {
                return false;
            }
            mpz_swap(x.cur, x.prev);
        }
        if (mpz_cmp_ui(x.prev, 1) != 0) {
            return false;
        }
    }
    return true;
}

long long log(mpz_t& C, long long num) {
    mpz_t to_log;
    mpz_init(to_log);
    mpz_sqrt(to_log, C);
    mpz_mul_ui(to_log, to_log, 2);
    mpz_add_ui(to_log, to_log, 1);
    mpz_add(to_log, to_log, C);
    long long deg = 0;
    mpz_t one;
    mpz_init_set_ui(one, 1);
    while (mpz_cmp(to_log, one) > 0) {
        deg += 1;
        mpz_mul_ui(one, one, num);
    }
    return deg;
}