#include <Ecm.hpp>
#include <cmath>
#include <iostream>

Lenstra_ECM::Lenstra_ECM() {
    complete = false;
    primes = {2, 3, 5};
    mpz_init_set_ui(n, 4);
    mpz_init_set_ui(result, 2);
    B = 5;
}

void Lenstra_ECM::set_uncompletness() {
    complete = false;
}

mpz_t* Lenstra_ECM::get_result() {
    return &result;
}

bool Lenstra_ECM::is_not_satisfied(gmp_randstate_t& state, mpz_t& m) {
    if (is_prime(state, m, 100)) {
        mpz_init_set_si(result, -1);
        return true;
    }
    mpz_t mod;
    mpz_init(mod);
    if (mpz_mod_ui(mod, m, 2) == 0) {
        mpz_init_set_si(result, 2);
        return true;
    }
    if (mpz_mod_ui(mod, m, 3) == 0) {
        mpz_init_set_si(result, 3);
        return true;
    }
    return false;
}

bool Lenstra_ECM::multiply_point(
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
) {
    bool o1 = true;
    bool o2 = false;
    mpz_set(pow_x, x);
    mpz_set(pow_y, y);
    for (long long i = 0; i < steps + 1; i++) {
        if (num % 2 == 1) {
            if (add_points(n, a, x, y, o1, pow_x, pow_y, o2, inverse, gcd, inv, input, alpha, beta, x11, y11, g, r, s, t, q, u, v, w)) {
                return true;
            }
        }
        num /= 2;
        if (add_points(n, a, pow_x, pow_y, o2, pow_x, pow_y, o2, inverse, gcd, inv, input, alpha, beta, x11, y11, g, r, s, t, q, u, v, w)) {
                return true;
        }
    }
    return false;
}

void normalize(mpz_t& x, mpz_t& n) {
    if (mpz_cmp_si(x, 0) < 0) {
        mpz_add(x, x, n);
    }
}

bool Lenstra_ECM::add_points(
    mpz_t& n, 
    mpz_t& a, 
    mpz_t& x1,
    mpz_t& y1,
    bool& o1,
    mpz_t& x2,
    mpz_t& y2,
    bool& o2,
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
) {
    //if (mpz_cmp_si(y1, 0) != 0) gmp_printf("x1 y1 x2 y2 :  %Zd %Zd %Zd %Zd\n", x1, y1, x2, y2);
    if (o1) {
        mpz_set(x1, x2);
        mpz_set(y1, y2);
        o1 = o2;
        return false;
    } else if (o2) {
        o2 = o1;
        return false;
    }
    if (mpz_cmp(x1, x2) == 0) {
        if (mpz_cmp(y1, y2) != 0) {
            o1 = true;
            return false;
        } else {
            if (mpz_cmp_si(y1, 0) == 0) {
                o1 = true;
                return false;
            } else {
                mpz_mul_si(input, y1, 2);
                extended_euclid(input, n, gcd, inv, x11, y11, g, r, s, t, q, u, v, w);
                if (mpz_cmp_si(gcd, 1) != 0 && mpz_cmp(gcd, n) != 0) {
                    mpz_set(result, gcd);
                    gmp_printf("%Zd\n", result);
                    return true;
                } 
                mpz_set(inverse, inv);
                mpz_set_si(alpha, 3);
                mpz_mul(alpha, alpha, x1);
                mpz_mul(alpha, alpha, x1);
                mpz_add(alpha, alpha, a);
                mpz_mul(alpha, alpha, inverse);
                mpz_mod(alpha, alpha, n);
            }
        }
    } else {
        mpz_sub(input, x1, x2);
        //normalize(input, n);
        extended_euclid(input, n, gcd, inv, x11, y11, g, r, s, t, q, u, v, w);
        if (mpz_cmp_si(gcd, 1) != 0 && mpz_cmp(gcd, n) != 0) {
            mpz_set(result, gcd);
            gmp_printf("%Zd\n", result);
            return true;
        } 
        mpz_set(inverse, inv);
        mpz_sub(alpha, y1, y2);
        mpz_mul(alpha, alpha, inverse);
        mpz_mod(alpha, alpha, n);
    }
    //normalize(alpha, n);
    mpz_set(beta, y1);
    mpz_submul(beta, alpha, x1);
    mpz_mod(beta, beta, n);
    mpz_t x1s, x2s, ns, y1s, alphas, betas;
    mpz_init_set(x1s, x1);
    mpz_init_set(x2s, x2);
    mpz_init_set(y1s, y1);
    mpz_init_set(ns, n);
    mpz_init_set(alphas, alpha);
    mpz_init_set(betas, beta);
    //normalize(beta, n);
    mpz_neg(x1s, x1s);
    mpz_sub(x1s, x1s, x2s);
    mpz_addmul(x1s, alpha, alpha);
    mpz_mod(x1s, x1s, n);
    mpz_set(x1, x1s);
    //normalize(x1, n);
    mpz_set(y1s, betas);
    mpz_addmul(y1s, alphas, x1s);
    mpz_neg(y1s, y1s);
    mpz_mod(y1s, y1s, n);
    mpz_set(y1, y1s);
    //normalize(y1, n);
    o1 = false;
    //gmp_printf("%Zd %Zd %Zd %Zd\n", alpha, beta, x1, y1);
    mpz_clear(x1s);
    mpz_clear(x2s);
    mpz_clear(ns);
    mpz_clear(y1s);
    mpz_clear(alphas);
    mpz_clear(betas);
    return false;
}

bool Lenstra_ECM::try_ecm(
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
) {
    for (int prime : primes) {
        if (complete) {
            return true;
        }
        if (prime > B) {
            break;
        }
        long long steps = log(C, prime);
        //printf("%d\n", steps);
        if (prime < 100) {
            //gmp_printf("Steps: %Zd %d %d\n", C, prime, steps);
        }
        long long steps1 = static_cast<long long>(std::log2(static_cast<double>(prime)));
        for (int j = 0; j < steps; j++) {
            if (multiply_point(n, a, prime, x, y, pow_x, pow_y, inverse, gcd, inv, input, alpha, beta, x11, y11, g, r, s, t, q, u, v, w, steps1)) {
                return true;
            }
        }
    }
    return false;
}

void Lenstra_ECM::factor(mpz_t& n, int B, mpz_t& C) {
    std::mt19937 mt{std::random_device{}()};
    gmp_randstate_t state;
    gmp_randinit_mt(state);
    gmp_randseed_ui(state, mt());
    if (is_not_satisfied(state, n)) {
        return;
    }
    mtx.lock();
    if (B > this->B) {
        this->B = B;
        std::vector<int> prims = {};
        eratosthenes(prims, B);
        primes = prims;
    }
    mtx.unlock();
    mpz_t a;
    mpz_t u1;
    mpz_t v1;
    mpz_t pow_x;
    mpz_t pow_y;
    mpz_t inverse;
    mpz_t gcd;
    mpz_t inv;
    mpz_t input;
    mpz_t alpha;
    mpz_t beta;
    mpz_init(a);
    mpz_init(u1);
    mpz_init(v1);
    mpz_init(pow_x);
    mpz_init(pow_y);
    mpz_init(inverse);
    mpz_init(gcd);
    mpz_init(inv);
    mpz_init(input);
    mpz_init(alpha);
    mpz_init(beta);

    mpz_t x11;
    mpz_t y11;
    mpz_t g;
    mpz_t r;
    mpz_t s;
    mpz_t t;
    mpz_t q;
    mpz_t u;
    mpz_t v;
    mpz_t w;
    mpz_init(x11);
    mpz_init(y11);
    mpz_init(g);
    mpz_init(r);
    mpz_init(s);
    mpz_init(t);
    mpz_init(q);
    mpz_init(u);
    mpz_init(v);
    mpz_init(w);
    int k = 0;
    while (true) { 
        k += 1;
        if (complete) {
            break;
        }
        mpz_urandomm(a, state, n);
        mpz_urandomm(u1, state, n);
        mpz_urandomm(v1, state, n);
        // mpz_init_set_str(a, "100", 10);
        // mpz_init_set_str(u1, "200", 10);
        // mpz_init_set_str(v1, "300", 10);
        //gmp_printf("Randoms: %Zd %Zd %Zd\n", a, u, v);
        // mpz_class b = (v * v - u * u * u - a * u) % n;
        // normalize(b, n);
        // mpz_class D = (4 * a * a * a + 27 * b * b) % n;
        // std::tuple<mpz_class, mpz_class, mpz_class> gcds = extended_euclid(D, n);
        // mpz_class gcd = std::get<0>(gcds);
        // if (gcd > 1 && gcd < n) {
        //     result = gcd;
        //     complete = true;
        //     break;
        // } else if (gcd == n) {
        //     continue;
        // }
        if (try_ecm(n, a, B, C, u1, v1, pow_x, pow_y, inverse, gcd, inv, input, alpha, beta, x11, y11, g, r, s, t, q, u, v, w)) {
            complete = true;
            break;
        }
    }
    printf("Curves: %d\n", k);
}