// #include <gmp.h>
#include <Ecm.cpp>
#include <random>
#include <iostream> 
#include <set>
// #include <functional>
#include <string>
#include <algorithm>
#include <thread> 

// int main(int argc, char* argv[]) {
//     int n_threads = 8;
//     int B = 2000;
//     mpz_t C;
//     mpz_init_set_str(C, "100000", 10);
    // mpz_t n;
    // mpz_init_set_str(n, "161", 10);
    // mpz_t m;
    // mpz_init_set_str(m, "239", 10);
    // gcds result;
    // mpz_init(result.g);
    // mpz_init(result.x);
    // mpz_init(result.y);
    // extended_euclid(m, n, result);
    // gmp_printf("%Zd %Zd %Zd \n", result.g, result.x, result.y);
    // gmp_randstate_t state;
    // gmp_randinit_default(state);
    // gmp_randseed(state, m);
    // mpz_t p;
    // mpz_init_set_str(p, "15616700161111210676471", 10);
    // gmp_printf("%d\n", is_prime(state, p, 100));
    // mpz_t res;
    // mpz_init(res);
    // fast_power_mod(m, n, C, res);
    // gmp_printf("%Zd\n", res);

    // mpz_t mm;
    // mpz_init_set_str(mm, "17", 10);
    // mpz_t tt;
    // mpz_init_set_str(tt, "29", 10);
    // mpz_t cc;
    // mpz_init(cc);
    // mpz_mod(cc, mm, tt);
    // gmp_printf("%Zd\n", cc);
    // mpz_t num;
    // mpz_t x;
    //mpz_init_set_str(x, "100", 10);
    // mpz_init_set_str(num, "12327723467569", 10);
    //mpz_init_set_str(num, "22121212111111111111111111111111111", 10);
    // mpz_init_set_str(num, "216295719507710942400449895291244347133", 10);
    // Lenstra_ECM ecm = Lenstra_ECM();
    // ecm.factor(num, B, C);
    // gmp_printf("Result: %Zd\n", *ecm.get_result());

    // mpz_t x11;
    // mpz_t y11;
    // mpz_t g;
    // mpz_t r;
    // mpz_t s;
    // mpz_t t;
    // mpz_t q;
    // mpz_t u;
    // mpz_t v;
    // mpz_t w;
    // mpz_t gcd;
    // mpz_t inv;
    // mpz_init(x11);
    // mpz_init(y11);
    // mpz_init(g);
    // mpz_init(r);
    // mpz_init(s);
    // mpz_init(t);
    // mpz_init(q);
    // mpz_init(u);
    // mpz_init(v);
    // mpz_init(w);
    // mpz_init(gcd);
    // mpz_init(inv);
    // extended_euclid(mm, tt, gcd, inv, x11, y11, g, r, s, t, q, u, v, w);
    // gmp_printf("See: %Zd %Zd\n", gcd, inv);
    // mpz_t x1, x2, n, y1, alpha, beta;
    // mpz_init_set_ui(x1, 200);
    // mpz_init_set_ui(x2, 200);
    // mpz_init_set_ui(y1, 300);
    // mpz_init_set_ui(n, 323);
    // mpz_init_set_ui(alpha, 254);
    // mpz_init_set_ui(beta, 211);
    // mpz_neg(x1, x1);
    // mpz_sub(x1, x1, x2);
    // mpz_addmul(x1, alpha, alpha);
    // mpz_mod(x1, x1, n);
    // //normalize(x1, n);
    // mpz_set(y1, beta);
    // mpz_addmul(y1, alpha, x1);
    // mpz_neg(y1, y1);
    // mpz_mod(y1, y1, n);
    // //normalize(y1, n);
    // gmp_printf("%Zd %Zd %Zd %Zd\n", alpha, beta, x1, y1);
// }

void solve(
    Lenstra_ECM& ecm, 
    mpz_t& C, 
    int B, 
    mpz_t& num
) {
    ecm.factor(num, B, C);
}

void factor(
    Lenstra_ECM& ecm, 
    int n_threads, 
    mpz_t& C, 
    int B, 
    mpz_t& num
) {
    std::vector<std::thread> threads;
    for (int i = 0; i < n_threads; i++) {
        std::thread th(solve, std::ref(ecm), std::ref(C), B, std::ref(num));
        threads.push_back(std::move(th));
    }
    for (int i = 0; i < n_threads; i++) {
        threads[i].join();
    }
    //ecm.set_uncompletness();
    gmp_printf("%Zd\n", *ecm.get_result());
}

int main(int argc, char* argv[]) {
    auto ecm = Lenstra_ECM();
    int n_threads = 16;
    int B = 50000;
    mpz_t C;
    mpz_init_set_ui(C, 10000000);
    mpz_t n;
    //mpz_init_set_str(n, "10000000000000000000000083000000000000000000000091", 10);
    //mpz_init_set_str(n, "2074757784440496479256203931845580575506223116121218449997828664845326405706454073199853524473551897144098943305650394591197575537705887653943437417056981843530590901700754761842687", 10);
    //mpz_init_set_str(n, "10000000089000000133", 10);
    //mpz_init_set_str(n, "10000000000000431000000000002257", 10);
    //mpz_init_set_str(n, "100000000000000003300000000000000009", 10);
    //mpz_init_set_str(n, "86347590191400517636741183471745339", 10);
    //mpz_init_set_str(n, "11631287175507576935276170884437407075704861820681", 10);
    //mpz_init_set_str(n, "100000000000000000050700000000000000004563", 10);
    //mpz_init_set_str(n, "10000000000000000081000000000000000153", 10);
    mpz_init_set_str(n, "107095543673036920778818309058867449715781454676669646697666415811636161270588853460929436069818260743444352767", 10);
    factor(ecm, n_threads, C, B, n);
}
