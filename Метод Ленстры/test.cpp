#include <gmpxx.h>
#include <random>
#include <iostream> 
#include <Ecm.cpp>

int main() {
  // mpz_t x;
  // mpz_t y;
  // mpz_init_set_si(x, 1);
  // mpz_init_set_si(y, 1);
  // for (int i = 1; i < 10; i++) {
  //   mpz_t num;
  //   mpz_init_set_si(num, i);
  //   mpz_mul(x, x, num);
  // }
  // mpz_out_str(stdout, 10, x);
    // mpz_class ran;
    // gmp_randclass rr(gmp_randinit_mt);
    // std::mt19937 mt{ std::random_device{}() };
    // rr.seed(mt());
    // mpz_class x = mpz_class("3523523523235235333346274267272472423616616116161412414274");
    // for (int i = 2000; i < 2100; i++) {
    //   ran = rr.get_z_range(x);
    //   //std::cout << (ran + 1) % 43452 << "\n";
    //   mpz_class num = mpz_class(i); 
    //   std::cout << i << " " << is_prime(rr, num, 100) << "\n";
    // }
    // Point p = Point(mpz_class(3), mpz_class(5), true);
    // std::cout << p.get_x() << " " << p.get_y() << " " << p.get_is_o() << "\n";
    // mpz_class z = mpz_class("9099999999999999999099999999999999999999");
    // mpf_class n = mpf_class(z, 1000);
    // mpf_class u = floor(sqrt(n));
    // u.set_prec(100000);
    // mpz_class z = mpz_class(5);
    // mpz_class t = mpz_class(3);
    // auto res = extended_euclid(z, t);
    // std::cout << std::get<0>(res) << " " << std::get<1>(res) << " " << std::get<2>(res) << "\n";
    auto ecm = Lenstra_ECM();
    mpz_class n = mpz_class("87178204020795979291199");
    mpz_class C = mpz_class("1000000000");
    ecm.factor(n, 100000, C);
    std::cout << ecm.get_result() << "\n";
}
