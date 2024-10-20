#include <gmpxx.h>
#include <Ecm.cpp>
#include <random>
#include <iostream> 
#include <set>
#include <thread> 

void solve(
    Lenstra_ECM& ecm, 
    mpz_class& C, 
    int B, 
    mpz_class& num
) {
    ecm.factor(num, B, C);
}

mpz_class factor(
    Lenstra_ECM& ecm, 
    int n_threads, 
    mpz_class& C, 
    int B, 
    mpz_class& num
) {
    std::vector<std::thread> threads;
    for (int i = 0; i < n_threads; i++) {
        std::thread th(solve, std::ref(ecm), std::ref(C), B, std::ref(num));
        threads.push_back(std::move(th));
    }
    for (int i = 0; i < n_threads; i++) {
        threads[i].join();
    }
    ecm.set_uncompletness();
    return ecm.get_result();
}

std::multiset<mpz_class> factor_fully(
    Lenstra_ECM& ecm, 
    int n_threads,
    mpz_class& C, 
    int B, 
    mpz_class& n
) {
    mpz_class num = mpz_class(n);
    std::multiset<mpz_class> mset;
    while (num % 2 == 0) {
        mset.insert(2);
        num /= 2;
    }
    while (num % 3 == 0) {
        mset.insert(3);
        num /= 3;
    }
    while (num != 1) {
        mpz_class non_tr = factor(ecm, n_threads, C, B, num);
        if (non_tr == -1) {
            mset.insert(num);
            break;
        }
        mpz_class d = factor(ecm, n_threads, C, B, non_tr);
        if (d == -1) {
            while (num % non_tr == 0) {
                mset.insert(non_tr);
                num /= non_tr;
            }
        } else {
            std::multiset<mpz_class> mset1 = factor_fully(ecm, n_threads, C, B, non_tr);
            mset.insert(mset1.begin(), mset1.end());
            num /= non_tr;
        }
    }
    return mset;
}

int main() {
    auto ecm = Lenstra_ECM();
    int B = 200000;
    mpz_class C = mpz_class(10000000);
    // mpz_class n = mpz_class("258201039002499283020763059998770852617721519");
    // mpz_class n = mpz_class("235928351012000155533311111133333333");
    std::multiset<mpz_class> prs = factor_fully(ecm, 16, C, B, n);
    mpz_class pr = *prs.begin();
    int deg = 0;
    bool need = false;
    for (auto it = prs.begin(); it != prs.end(); it++) {
        if (pr == *it) {
            deg++;
        } else {
            if (need) {
                std :: cout << " * ";
            }
            if (deg > 1) {
                std::cout << pr << "^" << prs.count(pr);
            } else {
                std::cout << pr;
            }
            pr = *it;
            deg = 1;
            need = true;
        }
    }
    if (need) {
        std :: cout << " * ";
    }
    if (deg > 1) {
        std::cout << pr << "^" << prs.count(pr) << "\n";
    } else {
        std::cout << pr << "\n";
    }
}
