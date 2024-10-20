#include <gmpxx.h>
#include <Ecm.cpp>
#include <random>
#include <iostream> 
#include <set>
// #include <functional>
#include <string>
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

// template <typename T>
// void consider_option(
//     std::string& option, 
//     std::string& current, 
//     std::string& val,
//     T& where, 
//     std::function<void(T&, T)> set, 
//     std::function<T(std::string&)> convert
// ) {
//     if (option == current) {
//         T val1 = convert(val);
//         set(where, val1);
//     }
// }

int main(int argc, char* argv[]) {
    auto ecm = Lenstra_ECM();
    int n_threads = 8;
    int B = 200000;
    mpz_class C = mpz_class(10000000);
    mpz_class n = mpz_class("235928351012000155533311111133333333");
    // mpz_class n = mpz_class("258201039002499283020763059998770852617721519");
    // mpz_class n = mpz_class("29582395193111127373737311121212121");
    // mpz_class n = mpz_class("295823951931111273737373111212121211221");
    // mpz_class n = mpz_class("295823951931111273737373111212121211221212121");
    // mpz_class n = mpz_class("2958239519311112737373731112121212112212121212121");
    // mpz_class n = mpz_class("973895488946806697397832708035027090186481661");
    // mpz_class n = mpz_class("102401151515151515151515151515151515151515151515");
    // mpz_class n = mpz_class("10240115151515151515151515151515151515151515151515");
    // mpz_class n = mpz_class("1024011515151515151515151515151515151515151515151515");
    // mpz_class n = mpz_class("1024011515151515151515151515151515151515151515151515151515");
    bool fully = false;
    for (int i = 1; i < argc; i++) {
        if (std::strcmp(argv[i], "--threads") == 0 && i < argc - 1) {
            try {
                n_threads = std::stoi(argv[i + 1]);
            } catch (std::invalid_argument& e) {
                std::cout << "Wrong input of --threads option" << "\n";
                return 1;
            }
        }
        if (std::strcmp(argv[i], "--B") == 0 && i < argc - 1) {
            try {
                B = std::stoi(argv[i + 1]);
            } catch (std::invalid_argument& e) {
                std::cout << "Wrong input of --B option" << "\n";
                return 1;
            }
        }
        if (std::strcmp(argv[i], "--С") == 0 && i < argc - 1) {
            try {
                C = mpz_class(argv[i + 1]);
            } catch (std::invalid_argument& e) {
                std::cout << "Wrong input of --C option" << "\n";
                return 1;
            }
        }
        if (std::strcmp(argv[i], "--n") == 0 && i < argc - 1) {
            try {
                n = mpz_class(argv[i + 1]);
            } catch (std::invalid_argument& e) {
                std::cout << "Wrong input of --n option" << "\n";
                return 1;
            }
        }
        if (std::strcmp(argv[i], "--fully") == 0) {
            fully = true;
        }
    }
    if (!fully) {
        mpz_class result = factor(ecm, n_threads, C, B, n);
        if (result == -1) {
            std::cout << n << "\n";
        } else {
            std::cout << result << "\n";
        }
        return 0;
    }
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