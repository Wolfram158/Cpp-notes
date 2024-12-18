#include <Ecm.cpp>
#include <random>
#include <iostream> 
#include <set>
#include <string>
#include <algorithm>
#include <thread> 
#include <fstream>

void solve(
    Lenstra_ECM& ecm, 
    mpz_t& C, 
    int B, 
    mpz_t& num
) {
    ecm.factor(num, B, C);
}

mpz_t* factor(
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
    ecm.set_uncompletness();
    return ecm.get_result();
}

std::multiset<mpz_class> factor_fully(
    Lenstra_ECM& ecm, 
    int n_threads,
    mpz_t& C, 
    int B, 
    mpz_class& n
) {
    mpz_class num = mpz_class(n);
    mpz_t num1;
    mpz_t non_tr1;
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
        mpz_init_set(num1, num.get_mpz_t());
        mpz_class non_tr = mpz_class(*factor(ecm, n_threads, C, B, num1));
        if (non_tr == -1) {
            mset.insert(num);
            break;
        }
        mpz_init_set(non_tr1, non_tr.get_mpz_t());
        mpz_class d = mpz_class(*factor(ecm, n_threads, C, B, non_tr1));
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

void add_examples(std::vector<mpz_class>& examples) {
    examples.push_back(mpz_class("235928351012000155533311111133333333"));
    examples.push_back(mpz_class("258201039002499283020763059998770852617721519"));
    examples.push_back(mpz_class("29582395193111127373737311121212121"));
    examples.push_back(mpz_class("295823951931111273737373111212121211221"));
    examples.push_back(mpz_class("295823951931111273737373111212121211221212121"));
    examples.push_back(mpz_class("2958239519311112737373731112121212112212121212121"));
    examples.push_back(mpz_class("973895488946806697397832708035027090186481661"));
    examples.push_back(mpz_class("102401151515151515151515151515151515151515151515"));
    examples.push_back(mpz_class("10240115151515151515151515151515151515151515151515"));
    examples.push_back(mpz_class("1024011515151515151515151515151515151515151515151515"));
    examples.push_back(mpz_class("1024011515151515151515151515151515151515151515151515151515"));
    examples.push_back(mpz_class("15151515151515151515151515151515151515"));
    examples.push_back(mpz_class("1515151515151515151515151515151515151515151515"));
    examples.push_back(mpz_class("316343013244102577022636990019897602752593577"));
    examples.push_back(mpz_class("31634331111111111111111111111111111"));
    examples.push_back(mpz_class("141313511251251251251251522111123411114212211112221"));
    examples.push_back(mpz_class("14131351125125125125125152211112341111421221111222121"));
    examples.push_back(mpz_class("141313511251251251251251522111123411114212211112221211111111"));
    examples.push_back(mpz_class("141313511251251253222222222222000000222211111111111111111111"));
    examples.push_back(mpz_class("12521441122222222222222222222222222222222222222222222222"));
    examples.push_back(mpz_class("1204124151111111111111111111111111111111111111111111111111"));
    examples.push_back(mpz_class("1010101010101010101010101010101010110101010101010101010101001010110111111"));
    examples.push_back(mpz_class("10101010101010101010101010101010101101010101010101010101010010101101111"));
    std::sort(examples.begin(), examples.end());
}

int main(int argc, char* argv[]) {
    auto ecm = Lenstra_ECM();
    int n_threads = 16;
    int B = 70000;
    mpz_class C = mpz_class(100000);
    mpz_class n = mpz_class("235928351012000155533311111133333333");
    bool fully = false;
    bool example = false;
    bool file = false;
    std::string file_name{""};
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
        if (std::strcmp(argv[i], "--C") == 0 && i < argc - 1) {
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
        if (std::strcmp(argv[i], "--file") == 0 && i < argc - 1) {
            try {
                file = true;
                file_name = std::string(argv[i + 1]);
            } catch (std::invalid_argument& e) {
                std::cout << "Wrong input of --file option" << "\n";
                return 1;
            }
        }
        if (std::strcmp(argv[i], "--fully") == 0) {
            fully = true;
        }
        if (std::strcmp(argv[i], "--example") == 0) {
            example = true;
        }
    }
    if (file) {
        fully = true;
        std::ifstream fin(file_name);
        fin >> n;
        fin >> B;
        fin >> C;
    } else if (example) {
        std::vector<mpz_class> examples;
        add_examples(examples);
        while (true) {
            std::cout << "Choose one of the following examples:\n";
            for (std::vector<mpz_class>::size_type i = 1; i < examples.size() + 1; i++) {
                std::cout << "(" << i << ")" << " " << examples[i - 1] << "\n";
            }
            std::vector<mpz_class>::size_type where;
            std::cin >> where;
            if (where >= 0 && where - 1 < examples.size()) {
                n = examples[where - 1];
                break;
            }
        }
    }
    mpz_t C1;
    mpz_init_set(C1, C.get_mpz_t());
    if (!fully) {
        mpz_t n1;
        mpz_init_set(n1, n.get_mpz_t());
        mpz_class result = mpz_class(*factor(ecm, n_threads, C1, B, n1));
        if (result == -1) {
            std::cout << n << "\n";
        } else {
            std::cout << result << "\n";
        }
        return 0;
    }
    std::multiset<mpz_class> prs = factor_fully(ecm, n_threads, C1, B, n);
    mpz_class pr = *prs.begin();
    int deg = 0;
    bool need = false;
    gmp_printf("%Zd = ", n);
    for (auto it = prs.begin(); it != prs.end(); it++) {
        if (pr == *it) {
            deg++;
        } else {
            if (need) {
                std::cout << " * ";
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
        std::cout << " * ";
    }
    if (deg > 1) {
        std::cout << pr << "^" << prs.count(pr) << "\n";
    } else {
        std::cout << pr << "\n";
    }
}