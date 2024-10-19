#include <Ecm.hpp>
#include <Utils.cpp>
#include <cmath>

bool Lenstra_ECM::is_not_satisfied(mpz_class& m) {
    if (is_prime(m)) {
        this->result = -1;
        return true;
    }
    if (m % 2 == 0) {
        this->result = 2;
        return true;
    }
    if (m % 2 == 0) {
        this->result = 3;
        return true;
    }
    return false;
}

std::variant<bool, Point> Lenstra_ECM::multiply_point(
    mpz_class& n, 
    mpz_class& a, 
    long long num, 
    std::variant<bool, Point>& point
) {
    long long steps = (long long) std::log2((double) num);
    std::variant<bool, Point> result = Point(0, 0, true);
    std::variant<bool, Point> current_pow = point;
    for (int i = 0; i < steps + 1; i++) {
        if (num % 2 == 0) {
            result = this->add_points(n, a, result, current_pow);
            if (holds_alternative<bool>(result)) {
                return false;
            }
        }
        num /= 2;
        current_pow = this->add_points(n, a, current_pow, current_pow);
        if (holds_alternative<bool>(result)) {
                return false;
        }
    }
    return result;
}

std::variant<bool, Point> Lenstra_ECM::add_points(
    mpz_class& n, 
    mpz_class& a, 
    std::variant<bool, Point>& point1, 
    std::variant<bool, Point>& point2
) {
    Point p1 = get<Point>(point1);
    Point p2 = get<Point>(point2);
    mpz_class alpha;
    bool end = false;
    if (p1.get_x() == p2.get_x()) {
        if (p1.get_y() != p2.get_y()) {
            return Point(0, 0, true);
        } else {
            if (p1.get_y() == 0) {
                return Point(0, 0, true);
            } else {
                auto input = 2 * p1.get_y();
                auto gcds = extended_euclid(input, n);
                auto gcd = std::get<0>(gcds);
                if (gcd != 1) {
                    this->result = gcd;
                    return false;
                } else {
                    auto inverse = std::get<1>(gcds);
                    alpha = ((3 * p1.get_x() * p1.get_x() + a) * inverse) % n;
                }
            }
        }
    } else {
        auto input = p1.get_x() - p2.get_x();
        auto gcds = extended_euclid(input, n);
        auto gcd = std::get<0>(gcds);
        if (gcd != 1) {
            this->result = gcd;
            return false;
        }
        auto inverse = std::get<1>(gcds);
        alpha = ((p1.get_y() - p2.get_y()) * inverse) % n;
    }
    auto beta = (p1.get_y() - alpha * p1.get_x()) % n;
    auto x = (alpha * alpha - p1.get_x() - p2.get_x()) % n;
    auto y = (-(alpha * x + beta)) % n;
    return Point(x, y);
}

bool Lenstra_ECM::try_ecm(
    mpz_class& n, 
    mpz_class& a, 
    int B, 
    long long steps, 
    Point& point
) {
    std::variant<bool, Point> Q = point;
    for (prime : this->primes) {
        if (prime > B) {
            break;
        }
        for (int j = 0; j < steps; j++) {
            this->multiply_point(n, a, prime, Q);
        }
        if (holds_alternative<bool>(Q)) {
            return true;
        }
    }
    return false;
}

void Lenstra_ECM::factor(mpz_class& n, int B, mpz_class& C) {
    if (B > this->B) {
        this->B = B;
        this->primes = eratosthenes(B);
    }
    long long steps = log(C);
    mpz_class ran;
    gmp_randclass rr(gmp_randinit_mt);
    std::mt19937 mt{std::random_device{}()};
    rr.seed(mt());
    while (true) {
        mpz_class a = rr.get_z_range(n);
        mpz_class u = rr.get_z_range(n);
        mpz_class v = rr.get_z_range(n);
        Point point = Point(u, v);
        if (this->try_ecm(n, a, B, steps, point)) {
            break;
        }
    }
}