#include <Point.cpp>
#include <variant>

class Lenstra_ECM {
public:
    Lenstra_ECM() = default;
    void factor(mpz_class& n, int B, mpz_class& C);
    mpz_class get_result();
private:
    std::vector<int> primes = {2, 3, 5};
    mpz_class n = 4;
    mpz_class result = 2;
    int B = 5;

    bool is_not_satisfied(gmp_randclass& rr, mpz_class& m);
    std::variant<bool, Point> add_points(
        mpz_class& n, 
        mpz_class& a, 
        std::variant<bool, Point>& point1, 
        std::variant<bool, Point>& point2
    );
    std::variant<bool, Point> multiply_point(
        mpz_class& n, 
        mpz_class& a, 
        long long num, 
        std::variant<bool, Point>& point
    );
    bool try_ecm(
        mpz_class& n, 
        mpz_class& a, 
        int B, 
        mpz_class& C, 
        Point& point
    );
};