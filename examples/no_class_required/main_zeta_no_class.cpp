#pragma GCC optimize ("O3")
#pragma GCC target ("avx2")
#include<string>
#include<iostream>
#include<stack>
#include<utility>
#include<fstream>
#include<sstream>
#include<string>
#include<tuple>
#include<chrono>
#include<thread>
#include<mutex>
#include<unistd.h>
#include<ctime>
#include<set>
#include<boost/multiprecision/cpp_bin_float.hpp>
#include<boost/multiprecision/cpp_complex.hpp>
#include<boost/multiprecision/complex_adaptor.hpp>
#include "PCZERO.hpp"

using my_complex = boost::multiprecision::cpp_complex<334, boost::multiprecision::backends::digit_base_2, void, std::int32_t, -262142, 262143>;
using my_real = boost::multiprecision::number<boost::multiprecision::backends::cpp_bin_float<334, boost::multiprecision::backends::digit_base_2, void, std::int32_t, -262142, 262143>, boost::multiprecision::et_off>;

const int MAX = 41;

int len;                 // n
std::array<my_real, MAX> f_freqs;      // l_k is freqs[k]
std::array<my_real, MAX> f_coeff;      // a_k is coeffs[k]

std::array<my_real, MAX> df_freqs;      // l_k is freqs[k]
std::array<my_real, MAX> df_coeff;      // a_k is coeffs[k]

std::array<my_real, MAX> M_freqs;      // l_k is freqs[k]
std::array<my_real, MAX> M_coeff;      // a_k is coeffs[k]

void init(int n) {
    len = n;
    f_freqs[0] = 0;
    f_coeff[0] = 1;
    for(int k = 1; k < len; k++) {
        f_freqs[k] = log(k + 1);
        f_coeff[k] = 1;
    }
    for(int k = 0; k < len; k++) {
        df_freqs[k] = f_freqs[k];
        df_coeff[k] = -f_coeff[k] * f_freqs[k];
    }
    for(int k = 0; k < len; k++) {
        M_freqs[k] = f_freqs[k];
        M_coeff[k] = abs(f_coeff[k]) * f_freqs[k] * f_freqs[k];
    }
}

my_complex f(my_complex x) {
    my_complex ans{0, 0};
    for(int i = 0; i < len; i++){
        ans += f_coeff[i] * exp(-f_freqs[i] * x);
    }
    return ans;
}

my_complex df(my_complex x) {
    my_complex ans{0, 0};
    for(int i = 0; i < len; i++){
        ans += df_coeff[i] * exp(-df_freqs[i] * x);
    }
    return ans;
}

my_real M(my_complex s1, my_complex s2) {
    my_real x = s1.real() < s2.real() ? s1.real() : s2.real();
    my_complex ans{0, 0};
    for(int i = 0; i < len; i++){
        ans += f_coeff[i] * exp(-f_freqs[i] * x);
    }
    return abs(ans);
}

my_real N(my_complex s1, my_complex s2) {
    return M(s1, s2);
}

int main(int argc, char* argv[]){
    if(argc != 6) {
        std::cout << "Usage:" << std::endl;
        std::cout << "   ./a.out (lower-left real) (lower-left imag) (upper-right real) (upper-right imag) (n)" << std::endl << std::endl;
        std::cout << "Note: You must use a polynomial function library that supports:" << std::endl;
        std::cout << " - Evaluation on a complex point" << std::endl;
        std::cout << " - Computing the derivative of the function" << std::endl;
        std::cout << " - Computing the M bound according to the paper" << std::endl;
        std::cout << "In addition, you must implement the function \"base_function()\" that returns f" << std::endl;
        return 0;
    }
    int n = atoi(argv[5]);
    init(n);
    pczero<my_complex, my_real, 8, 0> Solver(f, df, M, N);
    auto start = std::chrono::high_resolution_clock::now();
    my_complex LD{atof(argv[1]), atof(argv[2])}, RU{atof(argv[3]), atof(argv[4])};
    Solver.solve(LD, RU);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    Solver.print_final_information(duration);
    Solver.download();
//	pczero::test_newton();
    return 0;
}
