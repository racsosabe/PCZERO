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
#include "ExponentialPolynomial.hpp"

using my_complex = boost::multiprecision::cpp_complex<334, boost::multiprecision::backends::digit_base_2, void, std::int32_t, -262142, 262143>;
using my_real = boost::multiprecision::number<boost::multiprecision::backends::cpp_bin_float<334, boost::multiprecision::backends::digit_base_2, void, std::int32_t, -262142, 262143>, boost::multiprecision::et_off>;

using polynomial = ExponentialPolynomial<my_real, my_complex>;

int main(int argc, char* argv[]){
    if(argc != 5) {
        std::cout << "Usage:" << std::endl;
        std::cout << "   ./a.out (lower-left real) (lower-left imag) (upper-right real) (upper-right imag)" << std::endl << std::endl;
        std::cout << "Note: You must use a polynomial function library that supports:" << std::endl;
        std::cout << " - Evaluation on a complex point" << std::endl;
        std::cout << " - Computing the derivative of the function" << std::endl;
        std::cout << " - Computing the M bound according to the paper" << std::endl;
        std::cout << "In addition, you must implement the function \"base_function()\" that returns f" << std::endl;
        return 0;
    }
    polynomial f = base_function();
    polynomial df = f.derivative();
    polynomial M = f.Mbound();
    polynomial N = M;
    pczero<my_complex, my_real, 8, 1> Solver(f, df, M, N);
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
