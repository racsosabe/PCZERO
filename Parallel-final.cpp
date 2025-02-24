/*
 * File: Parallel-final.cpp
 * Author: Oswaldo Velasquez, Manuel Toribio and Victor Galvan
 * Version: 1.0.0
 * Date: October 25, 2024
 *
 * Copyright (c) 2024 Oswaldo Velasquez, Manuel Toribio and Victor Galvan
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

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
 
 typedef boost::multiprecision::cpp_complex<334, boost::multiprecision::backends::digit_base_2, void, std::int32_t, -262142, 262143> my_complex;
 typedef boost::multiprecision::number<boost::multiprecision::backends::cpp_bin_float<334, boost::multiprecision::backends::digit_base_2, void, std::int32_t, -262142, 262143>, boost::multiprecision::et_off>  my_real;
 
 const my_real PI = acos(my_real(-1));
 const my_complex NIL{3, 0};
 
 const my_real EPS = pow(my_real(0.1), 50);
 const my_real eps1 = pow(my_real(0.1), 60);
 const my_real eps2 = pow(my_real(0.1), 90);
 const my_real eps3 = 0.5;
 const int limit = 8;
 
 std::string outputname;
 std::string outputnamex;
 std::string benchmarkname;
 
 const std::string red("\033[0;31m");
 const std::string green("\033[0;32m");
 const std::string reset("\033[0m");
 
 typedef std::tuple<my_complex, my_complex, my_real, my_real, my_real, my_real> data;
 typedef std::pair<my_complex, my_real> answer;
 
 // Defining f, its derivative df and its M bound for the second derivative
 
 // MAX must be greater than the value of n
 const int MAX = 41;
 // A Dirichlet polynomial f(s) = \sum_{k=0}^n a_k e^{-l_k s}, where 
 // 0 = l_0< l_1 < l_2 < ... < l_n is represented by the following structure
 struct dirichlet{
     int len;                 // n
     std::array<my_real, MAX> freqs;      // l_k is freqs[k]
     std::array<my_real, MAX> coeff;      // a_k is coeffs[k]
 };
 
 // Computing the derivative of the above function
 dirichlet mderivative(const dirichlet &f){
     dirichlet tmp;
     int le = f.len;
     tmp.len = le;
     for(int k = 0; k < le; k++) {
         tmp.freqs[k] = f.freqs[k];
         tmp.coeff[k] = -f.coeff[k] * f.freqs[k];
     }
     return tmp;
 }
 
 // Computing the M bound for the second derivative of the function f
 dirichlet Mbound(const dirichlet &f){
     dirichlet tmp;
     int le = f.len;
     tmp.len = le;
     for(int k = 0; k < le; k++) {
         tmp.freqs[k] = f.freqs[k];
         tmp.coeff[k] = abs(f.coeff[k]) * f.freqs[k] * f.freqs[k];
     }
     return tmp;
 }
 
 // Definition of the matrix associated to the Riemann Zeta Function partial sum zeta_n = \sum_{k=0}^n k^{-s}
 // It can be customized for other Dirichlet polynomials, and with slight modifications, for exponential polynomials (with real frequencies)
 dirichlet zetan_matrix(int v){
     dirichlet tmp;
     tmp.len = v;
     tmp.freqs[0] = 0;
     tmp.coeff[0] = 1;
     for(int k = 1; k < v; k++) {
         tmp.freqs[k] = log(k + 1);
         tmp.coeff[k] = 1;
     }
     return tmp;
 }
 
 void show(const dirichlet &a){
     std::cout << "n = " << a.len << "\n";
     std::cout << "Freqs: ";
     for(int l = 0; l < a.len; l++){
         std::cout << a.freqs[l] << " \n"[l + 1 == a.len];
     }
     std::cout << "Coeffs: ";
     for(int l = 0; l < a.len; l++){
         std::cout << a.coeff[l] << " \n"[l + 1 == a.len];
     }
 }
 
 my_complex eval(const dirichlet &A, const my_complex &x){
     my_complex ans{0, 0};
     for(int i = 0; i <= A.len; i++){
         ans += A.coeff[i] * exp(-A.freqs[i] * x);
     }
     return ans;
 }
 
 void print_progress(int cur, int total){
     int frac = cur * 80 / total;
     std::cout << green;
     std::cout << "|";
     for(int i = 0; i < 80; i++) std::cout << " ="[i < frac];
     std::cout << "|\n";
     std::cout << reset;
 }
 
 int cnt = 0;
 int total_zeros;
 int excluded = 0;
 int active_threads = limit;
 
 int n;
 int successful_newtons = 0;
 const int USE_EXCLUSION_TEST = 0; // Set to 1 to use exclusion test for rectangles
 dirichlet f_matrix, df_matrix, M_matrix;
 
 struct MutexRectangles{
     std::stack<data> S;
     mutable std::mutex m;
 
     void push(data val){
         std::lock_guard<std::mutex> lock(m);
         S.emplace(val);
     }
 
     data pop(){
         std::lock_guard<std::mutex> lock(m);
         if(S.empty()){
             return data(NIL, NIL, 0, 0, 0, 0);
         }
         data res = S.top();
         S.pop();
         return res;
     }
 
     int size(){
         std::lock_guard<std::mutex> lock(m);
         return S.size();
     }
 };
 
 struct MutexAnswer{
     std::stack<answer> S;
     mutable std::mutex m;
 
     void push(answer val){
         std::lock_guard<std::mutex> lock(m);
         S.emplace(val);
         cnt++;
         if(cnt % 10 == 0){
             std::cerr << "\033[2J\033[1;1H";
             std::cerr << "Now we've got " << cnt << "/" << total_zeros << " zeros" << '\n';
             std::cerr << "Number of active threads: " << active_threads << '\n';
             std::cerr << "Excluded rectangles: " << excluded << '\n';
             auto cur = std::chrono::system_clock::now();
             std::time_t cur_time = std::chrono::system_clock::to_time_t(cur);
             print_progress(cnt, total_zeros);
             std::cerr << "Registered at " << std::ctime(&cur_time) << '\n';
         }
     }
 
     answer pop(){
         std::lock_guard<std::mutex> lock(m);
         if(S.empty()){
             return answer(NIL, 0);
         }
         answer res = S.top();
         S.pop();
         return res;
     }
 
     int size(){
         std::lock_guard<std::mutex> lock(m);
         return S.size();
     }
 };
 
 MutexRectangles S;
 MutexAnswer Answers;
 
 
 // Initialization of custom functions
 void funcinit(int argc, char* argv[]){
     n = atoi(argv[1]);
     std::cout << "Main function: \n";
     f_matrix = zetan_matrix(n);
     show(f_matrix);
     std::cout << "Derivative function: \n";
     df_matrix = mderivative(f_matrix);
     std::cout << "M bound function: \n";
     M_matrix = Mbound(f_matrix);
     outputname = "zeros" + std::to_string(n) + ".txt";
     outputnamex = "zeros" + std::to_string(n) + ".gp";
     benchmarkname = "log" + std::to_string(n) + ".txt";
 }
 
 // Starting the definition of custom functions 
 my_complex f(const my_complex &x){
     return eval(f_matrix, x);
 }
 
 my_complex df(const my_complex &x){
     return eval(df_matrix, x);
 }
 
 my_real M(const my_complex &s1, const my_complex &s2){
     my_real ans = (s1.real() < s2.real()) ? abs(eval(M_matrix, s1.real())) : abs(eval(M_matrix, s2.real()));
     return ans;
 }
 
 my_real N(const my_complex &LD, const my_complex &RU){
     return M(LD, RU);
 }
 
 
 my_real get_min(const my_complex &a, my_complex b){
     my_complex prod = a * conj(b);
     if(prod.real() >= norm(a)) return abs(a);
     if(prod.real() >= norm(b)) return abs(b);
     return abs(prod.imag() / abs(b - a));
 }
 
 
 bool enrect(my_complex s1, my_complex s2, my_complex x){
     my_real Lreal = s1.real() < s2.real() ? s1.real() : s2.real();
     my_real Rreal = s1.real() + s2.real() - Lreal;
     my_real Limag = s1.imag() < s2.imag() ? s1.imag() : s2.imag();
     my_real Rimag = s1.imag() + s2.imag() - Limag;
     return Lreal <= x.real() and x.real() <= Rreal and Limag <= x.imag() and x.imag() <= Rimag;
 }
 
 bool var_arg(const my_complex &s1, const my_complex &s2, my_real &res){
     my_complex d = s2 - s1;
     my_complex s3 = s1;
     my_real t = 1;
     my_real t0 = 0;
     res = 0;
     while(t0 + eps1 < 1){
         my_complex s4 = s3 + t * d;
         my_complex fs3 = f(s3);
         my_real R0 = M(s3, s4) * norm(s4 - s3) / 8;
         my_real minP = get_min(fs3, f(s4));
         while(minP <= R0 + EPS and t > eps2){
             t /= 2;
             R0 /= 4;
             s4 = s3 + t * d;
             minP = get_min(fs3, f(s4));
         }
         if(t <= eps2){
             return false;
         }
         res += arg(f(s4) / f(s3));
         t0 += t;
         s3 += d * t;
         t = abs((s2 - s3) / (s2 - s1)) < 2 * t ? abs((s2 - s3) / (s2 - s1)) : 2 * t;
         if(t > 1 + eps1 - t0) t = 1 + eps1 - t0;
     }
     return true;
 }
 
 my_real certification_ratio(const my_complex &z, const my_complex &LD, const my_complex &RU){
     my_real eta = abs(f(z) / df(z));
     my_real K = N(LD, RU) / abs(df(z));
     my_real h = K * eta;
     if(h < 0.5){
         return 2 * eta / (1 + sqrt(1 - 2 * h));
     }
     else return RU.real() - LD.real() > RU.imag() - LD.imag() ? RU.real() - LD.real() : RU.imag() - LD.imag();
 }
 
 bool newton(const my_complex &LD, const my_complex &RU, my_complex &res){
     if(RU.real() - LD.real() > 3 and RU.imag() - LD.imag() > 3) return false;
     my_complex x0 = (LD + RU) / 2.0;
     my_complex xi;
     int iter = 10;
     do{
         xi = x0;
         x0 = xi - f(xi) / df(xi);
     }while(--iter > 0 and abs(x0 - xi) > eps2 and enrect(LD, RU, x0));
     if(!enrect(LD, RU, x0) or abs(x0 - xi) > eps2) return false;
     res = x0;
     return true;
 }
 
 bool exclusion_test(my_complex &LD, my_complex &RU){
     my_complex z0 = (LD + RU) / 2;
     my_real r = abs(RU - LD) / 2;
     my_real upper_bound = N(LD, RU);
     my_real val = abs(f(z0)) - abs(df(z0)) * r - upper_bound * r * r / 2;
     return val > eps2;
 }
 
 void add_answer(my_complex &s, my_complex &LD, my_complex &RU, int id){
     Answers.push(std::make_pair(s, certification_ratio(s, LD, RU)));
 }
 
 void process_one(std::tuple<my_complex, my_complex, my_real, my_real, my_real, my_real> &val, int id){
     my_complex LD, RU;
     my_real V0, V1, V2, V3;
     std::tie(LD, RU, V0, V1, V2, V3) = val;
     my_real zeros = round((V0 + V1 + V2 + V3) / 2 / PI);
     my_real dx = (RU.real() - LD.real()) * 0.01;
     my_real dy = (RU.imag() - LD.imag()) * 0.01;
     if(dx > eps3) dx = eps3;
     if(dy > eps3) dy = eps3;
     if(fabs(zeros - 1) < EPS){
         my_complex res;
         if(newton(LD, RU, res)){
             add_answer(res, LD, RU, id);
             return;
         }
     }
     if(RU.real() - LD.real() > RU.imag() - LD.imag()){
         my_complex M1{(LD.real() + RU.real()) / 2, LD.imag()};
         my_complex M2{(LD.real() + RU.real()) / 2, RU.imag()};
         my_real FM;
         while(not var_arg(M1, M2, FM)){
             M1 += dx;
             M2 += dx;
         }
         my_real r = abs(LD - RU) / 2;
         my_real FD = 0, FU = 0;
         bool exclude_left = USE_EXCLUSION_TEST? exclusion_test(LD, M2) : false;
         if(!exclude_left) var_arg(LD, M1, FD);
         bool exclude_right = USE_EXCLUSION_TEST ? exclusion_test(M1, RU) : false;
         if(!exclude_right) var_arg(RU, M2, FU);
         if(!USE_EXCLUSION_TEST or (!exclude_left and !exclude_right)){
             my_real v1 = FD + FM + V2 - FU + V3;
             my_real v2 = V0 - FD + V1 + FU - FM;
             if(abs(v1 / 2 / PI) < 0.1) exclude_left = true;
             if(abs(v2 / 2 / PI) < 0.1) exclude_right = true;
             if(abs(round((v1 + v2) / 2 / PI) - zeros) > 0.1){
                 std::cerr << "Fatal error. Sum of partitions isn't equal to total" << std::endl;
                 exit(0);
             }
         }
         if(!exclude_left) S.push(std::make_tuple(LD, M2, FD, FM, V2 - FU, V3));
         if(!exclude_right) S.push(std::make_tuple(M1, RU, V0 - FD, V1, FU, -FM));
     }
     else{
         my_complex M1{RU.real(), (LD.imag() + RU.imag()) / 2};
         my_complex M2{LD.real(), (LD.imag() + RU.imag()) / 2};
         my_complex X1{RU.real(), LD.imag()};
         my_complex X2{LD.real(), RU.imag()};
         my_real FM;
         while(not var_arg(M1, M2, FM)){
             M1 += dy;
             M2 += dy;
         }
         my_real r = abs(LD - RU) / 2;
         my_real FD = 0, FU = 0;
         bool exclude_down = USE_EXCLUSION_TEST ? exclusion_test(LD, M1) : false;
         if(!exclude_down) var_arg(X1, M1, FD);
         bool exclude_up = USE_EXCLUSION_TEST ? exclusion_test(M2, RU) : false;
         if(!exclude_up) var_arg(X2, M2, FU);
         if(!USE_EXCLUSION_TEST or (!exclude_down and !exclude_up)){
             my_real v1 = V0 + FD + FM + V3 - FU;
             my_real v2 = V1 - FD + V2 + FU - FM;
             if(abs(v1 / 2 / PI) < 0.1) exclude_down = true;
             if(abs(v2 / 2 / PI) < 0.1) exclude_up = true;
             if(abs(round((v1 + v2) / 2 / PI) - zeros) > 0.1){
                 std::cerr << "Fatal error. Sum of partitions isn't equal to total" << std::endl;
                 exit(0);
             }
         }
         if(!exclude_down) S.push(std::make_tuple(LD, M1, V0, FD, FM, V3 - FU));
         if(!exclude_up) S.push(std::make_tuple(M2, RU, -FM, V1 - FD, V2, FU));
     }
 }
 
 void solve_by_thread(int id){
     while(Answers.size() < total_zeros){
         std::tuple<my_complex, my_complex, my_real, my_real, my_real, my_real> cur = S.pop();
         if((std::get<0>(cur)).real() > my_real(2)) continue;
         process_one(cur, id);
     }
     active_threads -= 1;
 }
 
 void init(my_complex LD, my_complex RU, int n){
     if(USE_EXCLUSION_TEST and exclusion_test(LD, RU)){
         total_zeros = 0;
         return;
     }
     my_complex LU{LD.real(), RU.imag()};
     my_complex RD{RU.real(), LD.imag()};
     my_real P1, P2, P3, P4;
     var_arg(LD, RD, P1);
     var_arg(RD, RU, P2);
     var_arg(RU, LU, P3);
     var_arg(LU, LD, P4);
     total_zeros = (int)round((P1 + P2 + P3 + P4) / 2 / PI);
     S.push(std::make_tuple(LD, RU, P1, P2, P3, P4));
     std::cerr << "Initialized" << std::endl;
     std::cerr << "\033[2J\033[1;1H";
     std::cerr << "Now we've got " << cnt << "/" << total_zeros << " zeros" << '\n';
     std::cerr << "Number of active threads: " << active_threads << '\n';
     std::cerr << "Excluded rectangles: " << excluded << '\n';
     auto cur = std::chrono::system_clock::now();
     std::time_t cur_time = std::chrono::system_clock::to_time_t(cur);
     print_progress(cnt, total_zeros);
     std::cerr << "Registered at " << std::ctime(&cur_time) << '\n';
 
 }
 
 void solve(int argc, char* argv[]){
 ///	Custom rectangle for computing zeros
     my_complex LD{atof(argv[2]), atof(argv[3])}, RU{atof(argv[4]), atof(argv[5])};
 
     active_threads = limit;
     init(LD, RU, n);
     std::vector<std::thread> threads(limit);
     for(int i = 0; i < limit; i++){
         threads[i] = std::thread(solve_by_thread, i + 1);
     }
     for(int i = 0; i < limit; i++){
         threads[i].join();
     }
 }
 
 void download(){
     int cnt = 0;
     std::ofstream out; // zeros.txt
     std::ofstream outx; // zeros.gp
     out.open(outputname, std::ios_base::app);
     outx.open(outputnamex, std::ios_base::app);
     out.precision(std::numeric_limits<my_complex::value_type>::digits10);
     outx.precision(std::numeric_limits<my_complex::value_type>::digits10);
     outx << "zeroes = [";
     while(Answers.size()){
         my_complex s;
         my_real R;
         std::tie(s, R) = Answers.pop();
         out << s.real() << " " << s.imag() << " " << R << '\n';
         if(cnt) outx << ", ";
         outx << "[" << s.real() << " + " << s.imag() << " * I, " << R << "]";
         cnt = 1;
     }
     outx << "];\n";
     out.close();
     outx.close();
 }
 
 template<typename duration_type>
 void print_final_information(duration_type duration) {
     std::ofstream logger;
     logger.open(benchmarkname);
     std::cout << "Found " << Answers.size() << "/" << total_zeros << " zeros" << '\n';
     std::cout << "Excluded " << excluded << " rectangles during process" << '\n';
     print_progress(Answers.size(), total_zeros);
     std::cerr << "Process executed in " << duration.count() / 1000000000.0 << " s" << '\n';
     logger << "Process executed in " << duration.count() / 1000000000.0 << " s" << '\n';
     logger.close();
 
 }
 
 int main(int argc, char* argv[]){
     if(argc != 6) {
         std::cout << "Usage:" << std::endl;
         std::cout << "   ./a.out (n) (lower-left real) (lower-left imag) (upper-right real) (upper-right imag)" << std::endl;
         return 0;
     }
     funcinit(argc, argv);
     auto start = std::chrono::high_resolution_clock::now();
     solve(argc, argv);
     auto stop = std::chrono::high_resolution_clock::now();
     auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
     print_final_information(duration);
     download();
 //	test_newton();
     return 0;
 }
 