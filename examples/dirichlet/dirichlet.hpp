#include<array>
#include<cassert>
#include<fstream>
#include<sstream>
#include<boost/multiprecision/cpp_bin_float.hpp>

using my_complex = boost::multiprecision::cpp_complex<334, boost::multiprecision::backends::digit_base_2, void, std::int32_t, -262142, 262143>;
using my_real = boost::multiprecision::number<boost::multiprecision::backends::cpp_bin_float<334, boost::multiprecision::backends::digit_base_2, void, std::int32_t, -262142, 262143>, boost::multiprecision::et_off>;

// A Dirichlet polynomial f(s) = \sum_{k=0}^n a_k e^{-l_k s}, where 
// 0 = l_0< l_1 < l_2 < ... < l_n is represented by the following structure
template<typename my_real, typename my_complex, int MAX>
class dirichlet {
    int len;                 // n
    std::array<my_real, MAX> freqs;      // l_k is freqs[k]
    std::array<my_real, MAX> coeff;      // a_k is coeffs[k]
public:

    int getLength() const {
        return (*this).len;
    }

    my_real getFreq(int pos) const {
        assert(0 <= pos and pos < MAX);
        return (*this).freqs[pos];
    }

    my_real getCoef(int pos) const {
        assert(0 <= pos and pos < MAX);
        return (*this).coeff[pos];
    }

    void setLength(int len) {
        (*this).len = len;
    }

    void setFreq(int pos, my_real value) {
        assert(0 <= pos and pos < MAX);
        (*this).freqs[pos] = value;
    }

    void setCoef(int pos, my_real value) {
        assert(0 <= pos and pos < MAX);
        (*this).coeff[pos] = value;
    }
    
    my_complex operator () (const my_complex &x) const {
        my_complex ans = 0;
        int le = getLength();
        for(int i = 0; i < le; i++) {
            ans += getCoef(i) * exp(-getFreq(i) * x);
        }
        return ans;
    }
    
    my_real operator () (const my_complex &s1, const my_complex &s2) const {
    	my_real min_real = std::min(s1.real(), s2.real());
        return abs((*this)(min_real));
    }

    // Computing the derivative of the above function
    dirichlet<my_real, my_complex, MAX> derivative() {
        dirichlet<my_real, my_complex, MAX> tmp;
        int le = getLength();
        tmp.setLength(le);
        for(int k = 0; k < le; k++) {
            tmp.setFreq(k, getFreq(k));
            tmp.setCoef(k, -getCoef(k) * getFreq(k));
        }
        return tmp;
    }

    // Computing the M bound for the second derivative of the function f
    dirichlet<my_real, my_complex, MAX> Mbound() {
        dirichlet<my_real, my_complex, MAX> tmp = derivative().derivative();
        int le = tmp.getLength();
        for(int i = 0; i < le; i++) {
            my_real coef = tmp.getCoef(i);
            tmp.setCoef(i, abs(coef));
        }
        return tmp;
    }
    
    void operator = (const dirichlet<my_real, my_complex, MAX> &f) {
    	int le = f.getLength();
    	setLength(le);
    	for(int i = 0; i < le; i++) {
    	    setFreq(i, f.getFreq(i));
    	    setCoef(i, f.getCoef(i));
    	}
    }
    
    std::string toString() const {
    	std::ostringstream oss;
    	int le = getLength();
    	oss << "n = " << le << std::endl;
    	oss << "Freqs: ";
    	for(int i = 0; i < le; i++) {
    	    oss << getFreq(i) << " \n"[i + 1 == le];
    	}
    	oss << "Coeffs: ";
    	for(int i = 0; i < le; i++) {
    	    oss << getCoef(i) << " \n"[i + 1 == le];
    	}
    	return oss.str();
    }
};

// Definition of the matrix associated to the Riemann Zeta Function partial sum zeta_n = \sum_{k=0}^n k^{-s}
// It can be customized for other Dirichlet polynomials, and with slight modifications, for exponential polynomials (with real frequencies)
template<typename my_real, typename my_complex, int MAX>
dirichlet<my_real, my_complex, MAX> zetan_matrix(int v){
    dirichlet<my_real, my_complex, MAX> tmp;
    tmp.setLength(v);
    tmp.setFreq(0, my_real(0));
    tmp.setCoef(0, my_real(1));
    for(int k = 1; k < v; k++) {
        tmp.setFreq(k, log(k + 1));
        tmp.setCoef(k, 1);
    }
    return tmp;
}

dirichlet<my_real, my_complex, 5> base_function() {
    return zetan_matrix<my_real, my_complex, 5>(5);
}
