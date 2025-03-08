#include <iostream>
#include <vector>
#include <sstream>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
#include <boost/multiprecision/complex_adaptor.hpp>
#include <cmath>

using my_complex = boost::multiprecision::cpp_complex<334, boost::multiprecision::backends::digit_base_2, void, std::int32_t, -262142, 262143>;
using my_real = boost::multiprecision::number<boost::multiprecision::backends::cpp_bin_float<334, boost::multiprecision::backends::digit_base_2, void, std::int32_t, -262142, 262143>, boost::multiprecision::et_off>;

template<typename my_real, typename my_complex>
class ExponentialPolynomial {
    struct Term {
        my_real lambda;                          // Frecuencia lambda_k
        std::vector<my_complex> coefficients;    // Coeficientes a_{k,j}

        Term(const my_real& lambda_, const std::vector<my_complex>& coeffs)
            : lambda(lambda_), coefficients(coeffs) {}

        // Evalúa p_k(s) e^{-lambda_k s}
        my_complex evaluate(const my_complex& s) const {
            my_complex poly_value = 0.0;
            my_complex s_power = 1.0;
            for (const auto& coeff : coefficients) {
                poly_value += coeff * s_power;
                s_power *= s;
            }
            return poly_value * boost::multiprecision::exp(-lambda * s);
        }

        // Devuelve una representación en string del término
        std::string toString() const {
            std::ostringstream oss;
            oss << "(";
            for (size_t j = 0; j < coefficients.size(); ++j) {
                if (j > 0) oss << " + ";
                oss << "(" << coefficients[j] << ")";
                if (j > 0) oss << "*s^" << j;
            }
            oss << ") * exp(-" << lambda << " * s)";
            return oss.str();
        }

        // Calcula la derivada del término
        Term derivative() const {
            std::vector<my_complex> deriv_coeffs(coefficients.size(), my_complex(0, 0));

            // Derivada del polinomio p_k(s)
            for (size_t j = 1; j < coefficients.size(); ++j) {
                deriv_coeffs[j - 1] += coefficients[j] * static_cast<double>(j);
            }

            // Resta lambda_k * p_k(s)
            for (size_t j = 0; j < coefficients.size(); ++j) {
                deriv_coeffs[j] -= lambda * coefficients[j];
            }

            return Term(lambda, deriv_coeffs);
        }

        // Evalúa e^{-lambda_k * min(Re(s1), Re(s2))} * |p_k|(max(|s1|, |s2|))
        my_real evaluate_as_m(const my_complex& s1, const my_complex& s2) const {
            my_real min_real = std::min(s1.real(), s2.real());
            my_real max_abs = std::max(abs(s1), abs(s2));

            // Calcula |p_k|(max(|s1|, |s2|))
            my_real poly_abs = 0.0;
            my_real s_power = 1.0;
            for (const auto& coeff : coefficients) {
                poly_abs += abs(coeff) * s_power;
                s_power *= max_abs;
            }

            return boost::multiprecision::exp(-lambda * min_real) * poly_abs;
        }
    };

    std::vector<Term> terms;

public:

    // Agrega un término p_k(s) e^{-lambda_k s}
    void addTerm(const my_real& lambda, const std::vector<my_complex>& coeffs) {
        terms.emplace_back(lambda, coeffs);
    }

    std::vector<Term> getTerms() const {
        return terms;
    }

    // Evalúa P(s) en s ∈ ℂ
    my_complex operator () (const my_complex& s) const {
        my_complex result = 0.0;
        for (const auto& term : terms) {
            result += term.evaluate(s);
        }
        return result;
    }

    // Devuelve una representación en string de P(s)
    std::string toString() const {
        std::ostringstream oss;
        oss << "P(s) = ";
        for (size_t k = 0; k < terms.size(); ++k) {
            if (k > 0) oss << " + ";
            oss << terms[k].toString();
        }
        return oss.str();
    }

    // Calcula la derivada del ExponentialPolynomial
    ExponentialPolynomial<my_real, my_complex> derivative() const {
        ExponentialPolynomial<my_real, my_complex> deriv_poly;
        for (const auto& term : terms) {
            deriv_poly.terms.push_back(term.derivative());
        }
        return deriv_poly;
    }

    ExponentialPolynomial<my_real, my_complex> Mbound() const {
        ExponentialPolynomial<my_real, my_complex> tmp = (*this).derivative().derivative();
        std::vector<Term> current_terms = tmp.getTerms();
        ExponentialPolynomial<my_real, my_complex> result;
        for(auto term : current_terms) {
            for(auto &values : term.coefficients) {
                values = abs(values);
            }
            result.addTerm(term.lambda, term.coefficients);
        }
        return result;
    }

    // Evalúa la suma personalizada para s1 y s2
    my_real operator () (const my_complex& s1, const my_complex& s2) const {
        my_real result = 0.0;
        for (const auto& term : terms) {
            result += term.evaluate_as_m(s1, s2);
        }
        return result;
    }

    void operator = (const ExponentialPolynomial<my_real, my_complex> &f) {
        (*this).terms = f.getTerms();
    }
};

ExponentialPolynomial<my_real, my_complex> base_function() {
    ExponentialPolynomial<my_real, my_complex> poly;
    // Implementación de f(s) de la imagen
    poly.addTerm(my_real(0.0), {my_complex(-1, 0)});                                       // -1
    poly.addTerm(my_real(2.0), {my_complex(-6, 0), my_complex(-5, 0)});                    // (-6 - 5s)e^{-2s}
    poly.addTerm(my_real(2.0 * acos(my_real(-1)) / 3.0), {my_complex(-2, 0), my_complex(3, 0)}); // (-2 + 3s)e^{-2pi/3 s}
    poly.addTerm(my_real(5.0 / 2.0), {my_complex(2, 0)});                                  // 2e^{-5/2 s}
    poly.addTerm(my_real(3.0), {my_complex(0, 0), my_complex(1, 0)});                      // se^{-3s}
    poly.addTerm(my_real(4.0), {my_complex(3, 0)});                                        // 3e^{-4s}
    return poly;
}
