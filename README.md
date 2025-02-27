# PCZERO

## Description

This code is an implementation of the algorithm described in the paper "A Reliable and Parallel Algorithm for the Computation of Zeros of Analytic Functions", called PCZERO.

## Requirements

This code requires the following:

-   **Boost Multiprecision Library**: Version 1.74 for arbitrary precision (Optional if high precision is not needed).
-   **g++ Compiler**: Must support the `-std=c++17` flag
-   **64-bit Operating System**: The current version has been tested on Linux/PC.

## Compilation and Usage

To use the `PCZERO.hpp`, just include it in your main. The PCZERO class requires the following information:

- `my_complex` class: The class that will be our complex number representation.
- `my_real` class: The class that will be our real number representation
- `polynomial` class: The class that will be our function type. This class must have implemented:
    - `my_complex evaluate(my_complex x)`: Must evaluate the function on the given complex.
    - `std::string toString()`: Returns a string representation of the function (May be left as returning an empty string).
    - `polynomial derivative()`: Must return the derivative of the function instance.
    - `polynomial Mbound()`: Must return the M function (as explained in the article).
    - `my_real evaluate_as_m(my_complex s1, my_complex s2)`: Evaluate the M value of the given complex numbers.
    - `=` operator overloaded. 
- `USE_EXCLUSION_TEST` flag: Set to 1 if the exclusion test will be used on rectangles.
- `thread_limit` value: Number of threads that will be used during the execution.

**Note**: Please consider that the given rectangle in the complex plane might be modified if there is any zero on its border, the final rectangle will be printed once the computation starts.

This program writes to the following files. If they already exist, the new information is appended to the end:

-   `zeros.txt`: Contains the computed zeros s in the format `Real Imaginary |f(s)|`.
-   `zeros.gp`: Contains the computed zeros as an array of tuples `[s, |f(s)|]` in PARI/GP format.
-   `log.txt`: Contains the execution time for processing the zeros.

## Examples

You will find two usages of the `PCZERO.hpp` header in the `examples` folder. It contains:

- Dirichlet polynomial example (Partial sum of Riemann zeta function): It is shown how to include the header of the Dirichlet polynomial class and usage of PCZERO.

You can compile it (assuming that the headers are correctly imported) with the following command:

```bash
g++ main_zeta.cpp -pthread -o main_zeta
```

Then, you can execute it with the following command:

```bash
./main_zeta (lower-left real part) (lower-left imaginary part) (upper-right real part) (upper-right imaginary part)
```

- Exponential polynomial example: It is shown how to include the header of the Exponential polynomial class and usage of PCZERO.

You can compile it (assuming that the headers are correctly imported) with the following command:

```bash
g++ main_exponential.cpp -pthread -o main_exponential
```

Then, you can execute it with the following command:

```bash
./main_exponential (lower-left real part) (lower-left imaginary part) (upper-right real part) (upper-right imaginary part)
```

## Authors

-   Oswaldo Velasquez
-   Manuel Toribio
-   Victor Galvan