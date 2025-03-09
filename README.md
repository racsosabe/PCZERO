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
- `base_function` typename: The function type that will be our f and df type. This can be:
    - A normal function that receives a `my_complex x` as an argument and returns a `my_complex`. For example, `my_complex f(my_complex x) { return x * x; }`, or
    - A class that overloads the operator `()` to receive a `my_complex` and return a `my_complex`. For example, `my_complex operator () (my_complex x) { return x * x; }`. This class must have the `=` operator overloaded as well.
- `bound_function` typename: The function type that will be our M and N type. This can be either:
    - A normal function that receives two `my_complex s1, my_complex s2` as arguments and returns a `my_real`. For example, `my_complex M(my_complex s1, my_complex s2) { return abs(s1 * s2); }`, or
    - A class that overloads the operator `()` to receive two `my_complex` varaibles and return a `my_real`. For example, `my_real operator () (my_complex s1, my_complex s2) { return abs(s1 * s2); }`. This class must have the `=` operator overloaded as well.
- `USE_EXCLUSION_TEST` flag: Set to 1 if the exclusion test will be used on rectangles.
- `thread_limit` value: Number of threads that will be used during the execution.

And the usage of the constructor is:

`pczero<my_complex, my_real, thread_limit, USE_EXCLUSION_TEST> Solver(f, df, M, N);`

Where $f$, $df$, $M$ and $N$ are the implemented functions/classes.

**Note**: Please consider that the given rectangle in the complex plane might be modified if there is any zero on its border, the final rectangle will be printed once the computation starts.

This program writes to the following files. If they already exist, the new information is appended to the end:

-   `zeros.txt`: Contains the computed zeros s in the format `Real Imaginary |f(s)|`.
-   `log.txt`: Contains the execution time for processing the zeros.

## Examples

You will find three usages of the `PCZERO.hpp` header in the `examples` folder. It contains:

 - Dirichlet polynomial example with class (Partial sum of Riemann zeta function)
 - Exponential polynomial example with class (no exclusion test)
 - Exponential polynomial example with class (with exclusion test)
 - Dirichlet polynomial example without class (Partial sum of Riemann zeta function)

The Dirichlet polynomial examples do not use exclusion test.

### Dirichlet polynomial example with class (Partial sum of Riemann zeta function)

It is shown how to include the header of the Dirichlet polynomial class and usage of PCZERO.

You can compile it (assuming that the headers are correctly imported) with the following command:

```bash
g++ main_zeta.cpp -pthread -o main_zeta
```

**If you are using Ubuntu 24.04, please add the flag `-march=native` to the compilation command.**

Then, you can execute it with the following command:

```bash
./main_zeta (lower-left real part) (lower-left imaginary part) (upper-right real part) (upper-right imaginary part)
```

The given example uses the function:

$$ \zeta_{5}(s)=1+2^{-s}+3^{-s}+4^{-s}+5^{-s} $$

However, the value of the iteration limit $n$ can be modified while declaring the class instances.

### Exponential polynomial example with class (no exclusion test)

It is shown how to include the header of the Exponential polynomial class and usage of PCZERO.

You can compile it (assuming that the headers are correctly imported) with the following command:

```bash
g++ main_exponential_no_exclusion.cpp -pthread -o main_exponential_no_exclusion
```

**If you are using Ubuntu 24.04, please add the flag `-march=native` to the compilation command.**

Then, you can execute it with the following command:

```bash
./main_exponential_no_exclusion (lower-left real part) (lower-left imaginary part) (upper-right real part) (upper-right imaginary part)
```

The given example uses the function:

$$ f(s)=-1+(-6-5s)e^{-2s} + (-2 + 3s)e^{-2\pi/3s} + 2e^{-5/2s} + se^{-3s} + 3e^{-4s} $$

### Exponential polynomial example with class (with exclusion test)

It is shown how to include the header of the Exponential polynomial class and usage of PCZERO.

You can compile it (assuming that the headers are correctly imported) with the following command:

```bash
g++ main_exponential_with_exclusion.cpp -pthread -o main_exponential_with_exclusion
```

**If you are using Ubuntu 24.04, please add the flag `-march=native` to the compilation command.**

Then, you can execute it with the following command:

```bash
./main_exponential_with_exclusion (lower-left real part) (lower-left imaginary part) (upper-right real part) (upper-right imaginary part)
```

The given example uses the function:

$$ f(s)=-1+(-6-5s)e^{-2s} + (-2 + 3s)e^{-2\pi/3s} + 2e^{-5/2s} + se^{-3s} + 3e^{-4s} $$

Notice that this implementation differs with the previous one just in the declaration of the pczero class `pczero<my_complex, my_real, 8, 1>` instead of `pczero<my_complex, my_real, 8, 0>`.

### Dirichlet polynomial example without class (Partial sum of Riemann zeta function)

It is shown how to include the header of the f, df, M and N functions and usage of PCZERO.

You can compile it (assuming that the headers are correctly imported) with the following command:

```bash
g++ main_zeta_no_class.cpp -pthread -o main_zeta_no_class
```

**If you are using Ubuntu 24.04, please add the flag `-march=native` to the compilation command.**

Then, you can execute it with the following command:

```bash
./main_zeta_no_class (lower-left real part) (lower-left imaginary part) (upper-right real part) (upper-right imaginary part) (n)
```

And it will compute the zeros of the function:

$$ \zeta_{n}(s)=1+2^{-s}+3^{-s}+\cdots +n^{-s} $$

An example of the execution is:

```bash
./main_zeta_no_class -4 0 1.74 100 5
```

This example uses the function:

$$ \zeta_{5}(s)=1+2^{-s}+3^{-s}+4^{-s} +5^{-s} $$

Which will produce the 26 zeros inside the given rectangle over $\zeta_{5}$ in `zeros.txt`. The corresponding plot is the following:

![](https://i.imgur.com/khfrUmI.jpeg)

**Note**: For all the examples that use classes, the base function is defined in their own headers.

## Authors

-   Oswaldo Velasquez
-   Manuel Toribio
-   Victor Galvan