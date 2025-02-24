# PCZERO

## Description

This code is an implementation of the algorithm described in the paper "A Reliable and Parallel Algorithm for the Computation of Zeros of Analytic Functions", called PCZERO.

## Requirements

This code requires the following:

-   **Boost Multiprecision Library**: Version 1.74
-   **g++ Compiler**: Must support the `-std=c++17` flag
-   **64-bit Operating System**: The current version has been tested on Linux/PC.

## Compilation and Usage

In the code, only the following variables need to be modified:

-   `USE_EXCLUSION_TEST`: Boolean to use exclusion test for rectangles or not (Line 157).

After making these changes, compile and execute the code as follows:

```bash
g++ Parallel-final.cpp -o Parallel-final -pthread -O3
./Parallel-final (n) (lower-left real part) (lower-left imaginary part) (upper-right real part) (upper-right imaginary part)
```

This program writes to the following files. If they already exist, the new information is appended to the end:

-   `zeros{n}.txt`: Contains the computed zeros ss in the format `Real Imaginary |f(s)|`.
-   `zeros{n}.gp`: Contains the computed zeros as an array of tuples `[s, |f(s)|]` in PARI/GP format.
-   `log{n}.txt`: Contains the execution time for processing the zeros.

## Authors

-   Oswaldo Velasquez
-   Manuel Toribio
-   Victor Galvan