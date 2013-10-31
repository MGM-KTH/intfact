#ifndef MAIN_H
#define MAIN_H

#define FACTORS_ARRAY_SIZE 100

#include <stdio.h>
#include <gmp.h>
#include "constants.h"

int main(int argc, char *argv[]);
void factorize(mpz_t N, unsigned long int factors[FACTORS_ARRAY_SIZE]);
int find_trivial_primes(mpz_t N, unsigned long int factors[FACTORS_ARRAY_SIZE]);
void print_factors(unsigned long int factors[FACTORS_ARRAY_SIZE], int num_factors);


#endif // end header guard