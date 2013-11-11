#ifndef MAIN_H
#define MAIN_H

#define FACTORS_ARRAY_SIZE 100

#include <stdio.h>
#include <gmp.h>
#include <time.h>
#include "constants.h"

int main(int argc, char *argv[]);

void factorize(mpz_t N, mpz_t factors[]);
void print_factors(mpz_t factors[], int num_factors);

int find_trivial_factors(mpz_t N, mpz_t factors[]);
int fermat(mpz_t N, mpz_t factors[], int num_factors);


#endif // end header guard
