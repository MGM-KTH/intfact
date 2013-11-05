#ifndef MAIN_H
#define MAIN_H

#define FACTORS_ARRAY_SIZE 100

#include <stdio.h>
#include <gmp.h>
#include <time.h>
#include "constants.h"
#include <math.h>

int main(int argc, char *argv[]);

void factorize(mpz_t N, mpz_t factors[]);
void print_factors(mpz_t factors[], int num_factors);
void next_in_seq(mpz_t next, mpz_t prev, mpz_t N);

int find_trivial_factors(mpz_t N, mpz_t factors[]);
int pollards(mpz_t N, mpz_t factors[], int num_factors);
long int min(long int a, long int b);


#endif // end header guard
