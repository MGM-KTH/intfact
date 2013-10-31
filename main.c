#include "main.h"

#define TRUE 1
#define FALSE 0
#define MAX_LINE_LENGTH 256
#define FAIL_STRING "fail"

int main(int argc, char *argv[]) {
    static char buffer[MAX_LINE_LENGTH];
    const char *line;
    while(( line = fgets(buffer, sizeof(buffer), stdin)) != NULL) {
        mpz_t current_N;
        mpz_init(current_N);
        mpz_set_str(current_N, line, 10);
        unsigned long int factors[FACTORS_ARRAY_SIZE];
        factorize(current_N, factors);
    }
    return 0;
}

/*
 * Factorizes N using Pollard's rho method
 */
void factorize(mpz_t N, unsigned long int factors[FACTORS_ARRAY_SIZE]) {
    // gmp_printf("Factorizing %Zd\n", N);

    int num_factors = find_trivial_primes(N, factors);
    // printf("Number of factors found: %d\n", num_factors);
    // gmp_printf("New N after trivial pruning: %Zd\n", N);

    if (mpz_cmp_si(N, 1)) {
        printf("%s\n", FAIL_STRING); // Print this if you can't factorize
    }
    else {
        print_factors(factors, num_factors);
    }
    printf("\n");
}

int find_trivial_primes(mpz_t N, unsigned long int factors[FACTORS_ARRAY_SIZE]) {
    int i;
    int factor_index = 0;
    unsigned long int remainder;
    for (i = 0; i < FIRST_PRIMES_SIZE; ++i) {
        remainder = mpz_fdiv_ui(N, first_primes[i]);
        if (remainder == 0) {
            while (remainder == 0) {
                factors[factor_index] = first_primes[i];
                ++factor_index;
                mpz_fdiv_q_ui(N, N, first_primes[i]);
                remainder = mpz_fdiv_ui(N, first_primes[i]);
            }
        } 
    }
    return factor_index;
}



void print_factors(unsigned long int factors[FACTORS_ARRAY_SIZE], int num_factors) {
    int i;
    for (i = 0; i < num_factors; ++i) {
        if (factors[i] == 0) {
            break;
        }
        printf("%lu\n", factors[i]); 
    }
}