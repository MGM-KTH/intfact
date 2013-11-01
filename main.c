#include "main.h"

#define TRUE 1
#define FALSE 0
#define MAX_LINE_LENGTH 256
#define FAIL_STRING "fail"

int main(int argc, char *argv[]) {
    int i;
    static char buffer[MAX_LINE_LENGTH];
    const char *line;
    mpz_t factors[FACTORS_ARRAY_SIZE];
    for (i = 0; i < FACTORS_ARRAY_SIZE; ++i) {
        mpz_init(factors[i]);
    }
    while(( line = fgets(buffer, sizeof(buffer), stdin)) != NULL) {
        mpz_t current_N;
        mpz_init(current_N);
        mpz_set_str(current_N, line, 10);
        factorize(current_N, factors);
    }
    return 0;
}

/*
 * Factorizes N using Pollard's rho method
 */
void factorize(mpz_t N, mpz_t factors[]) {
    // gmp_printf("Factorizing %Zd\n", N);
    int num_factors;

    num_factors = find_trivial_factors(N, factors);
    
    // If N != 1 (not fully factorized)
    if(mpz_cmp_si(N, 1) != 0) {    
        num_factors = pollards(N, factors, num_factors);
    }

    // gmp_printf("New N after trivial pruning: %Zd\n", N);

    // If N != 1 (not fully factorized)
    if (mpz_cmp_si(N, 1)) {
        printf("%s\n", FAIL_STRING); // Print this if you can't factorize
    }
    else {
        print_factors(factors, num_factors);
    }
    printf("\n");
}

int find_trivial_factors(mpz_t N, mpz_t factors[]) {
    int i;
    int factor_index = 0;
    unsigned long int remainder;
    for (i = 0; i < FIRST_PRIMES_SIZE; ++i) {
        remainder = mpz_fdiv_ui(N, first_primes[i]);
        if (remainder == 0) {
            // Remove all occurrences of this factor from N.
            // All occurrences must be stored.
            while (remainder == 0) {
                mpz_t factor;
                mpz_init(factor);
                mpz_set_ui(factor, first_primes[i]);
                mpz_set(factors[factor_index], factor);
                ++factor_index;
                mpz_fdiv_q_ui(N, N, first_primes[i]);
                remainder = mpz_fdiv_ui(N, first_primes[i]);
            }
        } 
    }
    return factor_index;
}


int pollards(mpz_t N, mpz_t factors[], int num_factors) {
    
    // Fix random number container
    mpz_t rand;
    mpz_init(rand);

    // Initialize and seed a randstate
    gmp_randstate_t rand_state;
    gmp_randinit_default(rand_state);
    gmp_randseed_ui(rand_state, time(NULL));

    // Get random number
    mpz_urandomm(rand, rand_state, N);

    //gmp_printf("number N: %Zd random number %Zd\n", N, rand);

    //next_in_seq(next, prev, n);
    //gmp_printf("next in sequence: %Zd\n", next);

    // TDOO: decide how to loop.
    //fix calculate_sequence function
    // Do GCD in gmp.
    // put into algorithm!
    return num_factors;
}

void next_in_seq(mpz_t next, mpz_t prev, mpz_t N) {
    mpz_pow_ui(next, prev, 2); // X^2
    mpz_add_ui(next, next, 1); // X^2 + 1
    mpz_mod(next, next, N);    // (X^2 + 1) mod N
}

void print_factors(mpz_t factors[], int num_factors) {
    int i;
    for (i = 0; i < num_factors; ++i) {
        if (factors[i] == 0) {
            break;
        }
        gmp_printf("%Zd\n", factors[i]);
    }
}
