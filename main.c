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
    // Main loop
    while(( line = fgets(buffer, sizeof(buffer), stdin)) != NULL) {
        mpz_t current_N;
        mpz_init(current_N);
        mpz_set_str(current_N, line, 10);
        factorize(current_N, factors);
    }
    return 0;
}

/*
 * Factorizes N using Fermat's method
 */
void factorize(mpz_t N, mpz_t factors[]) {
    mpz_t sqrt_N;
    mpz_init(sqrt_N);

    int num_factors;
    num_factors = find_trivial_factors(N, factors);

    if(mpz_cmp_si(N,1)) {
        // Primality testing
        if(mpz_probab_prime_p(N, 5)) {
            gmp_printf("%Zd\n", N);
            mpz_set_ui(N, 1);

        }else if(mpz_perfect_square_p(N) != 0) {
                // Check for perfect squares once
                mpz_root(sqrt_N, N, 2);
                gmp_printf("%Zd\n", sqrt_N);
                gmp_printf("%Zd\n", sqrt_N);
                mpz_set_ui(N, 1);
        }else if(mpz_sizeinbase(N,2) < 80){

            int result;
            while (mpz_cmp_si(N, 1) != 0) {
                if(mpz_probab_prime_p(N, 5)) {
                    gmp_printf("%Zd\n", N);
                    mpz_set_ui(N, 1);
                    break;
                }
                result = fermat(N, factors, num_factors);
                if (!result) {
                    break;
                }
                ++num_factors;
            }
        }
    }

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
    mpz_t factor;
    mpz_init(factor);
    for (i = 0; i < FIRST_PRIMES_SIZE; ++i) {
        remainder = mpz_fdiv_ui(N, first_primes[i]);
        if (remainder == 0) {
            // Remove all occurrences of this factor from N.
            // All occurrences must be stored.
            while (remainder == 0) {
                mpz_set_ui(factor, first_primes[i]);
                mpz_set(factors[factor_index], factor);
                ++factor_index;
                mpz_fdiv_q_ui(N, N, first_primes[i]);
                remainder = mpz_fdiv_ui(N, first_primes[i]);
            }
        } 
    }
    mpz_clear(factor);
    return factor_index;
}


int fermat(mpz_t N, mpz_t factors[], int num_factors) {
    
    mpz_t a;
    mpz_init(a);
    // Ugly fix for a = ceil(sqrt(N)), does a = trunc_sqrt(N-1) + 1 instead
    mpz_sub_ui(a, N, 1);
    mpz_sqrt(a, a);
    mpz_add_ui(a, a, 1);

    // Set b2 to a^2 - N
    mpz_t b2;
    mpz_init(b2);
    mpz_mul(b2, a, a);
    mpz_sub(b2, b2, N);

    mpz_t count;
    mpz_init_set_ui(count, 0);
    mpz_t limit;
    mpz_init(limit);
    
    // Set the iteration limit
    mpz_set_ui(limit,3000000);

    while(mpz_cmp(count, limit)<0) {

        if(mpz_perfect_square_p(b2) != 0) {
            // Get the factor from d = a - sqrt(b2)
            mpz_sqrt(b2, b2);
            mpz_sub(b2, a, b2);
            mpz_set(factors[num_factors],b2);
            mpz_fdiv_q(N, N, b2);

            // Clear variables
            mpz_clear(a);
            mpz_clear(b2);
            mpz_clear(count);
            mpz_clear(limit);
            return 1; 
        }

        // a = a + 1
        mpz_add_ui(a, a, 1);
        // b2 = a^2 - N
        mpz_mul(b2, a, a);
        mpz_sub(b2, b2, N);

        mpz_add_ui(count,count,1);
    }

    // Clear variables
    mpz_clear(a);
    mpz_clear(b2);
    mpz_clear(count);
    mpz_clear(limit);
    return 0;
}

void print_factors(mpz_t factors[], int num_factors) {
    int i;
    for (i = 0; i < num_factors; ++i) {
        // This if-check is not needed anymore.
        if (factors[i] == 0) {
            break;
        }
        gmp_printf("%Zd\n", factors[i]);
    }
}
