#include "main.h"

#define TRUE 1
#define FALSE 0
#define MAX_LINE_LENGTH 256
#define FAIL_STRING "fail"

int C = 1;

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
 * Factorizes N using Pollard's rho method (Including Brent's improvement)
 */
void factorize(mpz_t N, mpz_t factors[]) {
    mpz_t sqrt_N;
    mpz_init(sqrt_N);
    int num_factors;
    num_factors = find_trivial_factors(N, factors);

    if (mpz_cmp_si(N, 1)) {
        // Primality testing
        if(mpz_probab_prime_p(N, 5)) {
            gmp_printf("%Zd\n", N);
            mpz_set_ui(N, 1);
        }   
        // Check for perfect squares once
        else if (mpz_perfect_square_p(N) != 0) {
            mpz_root(sqrt_N, N, 2);
            gmp_printf("%Zd\n", sqrt_N);
            gmp_printf("%Zd\n", sqrt_N);
            mpz_set_ui(N, 1);
        }
        // Bit-limit
        else if (mpz_sizeinbase(N,2) < 99) {
            int result;
            while (mpz_cmp_si(N, 1) != 0) {
                // Pollard's
                result = pollards(N, factors, num_factors);
                if (!result) {
                    break;
                }
                ++num_factors;
                // Primality testing
                if(mpz_probab_prime_p(N, 5)) {
                    gmp_printf("%Zd\n", N);
                    mpz_set_ui(N, 1);
                    break;
                }
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

int pollards(mpz_t N, mpz_t factors[], int num_factors) {
    mpz_t x;
    mpz_t y;
    mpz_t d;
    mpz_t q;
    mpz_t temp_var;
    mpz_t ys;
    mpz_init_set_ui(x, 2);
    mpz_init_set_ui(y, 2);
    mpz_init_set_ui(q, 1);
    mpz_init(temp_var);
    mpz_init(ys);
    mpz_init(d);

    long int r = 1;
    long int m = 100; //((size_in_base_two) < 67) ? 67 : (size_in_base_two); 
    long int k = 0;
    long int count = 0;

    long int limit; // = 130000; // Fixed limit
    int size_in_base_two = mpz_sizeinbase(N,2);
    // Custom iteration limits
    if (size_in_base_two > 92) {
        limit = 950;
    }
    else if (size_in_base_two > 88) {
        limit = 500;
    }
    else {
        limit = 90000;
    }

    // Pollard's
    do {
        mpz_set(x, y);
        long int i = 0;
        while (i < r) {
            next_in_seq(y, y, N); // y = f(y)
            ++i;
        }
        k = 0;
        // Brent's
        do {
            long int j = 0;
            long int range = min(m, r-k);
            mpz_set(ys, y);
            while(j < range) {
                next_in_seq(y, y, N); // y = f(y)
                // q = q*|x-y|
                mpz_sub(temp_var, x, y);
                mpz_abs(temp_var, temp_var);
                mpz_mul(q, q, temp_var);
                mpz_mod(q, q, N);
                ++j;
            }
            mpz_gcd(d, q, N);
            k = k + m;
            ++count;
        } while ( k < r && mpz_cmp_si(d, 1) == 0 && count < limit);
        r = r*2;
    } while(count < limit && mpz_cmp_si(d, 1) <= 0);
    // This backtracking loop never finds any factors.
    // if (mpz_cmp(d, N) == 0) {
    //     count = 0;
    //     do {
    //         next_in_seq(ys, ys,  N);
    //         mpz_sub(temp_var, x, ys);
    //         mpz_abs(temp_var, temp_var);
    //         mpz_gcd(d, temp_var, N);
    //         ++count;
    //     } while (mpz_cmp_si(d,1) == 0 && count < limit);
    // }
    mpz_clear(x);
    mpz_clear(y);
    mpz_clear(q);
    mpz_clear(ys);
    mpz_clear(temp_var);
    // d == N
    if (mpz_cmp(d, N) == 0) {
        mpz_clear(d);
        return 0;
    }
    // d == 1 or d == 0, Invalid factors.
    else if (mpz_cmp_si(d, 1) == 0 || mpz_cmp_si(d, 0) == 0) {
        mpz_clear(d);
        return 0;
    }
    else {
        // Add factor to array
        mpz_set(factors[num_factors],d);
        // Update N
        mpz_divexact(N,N, d);
        mpz_clear(d);
        return 1;
    }
}

/*
 * Pseudo-random sequence generator for Pollard's rho algorithm
 */ 
void next_in_seq(mpz_t next, mpz_t prev, mpz_t N) {
    mpz_pow_ui(next, prev, 2); // X^2
    mpz_add_ui(next, next, C); // X^2 + C
    mpz_mod(next, next, N);    // (X^2 + 1) mod N
}

void print_factors(mpz_t factors[], int num_factors) {
    int i;
    for (i = 0; i < num_factors; ++i) {
        // This if-check is redundant
        if (factors[i] == 0) {
            break;
        }
        gmp_printf("%Zd\n", factors[i]);
    }
}

long int min(long int a, long int b) {
    return (a < b ? a : b);
}