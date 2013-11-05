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
    // TODO: We could clear factors and current_N here, but it's pointless, right?
    return 0;
}

/*
 * Factorizes N using Pollard's rho method
 */
void factorize(mpz_t N, mpz_t factors[]) {
    int num_factors;
    num_factors = find_trivial_factors(N, factors);

    int result;
    while (mpz_cmp_si(N, 1) != 0) {
        if(mpz_probab_prime_p(N, 5)) {
            gmp_printf("%Zd\n", N);
            mpz_set_ui(N, 1);
            break;
        }
        result = pollards(N, factors, num_factors);
        if (!result) {
            break;
        }
        ++num_factors;
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
    
    // Initialize random number container
    mpz_t xi_last;
    mpz_init_set_ui(xi_last, 1);

    // Initialize and seed a randstate
    // gmp_randstate_t rand_state;
    // gmp_randinit_default(rand_state);
    // gmp_randseed_ui(rand_state, time(NULL));

    // Get random number
    // mpz_urandomm(xi_last, rand_state, N);

    mpz_t x2i_last;
    // mpz_init(x2i_last);
    // next_in_seq(x2i_last, xi_last, N);
    mpz_init_set_ui(x2i_last, 2);

    mpz_t xi;
    mpz_t x2i;
    mpz_t diff;
    mpz_init(xi);
    mpz_init(x2i);
    mpz_init(diff);

    mpz_t d;
    mpz_init(d);
    // Pollard's + Brent's improvement. Only calculate GCD every now and then on the product of the diffs.
    // mpz_t gMul;
    // mpz_init_set_ui(gMul, 1);
    // mpz_t r;
    // mpz_init(r);

    long int count = 0;
    long int limit;
    // gmp_printf("Number of bits: %lu\n", mpz_sizeinbase(N,2));
    if (mpz_sizeinbase(N,2) > 84) {
        return 0;
    }
    else if (mpz_sizeinbase(N,2) > 44) {
        limit = 140000;
    } else if (mpz_sizeinbase(N,2) > 20) {
        limit = 100000;
    } else {
        limit = 20000;
    }
    // gmp_printf("N is %Zd, limit is %Zd\n", N, limit);

    long int r = 1;
    long int m = 1000; // TODO: set appropriate value. log(N) << m << N^(1/4) ?
    long int k = 0;
    mpz_t q;
    mpz_init(q);
    mpz_t temp_x2i;
    mpz_init(temp_x2i);
    mpz_t temp_sub;
    mpz_init(temp_sub);

    mpz_t temp_abs;
    mpz_init(temp_abs);

    // printf("limit is: %lu\n", limit);
    while(count < limit && mpz_cmp_si(d, 1) > 0) {
        mpz_set(xi_last, x2i_last); // x = y
        long int i = 1;
        while (i < r) {
            next_in_seq(x2i_last, x2i_last, N); 
            ++i;
            k = 0; // why not outside this loop?
        }
        long int j = 1;
        long int range = min(m, r-k);
        while(j < range) {
            // y = f(y)    
            next_in_seq(x2i_last, x2i_last, N);
            // q = q*|x-y| mod N
            mpz_sub(temp_sub, xi_last, x2i_last);
            mpz_set(temp_abs, temp_sub);
            mpz_abs(temp_abs, temp_abs);
            mpz_mul(q, q, temp_abs);
            mpz_mod(q, q, N);
        }
        mpz_gcd(d, q, N);
        k = k + m;
        ++j;
        if ( k < r && mpz_cmp_si(d, 1) > 0) {
            r = r*2;
        }
        ++count;
    }
    if (mpz_cmp(d, N) == 0) {
        while(1) {
            next_in_seq(temp_x2i, temp_x2i,  N);
            mpz_sub(xi_last, xi_last, temp_x2i);
            mpz_set(q, xi_last);
            mpz_gcd(d, q, N);
            if (mpz_cmp_si(d, 1) > 0) {
                break;
            }
        }
    }
    // Clear variables
    mpz_clear(xi_last);
    mpz_clear(x2i_last);
    mpz_clear(xi);
    mpz_clear(x2i);
    mpz_clear(diff);
    mpz_clear(q);
    mpz_clear(temp_abs);
    mpz_clear(temp_x2i);
    if (mpz_cmp(d, N) == 0) {
        mpz_clear(d);
        return 0;
    }
    else if (mpz_cmp_si(d, 0) == 0) {
        mpz_clear(d);
        return 0;
    }
    else {
        mpz_set(factors[num_factors],d);
        mpz_fdiv_q(N, N, d);
        mpz_clear(d);
        return 1;
    }
}

void next_in_seq(mpz_t next, mpz_t prev, mpz_t N) {
    mpz_pow_ui(next, prev, 2); // X^2
    mpz_add_ui(next, next, 3); // X^2 + 1
    mpz_mod(next, next, N);    // (X^2 + 1) mod N
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

long int min(long int a, long int b) {
    return (a < b ? a : b);
}

