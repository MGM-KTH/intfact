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

    if(mpz_probab_prime_p(N, 5)) {
        gmp_printf("%Zd\n", N);
        mpz_set_ui(N, 1);
    }else if(mpz_sizeinbase(N,2) < 70){

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
    // Ugly fix for ceil(sqrt(N)), does trunc_sqrt(N-1) + 1 instead
    mpz_sub_ui(a, N, 1);
    mpz_sqrt(a, a);
    mpz_add_ui(a, a, 1);

    mpz_t b2;
    mpz_init(b2);
    mpz_mul(b2, a, a);
    mpz_sub(b2, b2, N);

    mpz_t count;
    mpz_init_set_ui(count, 0);
    mpz_t limit;
    mpz_init(limit);
    if (mpz_sizeinbase(N,2) > 40) {
        mpz_set_ui(limit, 100000);
    } else if (mpz_sizeinbase(N,2) > 20) {
        mpz_set_ui(limit, 90000);
    } else {
        mpz_set_ui(limit, 20000);
    }
    //gmp_printf("N is %Zd, limit is %Zd\n", N, limit);

    while(mpz_cmp(count, limit)<0) {
        //gmp_printf("a is %Zd, b2 is %Zd\n", a, b2);
        //gmp_printf("N is %Zd\n", N);

        if(mpz_perfect_square_p(b2) != 0) {
            //gmp_printf("b2 found: %Zd\n", b2);
            mpz_sqrt(b2, b2);
            mpz_sub(b2, a, b2);
            mpz_set(factors[num_factors],b2);
            mpz_fdiv_q(N, N, b2);

            mpz_clear(a);
            mpz_clear(b2);
            mpz_clear(count);
            mpz_clear(limit);
            return 1; 
        }
        // gmp_printf("numbers: xi = %Zd, x2i = %Zd\n", xi, x2i);
        // gmp_printf("d = %Zd\n", d);
        mpz_add_ui(a, a, 1);
        mpz_mul(b2, a, a);
        mpz_sub(b2, b2, N);
        mpz_add_ui(count,count,1);
    }

    //gmp_printf("number N: %Zd random number %Zd\n", N, rand);

    // Clear variables
    mpz_clear(a);
    mpz_clear(b2);
    mpz_clear(count);
    mpz_clear(limit);
    return 0;
}

void next_in_seq(mpz_t next, mpz_t prev, mpz_t N) {
    mpz_pow_ui(next, prev, 2); // X^2
    mpz_add_ui(next, next, 1); // X^2 + 1
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
