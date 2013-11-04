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
    // gmp_printf("Factorizing %Zd\n", N);
    int num_factors;

    if(mpz_probab_prime_p(N, 25) > 0) {
        gmp_printf("%Zd\n", N);
        printf("\n");
        return;
    }

    num_factors = find_trivial_factors(N, factors);
    
    // If N != 1 (not fully factorized)
    int result;
    while(mpz_cmp_si(N, 1) != 0) {    
        result = pollards(N, factors, num_factors);
        if (!result) {
            break;
        }
        ++num_factors;

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
    mpz_init(xi_last);

    // Initialize and seed a randstate
    gmp_randstate_t rand_state;
    gmp_randinit_default(rand_state);
    gmp_randseed_ui(rand_state, time(NULL));

    // Get random number
    mpz_urandomm(xi_last, rand_state, N);

    mpz_t x2i_last;
    mpz_init(x2i_last);
    next_in_seq(x2i_last, xi_last, N);

    mpz_t xi;
    mpz_t x2i;
    mpz_t diff;
    mpz_init(xi);
    mpz_init(x2i);
    mpz_init(diff);

    mpz_t d;
    mpz_init(d);

    mpz_t count;
    mpz_init_set_ui(count, 0);
    mpz_t limit;
    mpz_init(limit);
    if(mpz_sizeinbase(N,2) > 20) {
        printf("hit limit\n");
        mpz_set_ui(limit, 1000000);
    }else{
        mpz_sqrt(limit, N);
    }
    gmp_printf("N is %Zd, limit is %Zd\n", N, limit);

    while(mpz_cmp(count, limit)<0) {
        next_in_seq(xi, xi_last, N);

        // TODO: Same as above. Next 2i is simply next-next-in-seq? i -> i+1 and 2i -> 2i+2?
        next_in_seq(x2i, x2i_last, N);
        next_in_seq(x2i, x2i, N);

        mpz_sub(diff, x2i, xi);
        mpz_gcd(d, diff, N);

        if(mpz_cmp_si(d, 1) > 0) {
            // gmp_printf("factor found: %Zd\n", d);
            mpz_set(factors[num_factors],d);
            mpz_fdiv_q(N, N, d);
            mpz_clear(xi_last);
            mpz_clear(x2i_last);
            mpz_clear(xi);
            mpz_clear(x2i);
            mpz_clear(diff);
            return 1; 
        }
        // gmp_printf("numbers: xi = %Zd, x2i = %Zd\n", xi, x2i);
        // gmp_printf("d = %Zd\n", d);
        mpz_set(xi_last, xi);
        mpz_set(x2i_last, x2i);
        mpz_add_ui(count,count,1);
    }

    //gmp_printf("number N: %Zd random number %Zd\n", N, rand);

    // TODO? Call twice for prev 2i and once for prev i.

    // TDOO: decide how to loop.

    // Clear variables
    mpz_clear(xi_last);
    mpz_clear(x2i_last);
    mpz_clear(xi);
    mpz_clear(x2i);
    mpz_clear(diff);
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
