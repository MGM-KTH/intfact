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
        //printf("here\n");
        if(mpz_probab_prime_p(N, 5)) {
            gmp_printf("%Zd\n", N);
            mpz_set_ui(N, 1);
            break;
        }
        result = pollards(N, factors, num_factors);
        //printf("result = %d\n", result);
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
    mpz_t x;
    mpz_init_set_ui(x, 2);

    mpz_t y;
    mpz_init_set_ui(y, 2);

    mpz_t d;
    mpz_init(d);

    long int count = 0;
    // Same limit for all numbers
    long int limit = 160000;
    // gmp_printf("Number of bits: %lu\n", mpz_sizeinbase(N,2));
    // if (mpz_sizeinbase(N,2) > 78) {
    //     return 0;
    // }
    // else if (mpz_sizeinbase(N,2) > 44) {
    //     limit = 120000;
    // } else if (mpz_sizeinbase(N,2) > 20) {
    //     limit = 90000;
    // } else {
    //     limit = 14000;
    // }
    // // gmp_printf("N is %Zd, limit is %Zd\n", N, limit);

    long int r = 1;
    long int m = 100; // TODO: set appropriate value. log(N) << m << N^(1/4) ?
    long int k = 0;
    mpz_t q;
    mpz_t temp_var;
    mpz_t ys;
    mpz_t power_of_two;
    mpz_init(q);
    mpz_init(temp_var);
    mpz_init(ys);
    mpz_init(power_of_two);

    // printf("limit is: %lu\n", limit);
    
    // Tough primes to test:
    // 38947539857394857374958374
    // 16055930642404191022652830
    // 23574372351791222645014404
    // 30655492778990882331786734
    // 1140925584794963113241931
    // 734286514387144056614055
    // 3979929034358890127857254
    // 30126842333649376717015439
    // 8355749255324197201318662
    while(count < limit && mpz_cmp_si(d, 1) <= 0) {
        mpz_set(x, y); // x = y
        long int i = 0;
        while (i < r) {
            next_in_seq(y, y, N); 
            ++i;
        }
        k = 0;
        long int j = 0;
        long int range = min(m, r-k);
        while ( k < r && mpz_cmp_si(d, 1) <= 0) {
            mpz_set(ys, y);
            while(j < range) {
                // y = f(y)    
                next_in_seq(y, y, N);
                // q = q*|x-y| mod N
                mpz_sub(temp_var, x, y);
                mpz_abs(temp_var, temp_var);
                mpz_mul(q, q, temp_var);
                mpz_mod(q, q, N);
                ++j;
            }
            mpz_gcd(d, q, N);
            k = k + m;
            ++count;
            //printf("k = %ld, r = %ld\n", k, r);
        }
        r = r*2;
        // gmp_printf("d is %Zd\n", d);
    }
    //gmp_printf("d is %Zd\n", d);
    if (mpz_cmp(d, N) == 0) {
        while (mpz_cmp_si(d,1) <= 0) {
            next_in_seq(ys, ys,  N);
            mpz_sub(temp_var, x, ys);
            mpz_abs(temp_var, temp_var);
            mpz_gcd(d, temp_var, N);
            if (mpz_cmp_si(d, 1) > 0) {
                break;
            }
        }
    }
    // Clear variables
    mpz_clear(x);
    mpz_clear(y);
    mpz_clear(q);
    mpz_clear(temp_var);

    if (mpz_cmp(d, N) == 0) {
        // printf("d==N, fail!\n");
        mpz_clear(d);
        return 0;
    }
    else if (mpz_cmp_si(d, 1) == 0 || mpz_cmp_si(d, 0) == 0) {
        // gmp_printf("d is trivial: d=%Zd\n",d);
        mpz_clear(d);
        return 0;
    }
    else {
        // printf("d is factor!!!!\n");
        mpz_set(factors[num_factors],d);
        // mpz_fdiv_q(N, N, d);
        mpz_divexact(N,N, d);
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
