#include "main.h"

#define MAX_LINE_LENGTH 256
#define FAIL_STRING "fail"

int main() {
    static char buffer[MAX_LINE_LENGTH];
    const char *line;
    while(( line = fgets(buffer, sizeof(buffer), stdin)) != NULL) {
        mpz_t current_N;
        mpz_init(current_N);
        mpz_set_str(current_N, line, 10);
        factorize(current_N);
    }
    return 0;
}

/*
 * Factorizes N using Pollard's rho method
 */
void factorize(mpz_t N) {
    // gmp_printf("Factorizing %Zd\n", N);
    
    
    
    // printf("%s\n", FAIL_STRING); // Print this if you can't factorize
    printf("\n");
}