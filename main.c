#include "main.h"

int main() {
    mpz_t input[100];
    printf("array size %lu\n", sizeof(*input)/sizeof((*input)[0]));
    readInput(&input);
    gmp_printf("read %Zd\n", input[0]);
    printf("success");
    return 0;
}

void readInput(mpz_t (*array)[100]) {
    gmp_scanf("%Zd", array[0]);
}
