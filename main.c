#include "main.h"

int main() {
    mpz_t input[100];
    printf("array size %lu\n", sizeof(*input)/sizeof((*input)[0]));

    int input_count = readInput(&input);

    int i;
    for(i = 0; i < input_count; i++) {
        gmp_printf("read %Zd\n", input[i]);
    }

    printf("success");
    return 0;
}

int readInput(mpz_t (*array)[100]) {
    mpz_init((*array)[0]);
    int count = 0;
    printf("count is %d\n", count);
    while(count < 100 && gmp_scanf("%Zd", &((*array)[count])) != EOF) {
        count++;
        mpz_init((*array)[count]);
        //printf("count is %d\n", count);
        //printf("array size %lu\n", sizeof(*array)/sizeof((*array)[0]));
    }
    return count;
}
