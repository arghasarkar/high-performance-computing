//
// Created by ArghaWin10 on 03/02/2019.
//


// For printf
#include <stdio.h>
// For timing
#include <time.h>

// OpenMP
#include <omp.h>

int main() {

    clock_t begin = clock();

    int n = 2140000000;

    long sum = 0;
    int i = 0;
    int _temp = 0;

    #pragma omp parallel for reduction (+:sum) private(_temp)
    for (i = 0; i < n; i++) {
        _temp = i;
        sum = sum + i;
    }

    clock_t end = clock();

    printf("The sum is %ld.\n", sum);

    double time_taken = (double) (end - begin) / CLOCKS_PER_SEC;

    printf("The time taken is %f.\n", time_taken);

}