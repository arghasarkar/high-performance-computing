#include <iostream>
#include <string>
#include <cstdlib>

#include "InputFile.h"
#include "Driver.h"

// For printf
#include <stdio.h>
// For timing
#include <time.h>

// OpenMP
#include <omp.h>


int main(int argc, char *argv[])
{

    // set the number of threads
    omp_set_num_threads(8);

    if (argc != 2) {
        std::cerr << "Usage: deqn <filename>" << std::endl;
        exit(1);
    }


    const char* filename = argv[1];
    InputFile input(filename);

    std::string problem_name(filename);

    int len = problem_name.length();

    if(problem_name.substr(len - 3, 3) == ".in")
        problem_name = problem_name.substr(0, len-3);

    // Strip out leading path
    size_t last_sep = problem_name.find_last_of("/");

    if (last_sep != std::string::npos) {
        last_sep = last_sep + 1;
    } else {
        last_sep = 0;
    }

    problem_name = problem_name.substr(last_sep, problem_name.size());

    Driver driver(&input, problem_name);

    // CLOCK FOR OVER ALL LOOP
    double total_time_wtime = 0.0;
    double start_wtime = omp_get_wtime();

    // Main line of code
    int argha = 0;
    int max_runs = 1;
    for (argha = 0; argha < max_runs; argha++) {
        driver.run();
    }

    double end_wtime = omp_get_wtime();

    total_time_wtime = end_wtime - start_wtime;
    printf("Time taken: %f.\n", total_time_wtime);

    return 0;
}
