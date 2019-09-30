#include "Diffusion.h"

#include "ExplicitScheme.h"

#include <iostream>
#include <cstdlib>

// For printf
#include <stdio.h>
// For timing
#include <time.h>

// OpenMP
#include <omp.h>


Diffusion::Diffusion(const InputFile* input, Mesh* m) :
    mesh(m) 
{

    std::string scheme_str = input->getString("scheme", "explicit");

    if(scheme_str.compare("explicit") == 0) {
        scheme = new ExplicitScheme(input, mesh);
    } else {
        std::cerr << "Error: unknown scheme \"" << scheme_str << "\"" << std::endl;
        exit(1);
    }

    subregion = input->getDoubleList("subregion", std::vector<double>());

    if (subregion.size() != 0 && subregion.size() != 4) {
        std::cerr << "Error:  region must have 4 entries (xmin, ymin, xmax, ymax)" << std::endl;
        exit(1);
    }

    init();
}

Diffusion::~Diffusion()
{
    delete scheme;
}

void Diffusion::init()
{
    double* u0 = mesh->getU0();

    int x_max = mesh->getNx()[0];
    int y_max = mesh->getNx()[1];

    double* cellx = mesh->getCellX();
    double* celly = mesh->getCellY();

    int nx = x_max+2;

//     TODO: FOR LOOP: DIFFUSION 1 - CLOCK ADDED - OMP INNER-dynamic OUTER
    if(!subregion.empty()) {

        double start_time = omp_get_wtime();

        #pragma omp parallel for schedule(dynamic, 100)
        for (int j = 0; j < y_max+2; j++) {

            #pragma omp parallel for schedule(dynamic, 100)
            for (int i = 0; i < nx; i++) {
                // Change: i < nx from i < xmax+2 <-- Does fewer additions
                if (celly[j] > subregion[1] && celly[j] <= subregion[3] && cellx[i] > subregion[0] && cellx[i] <= subregion[2]) {
                    u0[i + j * nx] = 10.0;
                } else {
                    u0[i + j * nx] = 0.0;
                }

            }
        }

        double time_spent = omp_get_wtime() - start_time;
        //printf("diff1:%f\n", time_spent);
        //printf("Diff1: Time elapsed is %f seconds\n", time_spent);

    } else {

        // TODO FOR LOOP: DIFFUSION 2 - CLOCK ADDED - OMP Inner-dynamic Outer
        double start_time = omp_get_wtime();

        #pragma omp parallel for schedule(dynamic, 100)
        for (int j = 0; j < y_max+2; j++) {

            #pragma omp parallel for schedule(dynamic, 100)
            for (int i = 0; i < x_max+2; i++) {
                u0[i+j*nx] = 0.0;

            }
        }

        double time_spent = omp_get_wtime() - start_time;
        //printf("diff2:%f\n", time_spent);
        //printf("Diff2: Time elapsed is %f seconds\n", time_spent);

    }

    scheme->init();
}

void Diffusion::doCycle(const double dt)
{
    scheme->doAdvance(dt);
}
