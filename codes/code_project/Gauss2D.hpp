#pragma once
#include "iostream"
#include "math.h"
#include "Gauss.hpp"



struct Gauss2D_Output
{
    double * xi;
    double * eta;
    double * w;
};


Gauss2D_Output * Gauss2D(int N1, int N2)
{
    // preallocation
    double * xi = new double[N1 * N2] ();
    double * eta = new double[N1 * N2] ();
    double * w = new double[N1 * N2] ();


    // generate 1D rule
    GaussOutput * G1Dop_1 = new GaussOutput;
    GaussOutput * G1Dop_2 = new GaussOutput;

    G1Dop_1 = Gauss(N1, -1, 1);
    G1Dop_2 = Gauss(N2, -1, 1);

    // generate the mesh grid of quadrature points
    for (int ii = 0; ii < N1; ++ii)
    {
        for (int jj = 0; jj < N2; ++jj)
        {
            xi[jj * N1 + ii] = G1Dop_1->x[ii];
            eta[jj * N1 + ii] = G1Dop_2->x[jj];
            w[jj * N1 + ii] = G1Dop_1->w[ii] * G1Dop_2->w[jj];
        }
    }

    Gauss2D_Output * G2Dop = new Gauss2D_Output;

    G2Dop->xi = xi;
    G2Dop->eta = eta;
    G2Dop->w = w;


    return G2Dop; 
}


