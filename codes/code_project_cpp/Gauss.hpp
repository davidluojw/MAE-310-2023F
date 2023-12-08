#pragma once
#include "iostream"
#include "math.h"

struct GaussOutput
{
    double * x;
    double * w;

};



GaussOutput * Gauss(int N, int a, int b)
{
    long double PI = 3.14159265358979323846264338327950288419716939937510582;
    N = N - 1;
    int N1 = 0, N2 = 0;
    
    N1 = N + 1;
    N2 = N + 2;


    double xu[N1];
    int lin[N + 1];

    double dxu = 2.0 / (N1 - 1);
    xu[0] = -1.0; xu[N1 - 1] = 1.0;
    for (int ii = 1; ii < N1 - 1; ++ii)
    {
        xu[ii] = ii * dxu - 1;
    }


    for (int ii = 0; ii < N + 1; ++ii)
    {
        lin[ii] =  ii;
    }

    // Initial guess
    double y[N1];
    for (int ii = 0; ii < N1; ++ii)
    {
        y[ii] = cos((2*lin[ii] + 1) * PI / (2 * N + 2)) + (0.27 / N1) * sin(PI * xu[ii] * N / N2);
    }



    // Legendre-Gauss Vandermore Matrix
    double L[N1 * N2];

    for (int ii = 0; ii < N1; ++ii)
    {
        for (int jj = 0; jj < N2; ++jj)
        {
            L[ii * N2 + jj] = 0;
        }
    }


    

    // Derivative of LGVM
    double Lp[N1 * N2];

    for (int ii = 0; ii < N1; ++ii)
    {
        for (int jj = 0; jj < N2; ++jj)
        {
            Lp[ii * N2 + jj] = 0;
        }
    }
    
    // Compute the zeros of the N+1 Legendre Polynomial
    // using the recursion relation and the Newton-Raphson method

    int y0 = 2;

    // Iterate until new points are uniformly within epsilon of old points
    double max_value = 0.0;
    for (int ii = 0; ii < N1; ++ii)
    {
        if (abs(y[ii] - y0) > max_value)
        {
            max_value = abs(y[ii] - y0);
        }
    }

    // std::cout << "max_value = " << max_value << std::endl;

   double Lpp[N1];
    


    while (max_value > __DBL_EPSILON__)
    {
         for (int ii = 0; ii < N1; ++ii)
        {
            L[ii * N2] = 1;
        }

        for (int ii = 0; ii < N1; ++ii)
        {
            Lp[ii * N2] = 0;
        }

        for (int ii = 0; ii < N1; ++ii)
        {
            L[ii * N2 + 1] = y[ii];
        }

        for (int ii = 0; ii < N1; ++ii)
        {
            Lp[ii * N2 + 1] = 1;
        }

        for (int kk = 1; kk < N1; ++kk)
        {
            for (int ii = 0; ii < N1; ++ii)
            {
                L[ii * N2 + (kk + 1)] =( (2 * (kk+1) - 1) * y[ii] * L[ii * N2 + kk] - kk * L[ii * N2 + (kk - 1)] ) / (kk + 1);
            }
        }

        

        for (int ii = 0; ii < N1; ++ii)
        {
            Lpp[ii] = (N2) * (L[ii * N2 + (N1 - 1)] - y[ii] * L[ii * N2 + (N2 - 1)]) / (1.0 - y[ii] * y[ii]);
        }


        

        double y00[N1];

        for (int ii = 0; ii < N1; ++ii)
        {
            y00[ii] = y[ii];
        }

        for (int ii = 0; ii < N1; ++ii)
        {
            y[ii] = y00[ii] - L[ii * N2 + (N2 - 1)] / Lpp[ii];
        }

        max_value = 0.0;

        for (int ii = 0; ii < N1; ++ii)
        {
            if (abs(y[ii] - y00[ii]) > max_value)
            {
                max_value = abs(y[ii] - y00[ii]);
            }
        }

    }
    
    GaussOutput * Gop = new GaussOutput;
    

    // Linear map from [-1, 1] to [a,b]
    Gop->x = new double[N1];
    for (int ii =0; ii < N1; ++ii)
    {
        Gop->x[ii] = (a * (1 - y[ii]) + b * (1 + y[ii])) / 2.0; 
    }

    // Compute the weights
    Gop->w = new double[N1];
    for (int ii = 0; ii < N1; ++ii)
    {
        Gop->w[ii] = (b-a) / ((1 - y[ii] * y[ii]) * Lpp[ii] * Lpp[ii]) * ((double)N2 / (double)N1) * ((double)N2 / (double)N1);
    }
  

    return Gop;

}
