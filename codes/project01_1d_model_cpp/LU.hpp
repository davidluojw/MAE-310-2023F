#pragma once
#include "iostream"


struct LU_var
{
    double *mat;
    int * pp;
    bool is_fac;
};



// ----------------------------------------------------------------------
// perform LU-factorization for the matrix. The mat object will be replaced
// by the LU matrices. Only partial pivoting is performed. Complete pivoting
// is not used because the improvement of stability is marginal and the
// amount of time needed will increase.
// ----------------------------------------------------------------------
LU_var * LU_fac(double * mat, int N, int * pp)
{
    bool is_fac = false;
    for(int kk=0; kk<N-1; ++kk)
    {
        double max_value = std::abs(mat[kk*N+kk]);
        int max_index = kk;
        bool pivot_flag = false;
        for(int ii=kk+1; ii<N; ++ii)
        {
            if( max_value < std::abs(mat[ii*N+kk]) )
            {
                max_value = std::abs(mat[ii*N+kk]);
                max_index = ii;
                pivot_flag = true;
            }
        }

        if(pivot_flag)
        {
            std::swap( pp[kk] , pp[max_index] );

            for(int ii=0; ii<N; ++ii)
                std::swap( mat[kk*N+ii] , mat[max_index*N+ii] );
        }

        const double invAkk = 1.0 / mat[kk*N+kk];

        for(int ii=kk+1; ii<N; ++ii)
        {
            mat[ii*N+kk] = mat[ii*N+kk] * invAkk;
            for(int jj=kk+1; jj<N; ++jj)
                mat[ii*N+jj] -= mat[ii*N+kk] * mat[kk*N+jj];
        }
    }

    is_fac = true;

    LU_var * lu_var = new LU_var;
    lu_var->is_fac = is_fac;
    lu_var->mat = mat;
    lu_var->pp = pp;

    return lu_var;
}



// ----------------------------------------------------------------------
// with LU factorization performed, solve a linear problem with given RHS
// users are responsible for allocating the b and x arrays.
// ----------------------------------------------------------------------
double * LU_solve( double * bb, int N, int * pp, double * mat)
{
    double * xx = new double[N] ();
    for(int ii=0; ii<N; ++ii) xx[ii] = bb[pp[ii]];

    for(int ii=1; ii<N; ++ii)
        for(int jj=0; jj<ii; ++jj)
            xx[ii] -= mat[ii*N+jj] * xx[jj];

    for(int ii=N-1; ii>=0; --ii)
    {
        for(int jj=N-1; jj>ii; --jj)
            xx[ii] -= mat[ii*N+jj] * xx[jj];

        xx[ii] = xx[ii] / mat[ii*N+ii];
    }

    return xx;
}