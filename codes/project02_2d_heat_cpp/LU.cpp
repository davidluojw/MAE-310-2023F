#include "LU.hpp"


int main()
{
    // example 1
    // double * A = new double[4] ();
    // A[0] = 0.001; A[1] = 1.00; A[2] = 1.00; A[3] = 2.00;
    // int * pp = new int[2];
    // for (int ii = 0; ii < 2; ++ii)
    // {
    //     pp[ii] = ii;
    // }
    // LU_var * lu_var = new LU_var;

    // lu_var = LU_fac(A, 2, pp);

    // for (int ii = 0; ii < 4; ++ii)
    // {
    //     std::cout << "mat[" << ii << "] = " << lu_var->mat[ii] << std::endl;
    // }
    // for (int ii = 0; ii < 2; ++ii)
    // {
    //     std::cout << "pp[" << ii << "] = " << lu_var->pp[ii] << std::endl;
    // }

    // double * b = new double [2];
    // b[0] = 1.00; b[1] = 3.00;
    // double * x = new double[2] ();
    // x = LU_solve(b, 2, lu_var->pp, lu_var->mat);
    // for (int ii = 0; ii < 2; ++ii)
    // {
    //     std::cout << "x[" << ii << "] = " << x[ii] << std::endl;
    // }


    // example 2
    double * A = new double[16] ();
    A[0] = 4.0; A[1] = -2.0; A[2] = 4.0; A[3] = 2.0; 
    A[4] = -2.0; A[5] = 10.0; A[6] = -2.0; A[7] = -7.0; 
    A[8] = 4.0; A[9] = -2.0; A[10] = 8.0; A[11] = 4.0;
    A[12] = 2.0; A[13] = -7.0; A[14] = 4.0; A[15] = 7.0;
    int * pp = new int[4];
    for (int ii = 0; ii < 4; ++ii)
    {
        pp[ii] = ii;
    }
    LU_var * lu_var = new LU_var;

    lu_var = LU_fac(A, 4, pp);

    for (int ii = 0; ii < 16; ++ii)
    {
        std::cout << "mat[" << ii << "] = " << lu_var->mat[ii] << std::endl;
    }
    for (int ii = 0; ii < 4; ++ii)
    {
        std::cout << "pp[" << ii << "] = " << lu_var->pp[ii] << std::endl;
    }

    double * b = new double [4];
    b[0] = 8.0; b[1] = 2.0; b[2] = 16.0; b[3] = 6.0;
    double * x = new double[4] ();
    x = LU_solve(b, 4, lu_var->pp, lu_var->mat);
    for (int ii = 0; ii < 4; ++ii)
    {
        std::cout << "x[" << ii << "] = " << x[ii] << std::endl;
    }


    
    return 0;
}