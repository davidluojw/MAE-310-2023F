#include "iostream"
#include "math.h"



int Gauss(int N, int a, int b)
{
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
    
    

    for (int ii = 0; ii < N1; ++ii)
    {
        std::cout << "lin[" << ii << "] = " << lin[ii] << std::endl;
    }

    for (int ii = 0; ii < N1; ++ii)
    {
        std::cout << "xu[" << ii << "] = " << xu[ii] << std::endl; 
    }

    return 0;

    // // Initial guess
    // double y[N1];
    // y = cos( (2 * lin) )

}

int main()
{
    Gauss(5, 1, 1);

    return 0;
}