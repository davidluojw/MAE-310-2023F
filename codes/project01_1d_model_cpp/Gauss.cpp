#include "iostream"
#include "math.h"




int Gauss(int N, int a, int b)
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

     for (int ii = 0; ii < N1; ++ii)
    {
        std::cout << "y[" << ii << "] = " << y[ii] << std::endl;
    }

    return 0;

   

    

}

int main()
{
    Gauss(5, 1, 1);

    return 0;
}