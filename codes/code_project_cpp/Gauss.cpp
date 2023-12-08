#include "iostream"
#include "math.h"
#include "Gauss.hpp"

int main()
{
    GaussOutput * Gop = new GaussOutput;

    Gop = Gauss(5, -1, 1);

    for (int ii = 0; ii < 5; ++ii)
    {
        std::cout << "x[" << ii << "] = " << Gop->x[ii] << std::endl;
        std::cout << "w[" << ii << "] = " << Gop->w[ii] << std::endl;
    }

    delete Gop;

    return 0;
}