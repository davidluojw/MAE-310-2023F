#include "iostream"
#include "Gauss.hpp"
#include "PolyShape.hpp"



int main()
{
    GaussOutput * Gop = new GaussOutput;

    Gop = Gauss(5, -1, 1);

    for (int ii = 0; ii < 5; ++ii)
    {
        std::cout << "w[" << ii << "] = " << Gop->w[ii] << std::endl;
    }
    for (int ii = 0; ii < 5; ++ii)
    {
        std::cout << "xi[" << ii << "] = " << Gop->x[ii] << std::endl;
    }

    double val;



    val = Polyshape(3, 1, Gop->x[4], 0);

    std::cout << "val = " << val << std::endl;

    delete Gop;

    return 0;
}