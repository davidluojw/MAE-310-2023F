#include "iostream"
#include "Gauss.hpp"
#include "PolyShape.hpp"



int main()
{
    GaussOutput * Gop = new GaussOutput;

    Gop = Gauss(5, -1, 1);

    double val;

    val = Polyshape(2, 1, Gop->x[4], 1);

    std::cout << "val = " << val << std::endl;

    delete Gop;

    return 0;
}