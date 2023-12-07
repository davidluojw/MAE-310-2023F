#include "Quad.hpp"
#include "Gauss2D.hpp"



int main()
{
    Gauss2D_Output * G2Dop = new Gauss2D_Output;

    G2Dop = Gauss2D(3, 4);

    for (int ii = 0; ii < 3 * 4; ++ii)
    {
        std::cout << "xi[" << ii << "] = " << G2Dop->xi[ii] << std::endl;
    }
    for (int ii = 0; ii < 3 * 4; ++ii)
    {
        std::cout << "eta[" << ii << "] = " << G2Dop->eta[ii] << std::endl;
    }
    for (int ii = 0; ii < 3 * 4; ++ii)
    {
        std::cout << "w[" << ii << "] = " << G2Dop->w[ii] << std::endl;
    }

    double *val = new double[12]();
    for (int ii = 0; ii < 12; ++ii)
    {
        val[ii] = Quad(4, G2Dop->xi[ii], G2Dop->eta[ii]);
    }

    for (int ii = 0; ii < 12; ++ii)
    {
        std::cout << "val[" << ii << "] = " << val[ii] << std::endl;
    }

    delete G2Dop;
    delete [] val;

    return 0;
}