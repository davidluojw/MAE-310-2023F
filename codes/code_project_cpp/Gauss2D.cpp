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

    delete G2Dop;
    
    return 0;

}