#include "Quad_grad.hpp"
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

    double * val_xi = new double[12] ();
    double * val_eta = new double[12] ();
    Quad_grad_Oput * QdGdOp = new Quad_grad_Oput;
    for (int ii = 0; ii < 12; ++ii)
    {
        QdGdOp = Quad_grad(4, G2Dop->xi[ii], G2Dop->eta[ii]);
        val_xi[ii] = QdGdOp->val_xi;
        val_eta[ii] = QdGdOp->val_eta;
    }

    for (int ii = 0; ii < 12; ++ii)
    {
        std::cout << "val_xi[" << ii << "] = " << val_xi[ii] << std::endl;
    }
    for (int ii = 0; ii < 12; ++ii)
    {
        std::cout << "val_eta[" << ii << "] = " << val_eta[ii] << std::endl;
    }

    return 0;
}