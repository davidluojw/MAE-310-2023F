#pragma once
#include "iostream"


struct Quad_grad_Oput
{
    double val_xi;
    double val_eta;
};

Quad_grad_Oput * Quad_grad(int aa, double xi, double eta)
{
    Quad_grad_Oput * QdGdOp = new Quad_grad_Oput;

    if (aa == 1)
    {
        QdGdOp->val_xi = -0.25 * (1.0 - eta);
        QdGdOp->val_eta = -0.25 * (1.0 - xi);
    }
    else if (aa == 2)
    {
        QdGdOp->val_xi = 0.25 * (1.0 - eta);
        QdGdOp->val_eta = -0.25 * (1.0 + xi);
    }
    else if (aa == 3)
    {
        QdGdOp->val_xi = 0.25 * (1.0 + eta);
        QdGdOp->val_eta = 0.25 * (1.0 + xi);
    }
    else if (aa == 4)
    {
        QdGdOp->val_xi = -0.25 * (1.0 + eta);
        QdGdOp->val_eta = 0.25 * (1.0 - xi);
    }
    else fprintf(stderr, "Error: degree has to be 1, 2, 3 or 4.");


    return QdGdOp;
}
