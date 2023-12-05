#include "Gauss.hpp"
#include "PolyShape.hpp"

//--------------------------------------------------------
// Problem definition
// exact solution

double exact(double x)
{
    return pow(x, 5);
}

double exact_x(double x)
{
    return 5 * pow(x, 4);
}

double f(double x)
{
    return -20 * pow(x, 3);
}

double g = 1.0, h = 0.0;

//--------------------------------------------------------

int main()
{
    // Parameters of the FEM
    int n_el = 20;               // number of elements
    int n_en = 4;                // number of elements nodes
    int deg = n_en - 1;          // polynomial degree
    int n_np = n_el * deg + 1;   // number of points
    int n_eq = n_np - 1;         // number of equations
    int n_int = 3;               // number of quadrature points     

    //--------------------------------------------------------
    // Generate the mesh
    // nodal coordinates
    double hh = 1 / (double)n_el;
    double x_coor[n_np];

    for (int ii = 0; ii < n_np; ++ii)
    {
        x_coor[ii] = ii * hh / deg;
    }
    // for (int ii = 0; ii < n_np; ++ii)
    // {
    //     std::cout << "x_coor[" << ii << "] = " << x_coor[ii] << std::endl;
    // }


    return 0;
}

