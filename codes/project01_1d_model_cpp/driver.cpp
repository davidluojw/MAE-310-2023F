#include "Gauss.hpp"
#include "PolyShape.hpp"
#include "LU.hpp"

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
    return -20.0 * pow(x, 3);
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


    // IEN
    int IEN[n_en * n_el];

    for (int ee = 0; ee < n_el; ++ee)
    {
        for (int aa = 0; aa < n_en; ++aa)
        {
            IEN[aa * n_el + ee] = ee * deg + (aa+1);
        }
    }
    // for (int ii = 0; ii < n_en * n_el; ++ii)
    // {
    //     std::cout << "IEN[" << ii << "] = " << IEN[ii] << std::endl;
    // }

    //--------------------------------------------------------

    // ID and LM arrays are generated based on the BC info
    int ID[n_np];     // Space dimension is 1D
    for (int ii = 0; ii < n_np; ++ii)
    {
        ID[ii] = ii + 1;
    }
    ID[n_np - 1] = 0; // Modify ID according to the Dirichlet BC info
    // for (int ii = 0; ii < n_np; ++ii)
    // {
    //     std::cout << "ID[" << ii << "] = " << ID[ii] << std::endl;
    // }

    double LM[1 * n_en * n_el];
    for (int ee = 0; ee < n_el; ++ee)
    {
        for (int aa = 0; aa < n_en; ++aa)
        {
            LM[aa * n_el + ee] = ID[IEN[aa * n_el + ee] - 1];
        }
    }
    // for (int ii = 0; ii < n_en * n_el; ++ii)
    // {
    //     std::cout << "LM[" << ii << "] = " << LM[ii] << std::endl;
    // }

    //--------------------------------------------------------

    // generate the quadrature rule
    GaussOutput * Gop = new GaussOutput;
    Gop = Gauss(n_int, -1, 1);
    // for (int ii = 0; ii < n_int; ++ii)
    // {
    //     std::cout << "w[" << ii << "] = " << Gop->w[ii] << std::endl;
    // }
    // for (int ii = 0; ii < n_int; ++ii)
    // {
    //     std::cout << "xi[" << ii << "] = " << Gop->x[ii] << std::endl;
    // }

    //--------------------------------------------------------

    // Assembly of K and F
    double *K = new double[n_eq * n_eq] ();   // allocate the global stiffness matrix
    double *F = new double[n_eq] ();          // allocate the global load vector

    
    
    for (int ee =0; ee < n_el; ++ee)
    {
        double *k_e = new double[n_en * n_en] ();
        double *f_e = new double[n_en] ();
        double *x_ele = new double[n_en] ();
        
        // Global coordinates of x
        for (int aa = 0; aa < n_en; ++aa)
        {
            x_ele[aa] = x_coor[IEN[aa * n_el + ee] - 1];
        } 
        // for (int ii = 0; ii < n_en; ++ii)
        // {
        //     std::cout << "x_ele[" << ii << "] = " << x_ele[ii] << std::endl;
        // }

        // Compute every Gauss points
        for (int ll = 0; ll < n_int; ++ll)
        {
            double dx_dxi = 0.0;
            double x_l = 0.0;

            // push forward from parent frame Gop->x[ll] to physical frame x_l
            for (int aa = 0; aa < n_en; ++aa)
            {
                dx_dxi = dx_dxi + x_ele[aa] * Polyshape(deg, aa + 1, Gop->x[ll], 1);
                x_l = x_l + x_ele[aa] * Polyshape(deg, aa + 1, Gop->x[ll], 0);
            }
            // std::cout << "dx_dxi = " << dx_dxi << std::endl;
            // std::cout << "x_l = " << x_l << std::endl;
            

            // inverse of dx_dxi
            double dxi_dx = 1.0 / dx_dxi;

            // generate element stiffness matrix and element force vector
            for (int aa = 0; aa < n_en; ++aa)
            {
                f_e[aa] = f_e[aa] + Gop->w[ll] * f(x_l) * Polyshape(deg, aa + 1, Gop->x[ll], 0) * dx_dxi; 
                for (int bb = 0; bb < n_en; ++bb)
                {
                    k_e[aa * n_en + bb] = k_e[aa * n_en + bb] + Gop->w[ll] * Polyshape(deg, aa + 1, Gop->x[ll], 1) * Polyshape(deg, bb + 1, Gop->x[ll], 1) * dxi_dx;
                }
            }
        }
        // for (int ii = 0; ii < n_en * n_en; ++ii)
        // {
        //     std::cout << "k_e [" << ii << "] = " << k_e[ii] << std::endl;
        // }
        // for (int ii = 0; ii < n_en; ++ii)
        // {
        //     std::cout << "f_e [" << ii << "] = " << f_e[ii] << std::endl;
        // }

        // Now we need to put element k and f into global K and F
        for (int aa = 0; aa < n_en; ++aa)
        {
            int PP;
            PP = LM[aa * n_el + ee];
            if (PP > 0)
            {
                F[PP - 1] = F[PP - 1] + f_e[aa];
                for (int bb = 0; bb < n_en; ++bb)
                {
                    int QQ;
                    QQ = LM[bb * n_el + ee];
                    if (QQ > 0)
                    {
                        K[(PP - 1) * n_eq + QQ - 1] = K[(PP - 1) * n_eq + QQ - 1] + k_e[aa * n_en + bb];
                    } 
                    else 
                    {
                        F[PP - 1] = F[PP - 1] - k_e[aa * n_en + bb] * g;
                    }
                }
            }

        }

        if (ee == 0)
        {
            F[ID[IEN[ 0 * n_el + ee] - 1] - 1] = F[ID[IEN[ 0 * n_el + ee] - 1] - 1] + h;
        }



        delete [] k_e;
        delete [] f_e;
        delete [] x_ele;
    }

    // for (int ii = 0; ii < n_eq * n_eq; ++ii)
    // {
    //     if (K[ii] != 0)
    //     {
    //         std::cout << "ii = " << ii << std::endl;
    //         int a = ii / n_eq + 1, b = ii % n_eq + 1;
    //         std::cout << "K [" << a << ", " << b << "] = " << K[ii] << std::endl;
    //     }
            
    // }
    // for (int ii = 0; ii < n_eq; ++ii)
    // {
    //     std::cout << "F [" << ii << "] = " << F[ii] << std::endl;
    // }

    //--------------------------------------------------------

    // Now we have K and F assembled and we solve the linear system Kd = F
    // pp: permutation infomation of row of equation Kd=F generated from LU-fac
    // initialize:
    int * pp = new int[n_eq];
    for (int ii = 0; ii < n_eq; ++ii)
    {
        pp[ii] = ii;
    }

    // LU factorization
    LU_var * lu_var = new LU_var;
    lu_var = LU_fac(K, n_eq, pp);
    // for (int ii = 0; ii < n_eq * n_eq; ++ii)
    // {
    //     std::cout << "mat[" << ii << "] = " << lu_var->mat[ii] << std::endl;
    // }
    // for (int ii = 0; ii < n_eq; ++ii)
    // {
    //     std::cout << "pp[" << ii << "] = " << lu_var->pp[ii] << std::endl;
    // }

    // LU solve, solution x array
    double * x = new double[n_eq] ();
    x = LU_solve(F, n_eq, lu_var->pp, lu_var->mat);
    // for (int ii = 0; ii < n_eq; ++ii)
    // {
    //     std::cout << "x[" << ii << "] = " << x[ii] << std::endl;
    // }

    // --------------------------------------------------------

    // Generate the full solution vector by inserting back the Dirichlet value
    double *disp = new double[n_np]();
    for (int ii = 0; ii < n_eq; ++ii)
    {
        disp[ii] = x[ii];
    }
    disp[n_np - 1] = g;

    //-------------------------------------------------------
    
    // Calculate the error
    int nqp = 10;   // we need more points
    Gop = Gauss(nqp, -1, 1);

    // errors
    double L2_top = 0.0; double L2_bot = 0.0;
    double H1_top = 0.0; double H1_bot = 0.0;

    for (int ee = 0; ee < n_el; ++ee)
    {
        double *x_ele = new double[n_en] ();
        double *u_ele = new double[n_en] ();

        // Global coordinates of x
        for (int aa = 0; aa < n_en; ++aa)
        {
            x_ele[aa] = x_coor[IEN[aa * n_el + ee] - 1];
            u_ele[aa] = disp[IEN[aa * n_el + ee] - 1];
        } 


        for (int ll = 0; ll < n_int; ++ll)
        {
            double dx_dxi = 0.0;
            double x_l = 0.0;
            double du_dxi = 0.0;
            double uh = 0.0;
            // push forward from parent frame Gop->x[ll] to physical frame x_l
            for (int aa = 0; aa < n_en; ++aa)
            {
                x_l = x_l + x_ele[aa] * Polyshape(deg, aa + 1, Gop->x[ll], 0);
                uh  = uh  + u_ele[aa] * Polyshape(deg, aa + 1, Gop->x[ll], 0);
                dx_dxi = dx_dxi + x_ele[aa] * Polyshape(deg, aa + 1, Gop->x[ll], 1);
                du_dxi = du_dxi + u_ele[aa] * Polyshape(deg, aa + 1, Gop->x[ll], 1);
            }

            // inverse of dx_dxi
            double dxi_dx = 1.0 / dx_dxi;

            L2_top = L2_top + Gop->w[ll] * pow( uh - exact(x_l), 2) * dx_dxi;
            L2_bot = L2_bot + Gop->w[ll] * pow(exact(x_l), 2) * dx_dxi;

            H1_top = H1_top + Gop->w[ll] * pow(du_dxi * dxi_dx - exact_x(x_l), 2) * dx_dxi;
            H1_bot = H1_bot + Gop->w[ll] * pow(exact_x(x_l), 2) * dx_dxi;

        }
    }

    L2_top = sqrt(L2_top); L2_bot = sqrt(L2_bot);

    H1_top = sqrt(H1_top); H1_bot = sqrt(H1_bot);

    double L2_error = L2_top / L2_bot;
    double H1_error = H1_top / H1_bot;

    std::cout << "L2_error = " << L2_error << std::endl;
    std::cout << "H1_error = " << H1_error << std::endl;



    delete [] K;
    delete [] F;



    return 0;
}

