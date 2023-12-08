#include "Gauss2D.hpp"
#include "Quad.hpp"
#include "Quad_grad.hpp"
#include "LU.hpp"
#include <cmath>





//--------------------------------------------------------
// Problem definition
// exact solution

double kappa = 1.0;   // isotropic homogeneous heat conductivity
long double PI = 3.14159265358979323846264338327950288419716939937510582;

// manufactured solution and source term
double exact(double x, double y)
{
    return x * (1.0 - x) * y * (1.0 - y);
}

double exact_x(double x, double y)
{
    return (1.0 - 2.0 * x) * y * (1.0 - y);
}

double exact_y(double x, double y)
{
    return x * (1.0 - x) * (1.0 - 2.0 * y);
}

double f(double x, double y)
{
    return -2.0 * x * (x - 1.0) - 2.0 * y * (y - 1.0);
}

// Dirichlet BC
double g(double x, double y)
{
    return 0.1 * sin((x + y) * 2 * PI);
}

//--------------------------------------------------------

int main()
{
    //--------------------------------------------------------
    // Quadrature Rule
    int n_int_xi = 3;
    int n_int_eta = 3;
    int n_int = n_int_xi * n_int_eta;
    Gauss2D_Output * G2DOp = new Gauss2D_Output;
    G2DOp = Gauss2D(n_int_xi, n_int_eta);
    // for (int ii = 0; ii < n_int_xi * n_int_eta; ++ii)
    // {
    //     std::cout << "xi[" << ii << "] = " << G2DOp->xi[ii] << std::endl;
    // }
    // for (int ii = 0; ii < n_int_xi * n_int_eta; ++ii)
    // {
    //     std::cout << "eta[" << ii << "] = " << G2DOp->eta[ii] << std::endl;
    // }
    // for (int ii = 0; ii < n_int_xi * n_int_eta; ++ii)
    // {
    //     std::cout << "w[" << ii << "] = " << G2DOp->w[ii] << std::endl;
    // }
    
    //--------------------------------------------------------
    // Parameters of the FEM
    int n_en = 4;                // 4-node quadrilateral element

    int n_el_x = 100;            // number of element in x-direction
    int n_el_y = 100;            // number of element in y-direction
    int n_el = n_el_x * n_el_y;  // total number of element in 2D domain
    
    int n_np_x = n_el_x + 1;     // number of node points in x-direction
    int n_np_y = n_el_y + 1;     // number of node points in y-direction 
    int n_np = n_np_x * n_np_y;  // total number of node points in 2D domain   

    int n_eq = 0;                // number of equation

    //--------------------------------------------------------
    // Generate the mesh: coordinates of the nodal points
    // nodal coordinates
    double *x_coor = new double[n_np]();
    double *y_coor = new double[n_np]();

    double hh_x = 1.0 / n_el_x;         // mesh size in the x-direction
    double hh_y = 1.0 / n_el_y;         // mesh size in the y-direction

    for (int ny = 0; ny < n_np_y; ++ny)
    {
        for (int nx = 0; nx < n_np_x; ++nx)
        {
            int index = ny * n_np_x + nx;   // nodal index for the (nx, ny) node
            x_coor[index] = nx * hh_x;
            y_coor[index] = ny * hh_y;
        }
    }
    // for (int ii = 0; ii < n_np; ++ii)
    // {
    //     std::cout << "x_coor[" << ii << "] = " << x_coor[ii] << std::endl;
    // }
    // for (int ii = 0; ii < n_np; ++ii)
    // {
    //     std::cout << "y_coor[" << ii << "] = " << y_coor[ii] << std::endl;
    // }

    //--------------------------------------------------------

    // setup the IEN array for element with local node numbering as
    // a=4 ------- a=3
    // |           |
    // |           |
    // |           |
    // |           |
    // a=1 ------- a=2
    int IEN[n_en * n_el];
    for (int ey = 0; ey < n_el_y; ++ey)
    {
        for (int ex = 0; ex < n_el_x; ++ex)
        {
            int ee = ey * n_el_x + ex;
            IEN[0 * n_el + ee] = ey * n_np_x + ex + 1;
            IEN[1 * n_el + ee] = ey * n_np_x + ex + 2;
            IEN[2 * n_el + ee] = (ey + 1) * n_np_x + ex + 2;
            IEN[3 * n_el + ee] = (ey + 1) * n_np_x + ex + 1;
        }
    }
    // for (int ii = 0; ii < n_en * n_el; ++ii)
    // {
    //     std::cout << "IEN[" << ii << "] = " << IEN[ii] << std::endl;
    // }

    //--------------------------------------------------------

    // ID and LM arrays are generated based on the BC info
    int *ID = new int[n_np]();              // Space dimension is 2D
    int counter = 1;                        // degree of freedom is 1 (heat problem, only temperature)
    // Dirichlet BC is around the domain
    for (int ny = 1; ny < n_np_y - 1; ++ny)
    {
        for (int nx = 1; nx < n_np_x - 1; ++nx)
        {
            ID[ny * n_np_x + nx] = counter;
            counter = counter + 1;
        }
    }
    // for (int ii = 0; ii < n_np; ++ii)
    // {
    //     std::cout << "ID[" << ii << "] = " << ID[ii] << std::endl;
    // }

    // According the Dirichlet BC define the number of equation
    n_eq = n_np - 2 * n_np_x - 2 * n_np_y + 4;

    double LM[1 * n_en * n_el];            // 2D, heat conduction problem (only temp degree of freedom)
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

    //--------------------------------------------------------
    // Start the assembly procedure
    // Assembly of K and F
    double *K = new double[n_eq * n_eq] ();   // allocate the global stiffness matrix
    double *F = new double[n_eq] ();          // allocate the global load vector

    //--------------------------------------------------------
    // Loop of element number
    for (int ee =0; ee < n_el; ++ee)
    {
        double *k_ele = new double[n_en * n_en] ();
        double *f_ele = new double[n_en] ();
        double *x_ele = new double[n_en] ();
        double *y_ele = new double[n_en] ();
        
        // --------------------------------------------------------
        // Global coordinates of x and y of nodes in ee(th) element
        for (int aa = 0; aa < n_en; ++aa)
        {
            x_ele[aa] = x_coor[IEN[aa * n_el + ee] - 1];
            y_ele[aa] = y_coor[IEN[aa * n_el + ee] - 1];
        } 
        for (int ii = 0; ii < n_en; ++ii)
        {
            // std::cout << "x_ele[" << ii << "] = " << x_ele[ii] << std::endl;
            // std::cout << "y_ele[" << ii << "] = " << y_ele[ii] << std::endl;
        }

        //--------------------------------------------------------
        // Compute every Gauss points, loop over quadrature points
        for (int ll = 0; ll < n_int; ++ll)
        {
            double x_l = 0.0; double y_l = 0.0;
            double dx_dxi = 0.0; double dx_deta = 0.0;
            double dy_dxi = 0.0; double dy_deta = 0.0;
            double Na_xi = 0.0; double Na_eta = 0.0;
            double Nb_xi = 0.0; double Nb_eta = 0.0;
            

            //--------------------------------------------------------
            // push forward from parent frame x_ele * Quad (x_e * Na(ll)) to physical frame x_l
            // push forward from parent frame y_ele * Quad (y_e * Na(ll)) to physical frame y_l
            // push forward from parent frame x_ele * QdGdOp->val_xi[ll] (x_e * Na,xi(ll)) to physical frame dx_dxi
            // push forward from parent frame x_ele * QdGdOp->val_eta[ll] (x_e * Na,eta(ll)) to physical frame dx_deta
            // push forward from parent frame y_ele * QdGdOp->val_xi[ll] (y_e * Na,xi(ll)) to physical frame dy_dxi
            // push forward from parent frame y_ele * QdGdOp->val_eta[ll] (y_e * Na,eta(ll)) to physical frame dy_deta
            for (int aa = 0; aa < n_en; ++aa)
            {
                Quad_grad_Oput * QdGdOp = new Quad_grad_Oput; 

                QdGdOp = Quad_grad(aa + 1, G2DOp->xi[ll], G2DOp->eta[ll]);
                x_l = x_l + x_ele[aa] * Quad(aa + 1, G2DOp->xi[ll], G2DOp->eta[ll]);
                y_l = y_l + y_ele[aa] * Quad(aa + 1, G2DOp->xi[ll], G2DOp->eta[ll]);
                dx_dxi = dx_dxi + x_ele[aa] * QdGdOp->val_xi;
                dx_deta = dx_deta + x_ele[aa] * QdGdOp->val_eta;
                dy_dxi = dy_dxi + y_ele[aa] * QdGdOp->val_xi;
                dy_deta = dy_deta + y_ele[aa] * QdGdOp->val_eta;

                delete QdGdOp;
                
            }
            // std::cout << "x_l(" << ll << ") = " << x_l << std::endl;
            // std::cout << "y_l(" << ll << ") = " << y_l << std::endl;
            // std::cout << "dx_dxi(" << ll << ") = " << dx_dxi << std::endl;
            // std::cout << "dx_deta(" << ll << ") = " << dx_deta << std::endl;
            // std::cout << "dy_dxi(" << ll << ") = " << dy_dxi << std::endl;
            // std::cout << "dy_deta(" << ll << ") = " << dy_deta << std::endl;
            
            // Jacobian determinant
            double detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;

            

            // --------------------------------------------------------
            // generate element stiffness matrix and element force vector
            for (int aa = 0; aa < n_en; ++aa)
            {
                f_ele[aa] = f_ele[aa] + G2DOp->w[ll] * f(x_l, y_l) * Quad(aa + 1, G2DOp->xi[ll], G2DOp->eta[ll]) * detJ; 

                Quad_grad_Oput * QdGdOp_aa = new Quad_grad_Oput; 
                QdGdOp_aa = Quad_grad(aa + 1, G2DOp->xi[ll], G2DOp->eta[ll]);
                Na_xi = QdGdOp_aa->val_xi; Na_eta = QdGdOp_aa->val_eta;
                double Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
                double Na_y = (Na_xi * (-dx_deta) + Na_eta * dx_dxi) / detJ;
                // std::cout << "Na_x = " << Na_x << std::endl;
                // std::cout << "Na_y = " << Na_y << std::endl;
                for (int bb = 0; bb < n_en; ++bb)
                {
                    Quad_grad_Oput * QdGdOp_bb = new Quad_grad_Oput; 
                    QdGdOp_bb = Quad_grad(bb + 1, G2DOp->xi[ll], G2DOp->eta[ll]);
                    Nb_xi = QdGdOp_bb->val_xi; Nb_eta = QdGdOp_bb->val_eta;
                    double Nb_x = (Nb_xi * dy_deta -  Nb_eta * dy_dxi) / detJ;
                    double Nb_y = (Nb_xi * (-dx_deta) + Nb_eta * dx_dxi) / detJ;
                    // std::cout << "Nb_x = " << Nb_x << std::endl;
                    // std::cout << "Nb_y = " << Nb_y << std::endl;
                    k_ele[aa * n_en + bb] = k_ele[aa * n_en + bb] + G2DOp->w[ll] * kappa * (Na_x * Nb_x + Na_y * Nb_y) * detJ;

                    delete QdGdOp_bb;
                }

                delete QdGdOp_aa;
            }
            // for (int ii = 0; ii < n_en * n_en; ++ii)
            // {
            //     std::cout << "k_ele [" << ii << "] = " << k_ele[ii] << std::endl;
            // }
            // for (int ii = 0; ii < n_en; ++ii)
            // {
            //     std::cout << "f_ele [" << ii << "] = " << f_ele[ii] << std::endl;
            // }
        }

        //--------------------------------------------------------
        // Now we need to put element k and f into global K and F, global assembly
        for (int aa = 0; aa < n_en; ++aa)
        {
            int PP;
            PP = LM[aa * n_el + ee];
            if (PP > 0)
            {
                F[PP - 1] = F[PP - 1] + f_ele[aa];
                for (int bb = 0; bb < n_en; ++bb)
                {
                    int QQ;
                    QQ = LM[bb * n_el + ee];
                    if (QQ > 0)
                    {
                        K[(PP - 1) * n_eq + QQ - 1] = K[(PP - 1) * n_eq + QQ - 1] + k_ele[aa * n_en + bb];
                    } 
                    else 
                    {
                        // F[PP - 1] = F[PP - 1] - k_ele[aa * n_en + bb] * g;
                        // do something for non-zero g boundary condition
                        F[PP - 1] = F[PP - 1] - k_ele[aa * n_en + bb] * g(x_ele[bb], y_ele[bb]);
                        
                    }
                }
            }

        }

        // do something for h boundary condition
        // if (ee == 0)
        // {
        //     F[ID[IEN[ 0 * n_el + ee] - 1] - 1] = F[ID[IEN[ 0 * n_el + ee] - 1] - 1] + h;
        // }



        delete [] k_ele;
        delete [] f_ele;
        delete [] x_ele;
        delete [] y_ele;
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
    double * d_temp = new double[n_eq] ();
    d_temp = LU_solve(F, n_eq, lu_var->pp, lu_var->mat);
    // for (int ii = 0; ii < n_eq; ++ii)
    // {
    //     std::cout << "d_temp[" << ii << "] = " << d_temp[ii] << std::endl;
    // }

    // --------------------------------------------------------
    // Generate the full solution vector by inserting back the Dirichlet value
    double *disp = new double[n_np]();
    for (int ii = 0; ii < n_np; ++ii)
    {
        int index = ID[ii];
        if (index > 0) disp[ii] = d_temp[index - 1];
    }
    // disp[n_np - 1] = g;   
    // non-zero Dirichlet BC
    // for (int ii = 0; ii < n_np; ++ii)
    // {
    //     std::cout << "disp[" << ii << "] = " << disp[ii] << std::endl;
    // }

    //--------------------------------------------------------
    // Plot the solution


    //-------------------------------------------------------
    // Calculate the error
    // Poetprocess the solution by calculating the error measured in L2 norm and H1 norm
    // errors
    double L2_top = 0.0; double L2_bot = 0.0;
    double H1_top = 0.0; double H1_bot = 0.0;

    //-------------------------------------------------------
    // Loop of element number
    for (int ee = 0; ee < n_el; ++ee)
    {
        double *x_ele = new double[n_en] ();
        double *y_ele = new double[n_en] ();
        double *u_ele = new double[n_en] ();

        //-------------------------------------------------------
        //Global coordinates of x and y of nodes in ee(th) element
        //Global displacement of nodes in ee(th) element
        for (int aa = 0; aa < n_en; ++aa)
        {
            x_ele[aa] = x_coor[IEN[aa * n_el + ee] - 1];
            y_ele[aa] = y_coor[IEN[aa * n_el + ee] - 1];
            u_ele[aa] = disp[IEN[aa * n_el + ee] - 1];
        } 
        // for (int ii = 0; ii < n_en; ++ii)
        // {
        //     std::cout << "x_ele[" << ii << "] = " << x_ele[ii] << std::endl;
        //     std::cout << "y_ele[" << ii << "] = " << y_ele[ii] << std::endl;
        //     std::cout << "u_ele[" << ii << "] = " << y_ele[ii] << std::endl;
        // }


        //--------------------------------------------------------
        // Compute every Gauss points, loop over quadrature points
        for (int ll = 0; ll < n_int; ++ll)
        {
            double x_l = 0.0; double y_l = 0.0; double u_l = 0.0;
            double dx_dxi = 0.0; double dx_deta = 0.0; 
            double dy_dxi = 0.0; double dy_deta = 0.0; 
            double du_dxi = 0.0; double du_deta = 0.0;

            //--------------------------------------------------------
            // push forward from parent frame x_ele * Quad (x_e * Na(ll)) to physical frame x_l
            // push forward from parent frame y_ele * Quad (y_e * Na(ll)) to physical frame y_l
            // push forward from parent frame u_ele * Quad (u_e * Na(ll)) to physical frame u_l
            // push forward from parent frame x_ele * QdGdOp->val_xi[ll] (x_e * Na,xi(ll)) to physical frame dx_dxi
            // push forward from parent frame x_ele * QdGdOp->val_eta[ll] (x_e * Na,eta(ll)) to physical frame dx_deta
            // push forward from parent frame y_ele * QdGdOp->val_xi[ll] (y_e * Na,xi(ll)) to physical frame dy_dxi
            // push forward from parent frame y_ele * QdGdOp->val_eta[ll] (y_e * Na,eta(ll)) to physical frame dy_deta
            for (int aa = 0; aa < n_en; ++aa)
            {
                Quad_grad_Oput * QdGdOp = new Quad_grad_Oput; 

                QdGdOp = Quad_grad(aa + 1, G2DOp->xi[ll], G2DOp->eta[ll]);
                x_l = x_l + x_ele[aa] * Quad(aa + 1, G2DOp->xi[ll], G2DOp->eta[ll]);
                y_l = y_l + y_ele[aa] * Quad(aa + 1, G2DOp->xi[ll], G2DOp->eta[ll]);
                u_l = u_l + u_ele[aa] * Quad(aa + 1, G2DOp->xi[ll], G2DOp->eta[ll]);
                dx_dxi = dx_dxi + x_ele[aa] * QdGdOp->val_xi;
                dx_deta = dx_deta + x_ele[aa] * QdGdOp->val_eta;
                dy_dxi = dy_dxi + y_ele[aa] * QdGdOp->val_xi;
                dy_deta = dy_deta + y_ele[aa] * QdGdOp->val_eta;
                du_dxi = du_dxi + u_ele[aa] * QdGdOp->val_xi;
                du_deta = du_deta + u_ele[aa] * QdGdOp->val_eta;

                delete QdGdOp;
                
            }
            // std::cout << "x_l(" << ll << ") = " << x_l << std::endl;
            // std::cout << "y_l(" << ll << ") = " << y_l << std::endl;
            // std::cout << "dx_dxi(" << ll << ") = " << dx_dxi << std::endl;
            // std::cout << "dx_deta(" << ll << ") = " << dx_deta << std::endl;
            // std::cout << "dy_dxi(" << ll << ") = " << dy_dxi << std::endl;
            // std::cout << "dy_deta(" << ll << ") = " << dy_deta << std::endl;
            
            // Jacobian determinant
            double detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;
            

            //--------------------------------------------------------
            // generate the error L2 norm and H1 norm
            double du_dx = (du_dxi * dy_deta - du_deta * dy_dxi) / detJ;
            double du_dy = (du_dxi * (-dx_deta) + du_deta * dx_dxi) / detJ;
            // std::cout << "du_dx(" << ll << ") = " << du_dx << std::endl;
            // std::cout << "du_dy(" << ll << ") = " << du_dy << std::endl;

            L2_top = L2_top + G2DOp->w[ll] * pow( u_l - exact(x_l, y_l), 2) * detJ;
            L2_bot = L2_bot + G2DOp->w[ll] * pow(exact(x_l, y_l), 2) * detJ;

            H1_top = H1_top + G2DOp->w[ll] * (pow(du_dx - exact_x(x_l, y_l), 2) + pow(du_dy - exact_y(x_l, y_l), 2)) * detJ;
            H1_bot = H1_bot + G2DOp->w[ll] * (pow(exact_x(x_l, y_l), 2), pow(exact_y(x_l, y_l), 2)) * detJ;

        }

        delete [] x_ele;
        delete [] y_ele;
        delete [] u_ele;
    }

    L2_top = sqrt(L2_top); L2_bot = sqrt(L2_bot);

    H1_top = sqrt(H1_top); H1_bot = sqrt(H1_bot);

    double L2_error = L2_top / L2_bot;
    double H1_error = H1_top / H1_bot;

    std::cout << "L2_error = " << L2_error << std::endl;
    std::cout << "H1_error = " << H1_error << std::endl;



    delete [] K;
    delete [] F;
    delete [] ID;

    delete G2DOp;
    delete [] x_coor;
    delete [] y_coor;

    delete [] pp;
    delete lu_var;
    delete [] d_temp;
    delete [] disp;

    return 0;
}
