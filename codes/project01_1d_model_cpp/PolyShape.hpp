#pragma once
#include "iostream"

double Polyshape(int degree, int a, double xi, int der)
{
    double val=0.0;
    switch (degree)
    {
    // Linear basis function
    case 1:
        if (a == 1)
        {
            if (der == 0) val =0.5 * (1.0 - xi);
            else if (der == 1) val = -0.5;
        }
        else if (a == 2)
        {
            if (der == 0) val = 0.5 * (1.0 + xi);
            else if (der == 1) val = 0.5;
        }
        
        break;
    // Quadratic basis function
    case 2:
        if (a == 1)
        {
            if (der == 0) val = 0.5 * xi * (xi - 1.0);
            else if (der == 1) val = xi - 0.5;
        }
        else if (a == 2)
        {
            if (der == 0) val = 1.0 - xi * xi;
            else if (der == 1) val = -2.0 * xi;
        }
        else if (a == 3)
        {
            if (der == 0) val = 0.5 * xi * (xi + 1.0);
            else if (der == 1) val  = xi + 0.5; 
        }

        break;
    // Cubic basis function
    case 3:
        if (a == 1)
        {
            if (der == 0) val = -9.0 * (xi - (1.0/3.0)) * (xi + (1.0/3.0)) * (xi -1.0) / 16.0;
            else if (der == 1) val = -9.0 * (2.0 * xi * (xi - 1.0)+ xi * xi - (1.0/9.0)) / 16.0;
        }
        else if (a == 2)
        {
            if (der == 0) val = 27.0 * (xi * xi - 1.0) * (xi - (1.0/3.0)) / 16.0;
            else if (der == 1) val = 27.0 * (2 * xi * (xi - (1.0/3.0)) + xi * xi - 1.0) / 16.0;
        }
        else if (a == 3)
        {
            if (der == 0) val = -27.0 * (xi * xi - 1.0) * (xi + (1.0/3.0)) / 16.0;
            else if (der == 1) val = -27.0 * (2.0 * xi * (xi + (1.0/3.0)) + xi * xi - 1.0) / 16.0; 
        }
        else if (a == 4)    
        {           
            if (der == 0) val  = 9.0 * (xi + 1.0)*(xi * xi - (1.0/9.0)) / 16.0;
            else if (der == 1) val = 9.0 * (xi * xi - (1.0/9.0) + 2.0 * xi * (xi + 1.0)) / 16.0;
        }
        
        break;

    
    default:
        fprintf(stderr, "Error: degree has to be 1, 2, or 3.");
        break;
    }

    return val;
}
