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
            if (der == 0) val =0.5 * (1 - xi);
            else if (der == 1) val = -0.5;
        }
        else if (a == 2)
        {
            if (der == 0) val = 0.5 * (1 + xi);
            else if (der == 1) val = 0.5;
        }
        
        break;
    // Quadratic basis function
    case 2:
        if (a == 1)
        {
            if (der == 0) val = 0.5 * xi * (xi - 1);
            else if (der == 1) val = xi - 0.5;
        }
        else if (a == 2)
        {
            if (der == 0) val = 1 - xi * xi;
            else if (der == 1) val = -2 * xi;
        }
        else if (a == 3)
        {
            if (der == 0) val = 0.5 * xi * (xi + 1);
            else if (der == 1) val  = xi + 0.5; 
        }

        break;
    // Cubic basis function
    case 3:
        if (a == 1)
        {
            if (der == 0) val = -9 * (xi - (1/3)) * (xi + (1/3)) * (xi -1) / 16;
            else if (der == 1) val = -9 * (2 * xi * (xi - 1)+ xi * xi - (1/9)) / 16;
        }
        else if (a == 2)
        {
            if (der == 0) val = 27 * (xi * xi - 1) * (xi - (1/3)) / 16;
            else if (der == 1) val = 27 * (2 * xi * (xi - (1/3)) + xi * xi - 1) / 16;
        }
        else if (a == 3)
        {
            if (der == 0) val = -27 * (xi * xi - 1) * (xi + (1/3)) / 16;
            else if (der == 1) val = -27 * (2 * xi * (xi + (1/3)) + xi * xi - 1) / 16; 
        }
        else if (a == 4)    
        {           
            if (der == 0) val  = 9 * (xi + 1)*(xi * xi - (1/9)) / 16;
            else if (der == 1) val = 9 * (xi * xi - (1/9) + 2 * xi * (xi + 1)) / 16;
        }
        
        break;

    
    default:
        fprintf(stderr, "Error: degree has to be 1, 2, or 3.");
        break;
    }

    return val;
}
