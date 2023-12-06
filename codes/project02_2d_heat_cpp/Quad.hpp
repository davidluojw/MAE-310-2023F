#pragma once
#include "iostream"


double Quad(int aa, double xi, double eta)
{
    double val = 0.0;
    if (aa == 1) val = 0.25 * (1.0 - xi) * (1.0 - eta);
    else if (aa == 2) val = 0.25 * (1.0 + xi) * (1.0 - eta);
    else if (aa == 3) val = 0.25 * (1.0 + xi) * (1.0 + eta);
    else if (aa == 4) val = 0.25 * (1.0 - xi) * (1.0 + eta);
    else fprintf(stderr, "Error: degree has to be 1, 2, or 3.");

    return val;
}