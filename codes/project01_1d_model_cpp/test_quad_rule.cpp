#include "Gauss.hpp"
#include "typeinfo"

int main()
{
    int n_int = 3;

    GaussOutput * Gop = new GaussOutput;

    Gop = Gauss(n_int, -1, 1);

    int n = 6;

    double exact;

    exact = (1 - pow(-1, n+1)) / (n + 1);

    double appro = 0.0;

   

    for (int ll = 0; ll < n_int; ++ll)
    {
        appro = appro + Gop->w[ll] * pow(Gop->x[ll], n);
    }

    double error;

    error = exact - appro;

    std::cout << "error = " << error << std::endl;

    delete Gop;

    return 0;


}