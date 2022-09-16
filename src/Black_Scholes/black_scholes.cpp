#include <iostream>
#include <string>
#include "black_scholes.h"

using namespace std;

int main()
{
    int N = 25;
    double sigma = 0.2, L = 1000., h = L / N, T = 100.;
    // double sigma = 0.2, L = 10., h = L / N, k , T = 1.;
    vector<double> X(N + 1), U_0(N + 1);

    for (int i = 0; i < N + 1; i++)
    {
        X[i] = i * h;
    }

    // Solution initiale
    U_0 = u0(L, N);

    // cout << "U_0 \n"
    //      << U_0 << endl;

    vector<double> U_sol(N + 1);

    // Résolution
    U_sol = BS_implicite(L, T, sigma, N, U_0);

    return 0;
}