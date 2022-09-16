#if !defined(HEAT_H)
#define HEAT_H
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include "../vecteur_template.h"

using namespace std;

struct mat_profil // definition de la structure mat_profil
{
    vector<int> INDIAG;

    vector<double> ATAB;
};
// creation des fonctions uex
double uex(double t, double x)
{
    return sin(t) * cos(M_PI * x);
}
// la fonction derivee de uex par rapport a t
double uex_t(double t, double x)
{
    return cos(t) * cos(M_PI * x);
}

// la fonction derivee seconde de uex par rapport a x
double uex_xx(double t, double x)
{
    return -(M_PI * M_PI) * uex(t, x);
}

double f(double t, double x, double nu)
{
    double fg = uex_t(t, x) - nu * uex_xx(t, x);
    return fg;
}

// Solution discrète initiale

vector<double> u_0(int N)
{
    double h = 1. / (N + 1);
    vector<double> U_0(N + 2);

    for (int i = 0; i < N + 2; i++)
    {
        U_0[i] = uex(0, i * h);
    }
    return U_0;
}

// definition de la fonction matrice_1D_profil

mat_profil matrice_1D_profil(int N, double c0, double c1)
{
    mat_profil A;

    // On donne d'abord des tailles aux vecteurs ATAB et INDIAG

    A.INDIAG.resize(N);

    A.ATAB.resize(2 * N - 1);

    // Ainsi comme notre matrice profile A s'illustrera par les vecteurs ATAB et INDIAG; on remplit ces derniers comme suit:

    for (int j = 0; j < 2 * N - 1; j++)
    {
        if (j % 2 == 0)

            A.ATAB[j] = c1; // stockage des elements de la diagonale de A

        else

            A.ATAB[j] = c0; // stockage des elements de la sous diagonale de A
    }

    cout << endl;

    for (int i = 0; i < N; i++)
    {
        A.INDIAG[i] = 2 * i; // stokage  des indices des elements diagonaux de A
    }

    return A;
}

vector<double> chaleur_explicite(double nu, double T, const vector<double> &U_0)
{

    int N = U_0.size() - 2;

    double h = 1. / (N + 1);

    double k = h * h / nu;

    int K = ceil(T / k) + 1; // ceil pour avoir la partie entiere superieure du quotient

    vector<double> Tmp(U_0), U_sol(N + 2);

    // cout << "K =\n"
    //      << K << endl;

    U_sol = U_0;
    Tmp = U_sol;

    for (int n = 1; n < K; n++)
    {
        U_sol[0] = uex(n * k, 0);
        // remplissage du vecteur U et du vecteur F_t au temps t
        for (int i = 1; i < N + 1; i++)
        {
            U_sol[i] = Tmp[i] + nu * k / (h * h) * (Tmp[i + 1] - 2 * Tmp[i] + Tmp[i - 1]) + k * f(n * k, i * h, nu);
        }
        U_sol[N + 1] = uex(n * k, (N + 1) * h);
        Tmp = U_sol; //
        // cout << "U_n =\n"
        //      << U_sol << endl;
    }
    return U_sol;
}

void verif_chaleur_explicite(double nu, double T, int N, ofstream &file)
{
    double h = 1. / (N + 1);
    vector<double> U_0(N + 2), U_sol(U_0), U_ex(U_0);
    U_0 = u_0(N);
    U_sol = chaleur_explicite(nu, T, U_0);

    double x;
    for (int i = 0; i < N + 2; i++)
    {
        x = i * h;
        U_ex[i] = uex(T, x);
    }

    vector<double> Vect_err = U_sol - U_ex;
    double epsilon = norme_inf(U_sol - U_ex);

    // ecriture a l'ecran N et epsilon_N
    cout << "pour le nombre de point internes du maillage, N = " << N << endl;
    cout << "l'erreur epsilon_N est egale a " << epsilon << endl;
    // ecriture dans un fichier (file) N et epsilon_N
    file << N << "  " << epsilon << endl;
}

#endif // HEAT_H
