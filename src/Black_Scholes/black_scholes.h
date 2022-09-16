#if !defined(BLACK_SCHOLES_H)
#define BLACK_SCHOLES_H

#include <iostream>
#include <vector>
#include <algorithm> // std::max
#include <fstream>
#include <cmath>
#include "../vecteur_template.h"

using namespace std;

// Résolution de l'equation de Black-Scholes

struct mat_profil // definition de la structure mat_profil
{
    vector<int> INDIAG;

    vector<double> ATAB;
};

struct mat_bande // Strucutre mat_bande cf:question 2
{
    vector<int> IND;
    vector<double> ATAB;
};

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

// les fonctions factorisation_LDLt et resolution_LDLt
void factorisation_LDLt(mat_profil &A)
{
    int k, j, i, ji, jj;
    double S;
    int n = A.INDIAG.size();
    for (i = 0; i < n; i++)
    {
        if (i == 0)
            ji = 0;
        else
            ji = i - A.INDIAG[i] + A.INDIAG[i - 1] + 1;
        for (j = 0; j < i; j++)
        {
            if (j == 0)
                jj = 0;
            else
                jj = j - A.INDIAG[j] + A.INDIAG[j - 1] + 1;
            for (k = max(ji, jj); k < j; k++)

                A.ATAB[A.INDIAG[i] - i + j] -= A.ATAB[A.INDIAG[i] - i + k] * A.ATAB[A.INDIAG[j] - j + k];
        }
        for (j = ji; j < i; j++)
        {
            S = A.ATAB[A.INDIAG[i] - i + j] / A.ATAB[A.INDIAG[j]];
            A.ATAB[A.INDIAG[i]] -= S * A.ATAB[A.INDIAG[i] - i + j];
            A.ATAB[A.INDIAG[i] - i + j] = S;
        }
    }
}

vector<double> resolution_LDLt(mat_profil &A, const vector<double> &b)
{
    vector<double> x(b);
    int n = x.size();
    int ji;
    int i;
    int j;
    double z = 0;
    for (i = 1; i < n; i++)
    {
        ji = i - A.INDIAG[i] + A.INDIAG[i - 1] + 1;

        for (j = ji; j < i; j++)
            x[i] -= A.ATAB[A.INDIAG[i] - i + j] * x[j];
    }

    for (i = 0; i < n; i++)
        x[i] = x[i] / A.ATAB[A.INDIAG[i]];

    for (i = n - 1; i > 0; i--)
    {
        ji = i - A.INDIAG[i] + A.INDIAG[i - 1] + 1;

        for (int j = ji; j < i; j++)
        {
            x[j] -= A.ATAB[A.INDIAG[i] - i + j] * x[i];
        }
    }
    return x;
}

mat_bande factorisation_LU(const mat_bande &A, int n)
{
    int d = A.IND.size();
    mat_bande B;

    // B a les mêmes dimensionnements que A
    B.IND = A.IND;
    B.ATAB = A.ATAB;

    B.ATAB[1] = A.ATAB[1];
    B.ATAB[d] = A.ATAB[d] / A.ATAB[1];

    for (int i = 1; i < n - 1; ++i)
    {
        B.ATAB[(i - 1) * d + 2] = A.ATAB[(i - 1) * d + 2];
        B.ATAB[i * d + 1] = A.ATAB[i * d + 1] - B.ATAB[i * d] * B.ATAB[(i - 1) * d + 2];
        B.ATAB[(i + 1) * d] = A.ATAB[(i + 1) * d] / B.ATAB[i * d + 1];
    }
    B.ATAB[(n - 1) * d + 1] = A.ATAB[(n - 1) * d + 1] - B.ATAB[(n - 1) * d] * B.ATAB[((n - 1) - 1) * d + 2];

    return B;
}

vector<double> resolution_LU(const mat_bande &A, const vector<double> &b)
{
    int n = b.size();
    int d = A.IND.size();
    mat_bande B = factorisation_LU(A, n); // factorisation LU
    vector<double> x(b);                  // Initialisation

    // Descente
    for (int i = 1; i < n; i++)
        x[i] = x[i] - B.ATAB[i * d] * x[i - 1];
    x[n - 1] = x[n - 1] / B.ATAB[(n - 1) * d + 1];
    // Remontée
    for (int i = n - 2; i >= 0; i--)
        x[i] = (x[i] - B.ATAB[i * d + 2] * x[i + 1]) / B.ATAB[i * d + 1];

    return x;
}

// Résolution de l'equation de Black-Scholes

// 1 Définition des fonctions K(x) et u0(L;N)

// fonction K(x)
double K(double x)
{
    return 0.95 * x;
}

// fonction u0(L;N)

vector<double> u0(double L, int N)
{
    double h = L / N, hx; // on a N+1  points et donc d'apres la formule h=b-a/nb de pts-1, on a h=L/N
    vector<double> U0(N + 1);

    for (int i = 0; i < N + 1; i++)
    {
        hx = i * h - K(i * h);
        U0[i] = max(hx, 0.); // discretisation du vecteur X=> X=i*h
    }

    return U0;
}

// 2 Définition la fonction matrice bande

mat_bande matrice_BS_bande(vector<double> X, double sigma, double h, double k)
{
    mat_bande BS;
    int N = X.size() - 1;
    BS.ATAB.resize(3 * (N + 1));
    double c0, c1, c2;

    BS.IND.push_back(-1);
    BS.IND.push_back(0);
    BS.IND.push_back(1);

    BS.ATAB[1] = 1.;

    for (int i = 1; i < N + 1; i++)
    {
        c0 = 1 + 2 * pow(sigma * X[i], 2) * k / (2 * pow(h, 2));
        c1 = -pow(sigma * X[i - 1], 2) * k / (2 * pow(h, 2));
        c2 = -pow(sigma * X[i + 1], 2) * k / (2 * pow(h, 2));
        // Diagonale
        BS.ATAB[i * 3 + 1] = c0;

        // Sous diagonale
        if (i != 1)
            BS.ATAB[i * 3] = c1;

        // Sur diagonale
        if (i != N)
            BS.ATAB[i * 3 + 2] = c2;
    }

    // Sous diagonale
    BS.ATAB[N * 3] = -2 * pow(sigma * X[N - 1], 2) * k / (2 * pow(h, 2));

    return BS;
}

vector<double> BS_impl_n(double L, double sigma, int N, const vector<double> &Tmp)
{
    vector<double> U_sol, X;
    double h = L / N, k = h, x;

    for (int i = 0; i < N + 1; i++)
    {
        x = i * h;
        X.push_back(x);
    }

    mat_bande M = matrice_BS_bande(X, sigma, h, k);
    mat_bande B = factorisation_LU(M, N + 1);
    U_sol = resolution_LU(B, Tmp);

    return U_sol;
}

vector<double> BS_implicite(double L, double T, double sigma, int N, const vector<double> &U_0)
{
    vector<double> U_sol(N + 1), Tmp(N + 1);
    double h = L / N, k = h, t = 0, x;
    ofstream file;
    string nomfichier = "Erreur_black_scholes.txt";

    file.open(nomfichier, ios::out); // ouverture du fichier
    assert(!file.fail());            // verification si mon fichier est bien ouvert

    // Donne N

    file << N << endl;

    U_sol = U_0;
    Tmp = U_sol;

    for (int i = 0; i < N + 1; i++)
    {
        x = i * h;
        file << t << ", " << x << ", " << U_sol[i] << endl;
    }
    file << endl;

    while (t < T + k)
    {
        t += k;

        // condition de Dirichlet à gauche
        Tmp[0] = 0.;

        U_sol = BS_impl_n(L, sigma, N, Tmp); // revoie la solution U_sol à l'atape t

        Tmp = U_sol;

        // ecriture dans un fichier (file) N et epsilon_N
        for (int i = 0; i < N + 1; i++)
        {
            x = i * h;
            file << t << ", " << x << ", " << U_sol[i] << endl;
        }
        file << endl;
    }
    file.close();

    return U_sol;
}

#endif // BLACK_SCHOLES_H
