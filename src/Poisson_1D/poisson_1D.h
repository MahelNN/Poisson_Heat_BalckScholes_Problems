#if !defined(POISSON_1D_H)
#define POISSON_1D_H
#include <iostream>
#include <vector>
#include "../vecteur_template.h"

using namespace std;

struct mat_profil // definition de la structure mat_profil
{
    vector<int> INDIAG;
    vector<double> ATAB;
};

// 1) definition de la fonction matrice_1D_profil
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
        A.INDIAG[i] = 2 * i; // stokage des indices des elements diagonaux de A
    }

    return A;
}

// 2)
// les fonctions factorisation_LDLt et resolution_LDLt (deja definis dans le TP3)
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

// definition de fonctions f et g
double f(double x)
{
    double f = cos(x) - 2;
    return f;
}
double g(double x)
{
    double g = pow(x, 2) + cos(x);
    return g;
}

// definition de la fonction verif_poisson
void verif_poisson(double a, double b, int N, ofstream &file)
{
    mat_profil A;
    double h = abs((b - a)) / (N + 1); // ici on met la valeur absolue pour s'assurer de la positivité du pas h

    vector<double> vect(N, 0); //  pour contenir   le second membre de notre systeme d'equation
    vector<double> uex(N, 0);  // pour contenir   le vecteur  u exact.
    double x;
    // remplissage du vecteur vect par les f_i et f_1 + g(a) et f_n + g(b) et celui Uexact

    x = a + h;
    vect[0] = f(x) + g(a) / pow(h, 2);

    for (int i = 1; i < N - 1; i++)
    {
        int j = i + 1;
        x = a + j * h;
        vect[i] = f(x);
    }

    // remplissage du vecteur u exact
    for (int i = 0; i < N; i++)
    {
        double y = a + (i + 1) * h;
        uex[i] = g(y);
    }

    vect[N - 1] = f(a + N * h) + g(b) / pow(h, 2);
    uex[N - 1] = g(a + N * h);
    double c0 = -1. / pow(h, 2);
    double c1 = 2. / pow(h, 2);

    A = matrice_1D_profil(N, c0, c1); // pour recuperer la matrice en question en stockage profil.

    factorisation_LDLt(A); // factorisation_LDLt de A

    // cout<<"le vecteur ATAB apres factorisation"<<A.ATAB<<endl;

    vector<double> u_sol(N, 0);

    u_sol = resolution_LDLt(A, vect); // pour contenir la solution de notre systeme d'equation
                                      // Des instructions pour verifier si notre programmme fonctionne
    //  cout<<"le vecteur second membre:"<<endl;
    //  cout<<vect<<endl;
    //  cout<<"le vecetur solution:"<<endl;
    //  cout<<u_sol<<endl;
    //  cout<<"le vecteur u_exact:"<<endl;
    //  cout<<uex<<endl;

    // calcul de l'erreur epsilon_N
    vector<double> w = u_sol - uex;
    double epsilon_N = norme_inf(w);
    // ecriture a l'ecran N et epsilon_N
    cout << "pour le nombre de point internes du maillage, N= " << N << endl;
    cout << "l'erreur epsilon_N est egale a " << epsilon_N << endl;
    // ecriture dans un fichier (file) N et epsilon_N
    file << N << "    " << epsilon_N << endl;
}

#endif // POISSON_1D_H
