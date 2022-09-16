#include <iostream>
#include <string>
#include "heat.h"

using namespace std;

int main()
{
    double nu, T = 0.0014;
    int N = 4;
    vector<double> U_0(N + 2), U(N + 2);

    U_0 = u_0(N);

    U = chaleur_explicite(nu, T, U_0);

    cout << "U= \n"
         << U << endl;

    // ecriture dans un fichier texte les N et epsilon_N
    string nomfichier = "Erreur_heat.txt";
    ofstream file1;

    file1.open(nomfichier, ios::out); // ouverture du fichier
    assert(!file1.fail());            // verification si mon fichier est bien ouvert

    for (int i = 0; i < 8; i++)
    {
        nu = 0.25 / (N + 1);

        verif_chaleur_explicite(nu, T, N, file1); // appel de la fonction verif_poisson
        N += N;
    }
    file1.close();

    return 0;
}
