#include <iostream>
#include <string>
#include "poisson_1D.h"

using namespace std;

// 3)
int main()
{
    double a = 4, b = 10;
    int N = 8;
    ofstream fichier;
    double h = (b - a) / (N + 1);
    double c0 = -1. / pow(h, 2);
    double c1 = 2. / pow(h, 2);
    //  mat_profil B;
    //  B=matrice_1D_profil(N, c0, c1);
    //  cout<<"le vecteur ATAB est : "<<endl;
    //  cout<<B.ATAB<<endl;

    //  cout<<"le vecteur INDIAG est : "<<endl;
    //  cout<<B.INDIAG<<endl;

    // ecriture dans un fichier texte les N et epsilon_N
    string nomfichier = "Erreur_poisson_1D.txt";
    ofstream file;
    file.open(nomfichier, ios::out); // ouverture du fichier
    assert(!file.fail());            // verification si mon fichier est bien ouvert
    for (int i = 0; i < 9; i++)
    {
        verif_poisson(a, b, N, file); // appel de la fonction verif_poisson
        N += N;
    }
    file.close();

    return 0;
}
