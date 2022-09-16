# Reading the file as one file stream
import csv
from os import read
import numpy as np
import matplotlib.pyplot as plt

# Black Scholes problem
data_x = []
data_y = []
data_t = []

# Extraction des données

f = open('Erreur_black_scholes.txt', 'r')

Ni = int(f.readline())

d = np.genfromtxt(f, delimiter=",", autostrip=True)

Nf = np.size(d[:, 2])

# print("d\n", d)
print("Nf\t", Nf)

for i in range(Nf):
    if(i % (Ni+1) == 0):
        data_t.append(d[i][0])
    if(i < Ni+1):
        data_x.append(d[i][1])
    data_y.append(d[i][2])

Nj = np.size(data_t)-1

# print("Ni\t", Ni)
# print("Nj\t", Nj)


M_sol = np.array(data_y).reshape(Ni+1, Nj+1)

# print("\nx = \t", data_x)
# print("\nt = \t", data_t)
# print("\n", M_sol)


# Evolution de la solution au cours du temps
X = np.array(data_x)

T = np.array(data_t)

T, X = np.meshgrid(T, X)

# print(X, "\n\n ", T)

fig = plt.figure()

ax = fig.add_subplot(111, projection='3d')

ax.plot_surface(T, X, M_sol, cmap="plasma")

ax.set_xlabel(" t ")

ax.set_ylabel(" x ")

ax.set_zlabel("U_sol ")


plt.show()


# Evolution du stock au temps final (regarder ici pour le payoff et le prix de l'achat)

label = range(Nj+1)

plt.plot(data_x, M_sol[:, -1], label="t"+str(label[-1]))

plt.xlabel(' x ', fontsize=12)
plt.ylabel(' Evolution du stock', fontsize=12)

plt.title('Evolution du stock temps final', fontsize=20)


plt.legend()
plt.grid()
plt.show()
