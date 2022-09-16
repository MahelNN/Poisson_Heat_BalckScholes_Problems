# Reading the file as one file stream
import csv
from os import read
import numpy as np
import matplotlib.pyplot as plt

# Equation de la chaleur
data_x = []
data_y = []

f = open('Erreur_heat.txt', 'r')
for row in f:
    row = row.split()
    data_x.append(float(row[0]))
    data_y.append(float(row[1]))

plt.plot(np.log(data_x), np.log(data_y), color='g',
         marker="d", label='$\epsilon_N$')
plt.plot(np.log(data_x), 0.5-1*np.log(data_x),
         color='gray', linestyle='-.', label="$y = -x+1/2 $")

plt.xlabel('Nombre de points (log)', fontsize=12)
plt.ylabel('log(erreur)', fontsize=12)

plt.title('Erreur en fonction du nombre de points', fontsize=20)

# plt.xscale('log')
# plt.yscale('log')

plt.legend()
plt.grid()
plt.show()

print(data_x)
print(data_y)
