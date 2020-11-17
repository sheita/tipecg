# -*- coding: utf-8 -*-

"""
Lecture du fichier CSV
"""

import csv

Temps = []
Signal = []
  
fichier = csv.reader(open('EX8.csv'), delimiter=";")
for ligne in fichier:
    Temps.append(ligne[0])
    Signal.append(ligne[1])

Temps.pop(0)
for i in range(len(Temps)):
    Temps[i] = Temps[i].replace(',','.')
    Temps[i] = float(Temps[i])

Signal.pop(0)
for i in range(len(Signal)):
    Signal[i] = Signal[i].replace(',','.')
    Signal[i] = float(Signal[i])

print(Temps)
print(Signal)


"""
Affichage du graphe
"""
import matplotlib.pyplot as plt

# Ligne horizontale à y = 0
plt.axhline(y=0, color='black', linestyle='--')

# Tracé de l'électrocardiogramme
plt.title("Electrocardiogramme")
plt.xlabel("Temps (s)")
plt.ylabel("Signal (V)")
plt.plot(Temps,Signal,'r')

# Centrage autour de y = 0
maxabs = max(Signal+[-x for x in Signal])
print(maxabs)
plt.ylim([-maxabs-1,maxabs+1])

# Affichage du graphe en plein écran directement
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()

# Couleur de fond de la fenêtre d'affichage du graphe
plt.rcParams["figure.facecolor"] = 'w'

# Affichage d'une ligne verticale
# plt.axvline(x=1, color='black', linestyle='--')
# plt.axhline(y=1, color='black', linestyle='--')

plt.show()
