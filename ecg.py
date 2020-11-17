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


"""
Analyse de spectre
"""

import numpy as np
import scipy.fftpack

SignalNp = np.array(Signal)
TempsNp = np.array(Temps)

# Number of samplepoints
N = 600
# sample spacing
T = 0.003
x = TempsNp
y = SignalNp
yf = scipy.fftpack.fft(y)
xf = np.linspace(0.0, 1.0/(2.0*T), N/2)

fig, ax = plt.subplots()
ax.stem(xf, 2.0/N * np.abs(yf[:N//2]), markerfmt=' ')
#ax.title("Analyse de Fourier du signal ECG")
#fig.xlabel("Fréquence (Hz)")
#fig.ylabel("Amplitude")
plt.show()
plt.figure()
plt.plot([1],[1])
plt.show()
