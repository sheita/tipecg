# -*- coding: utf-8 -*-

"""
Lecture du fichier CSV
"""

import csv

def lireCSV(nomdufichier):
    Temps, Signal = [], []
    
    """
    Le fichier CSV exportés de Latis Pro sont de la forme :
    Temps;EX8
    0;-0,132368206279352
    0,003;-0,221940886462107
    0,006;-0,147296986309811 (...)
    """
    
    fichier = csv.reader(open(nomdufichier), delimiter=";")
    for ligne in fichier:
        Temps.append(ligne[0])
        Signal.append(ligne[1])
    
    Temps.pop(0) # On enlève le titre de la colonne
    for i in range(len(Temps)):
        Temps[i] = Temps[i].replace(',','.')
        Temps[i] = float(Temps[i])
    
    Signal.pop(0) # On enlève le titre de la colonne
    for i in range(len(Signal)):
        Signal[i] = Signal[i].replace(',','.')
        Signal[i] = float(Signal[i])
    
    return (Temps, Signal)
    
Temps, Signal = lireCSV("EX8.csv")

"""
Affichage du graphe
"""
import matplotlib.pyplot as plt


def tracerECG(Temps, Signal):
    plt.figure("Graphe ECG")
    
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
    
    # Couleur de fond de la fenêtre d'affichage du graphe en blanc
    plt.rcParams["figure.facecolor"] = 'w'
    
    plt.show()

tracerECG(Temps, Signal)

"""
Analyse de fourier
"""

import numpy as np
import scipy.fftpack

def AnalyseFourier(Temps, Signal):
    plt.figure("Analyse de Fourier")

    # L'analyse de Fourier utilise des tableaux numpy
    SignalNp = np.array(Signal)
    
    N = 500
    T = 0.003 # Temps[i]-Temps[i-1] n'est pas consistant
    yf = scipy.fftpack.fft(SignalNp)
    xf = np.linspace(0, 1.0/(2.0*T), N/2)
    
    plt.stem(xf, 2.0/N * np.abs(yf[:N//2]), markerfmt=' ')
    plt.title("Analyse de Fourier du signal ECG")
    plt.xlabel("Fréquence (Hz)")
    plt.ylabel("Amplitude")
    plt.show()
    
    # Couleur de fond de la fenêtre d'affichage du graphe en blanc
    plt.rcParams["figure.facecolor"] = 'w'

# AnalyseFourier(Temps, Signal)


"""
Découpage du signal
"""

# Affichage d'une ligne verticale au maximum d'amplitude de l'ECG
plt.axvline(x=max(Signal), color='black', linestyle='--')


"""
Autres
"""
# Affichage d'une ligne verticale ou horizontale :
# plt.axvline(x=1, color='black', linestyle='--')
# plt.axhline(y=1, color='black', linestyle='--')
