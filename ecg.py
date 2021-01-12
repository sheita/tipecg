# -*- coding: utf-8 -*-

"""
Configuration (pas important)
"""

import matplotlib.pyplot as plt
print("\n")

# Couleur de fond de la fenêtre d'affichage du graphe en blanc
plt.rcParams["figure.facecolor"] = 'w'

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
    
    print("Fichier lu : %s \nNombre de points : %s \nPériode d'échantillonage : %s s \nTemps d'acquisition total : %s s"%(nomdufichier,len(Signal),Temps[1]-Temps[0],Temps[-1]))
    
    return (Temps, Signal)

"""
Affichage du graphe
"""

import matplotlib.pyplot as plt
from matplotlib import tight_layout
import numpy as np

def tracerECG(Temps, Signal):
    plt.figure("Graphe ECG")
    
    # Ligne horizontale à y = 0
    # plt.axhline(y=0, color='black', linestyle='--')
    # plt.axhline(y=sum(Signal)/len(Signal), color='black', linestyle='--')
    plt.axhline(y=1, color='black', linestyle='--')
    
    # Tracé de l'électrocardiogramme
    plt.title("Electrocardiogramme")
    plt.xlabel("Temps (s)")
    plt.ylabel("Signal (V)")
    plt.grid()
    
    plt.plot(Temps,Signal,'r')
    
    # Centrage autour de y = 0
    maxabs = max(Signal+[-x for x in Signal])
    plt.ylim([-maxabs-1,maxabs+1])
    
    # Affichage du graphe en plein écran directement
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    
    decoupage = decoupageSignal(Temps, Signal)
    for i in decoupage[0]:
        plt.axvline(x=Temps[i], color='black', linestyle='--')
    
    plt.show()

"""
Analyse de fourier
"""

import numpy as np
import scipy.fftpack

def analyseFourier(Temps, Signal):
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

"""
Lissage du signal
"""

from scipy import signal

def lissage(Temps, Signal):
    
    fe = 1/0.003 # Fréquence d'échantillonnage
    f_nyq = fe / 2 # Fréquence de Nyquist
    fc = 30 # Fréquence de coupure
    
    # Filtre de Butterworh
    b, a = signal.butter(4, fc/f_nyq, 'low', analog=False)
    SignalFiltre = signal.filtfilt(b, a, Signal)
    
    # Tracé de l'électrocardiogramme lissé
    plt.figure("ECG lissé")
    
    plt.title("Electrocardiogramme lissé")
    plt.plot(Temps, Signal, color='silver', label="Signal")
    plt.plot(Temps, SignalFiltre, color='red', label='Signal filtré')
    plt.legend()
    plt.xlabel("Temps (s)")
    plt.ylabel("Signal (V)")
    plt.grid()
    plt.show()
    
    return(list(SignalFiltre))

"""
Dérivée du signal
"""

import numpy as np

def derivee(Temps, Signal):
    plt.figure("Dérivée du signal d'ECG")
    
    # Ligne horizontale à y = 0
    plt.axhline(y=0, color='black', linestyle='--')
    
    # Tracé de l'électrocardiogramme
    plt.title("Dérivée du signal de l'électrocardiogramme")
    plt.xlabel("Temps (s)")
    plt.ylabel("Dérivée du signal (V.s-1)")
    plt.grid()
    plt.plot(Temps, Signal, color='silver', label='Signal (V)')
    
    # Centrage autour de y = 0
    maxabs = max(Signal+[-x for x in Signal])
    plt.ylim([-maxabs-1,maxabs+1])
    
    plt.plot(Temps, np.diff([0]+Signal), color='dodgerblue', label='Dérivée du signal')
    plt.show()

# Ancienne méthode
# print(len(list(np.where(np.diff([0]+Signal)==0)[0])))

"""
Découpage du signal
"""

def decoupageSignal(Temps, Signal):
    indicesMax = []
    conditionComplexe = 1.0
    
    # Détection des complexes
    for i in range(1,len(Signal)-2):
        if (Signal[i]>Signal[i-1] and Signal[i]>Signal[i+1]):
            indicesMax.append(i)
    
    indicesComplexes = [i for i in indicesMax if Signal[i]>conditionComplexe]
    
    print("\n")
    for i in range(len(indicesComplexes)):
        print("Complexe n°", i ,"identifié à", Signal[indicesComplexes[i]], "V au temps", Temps[indicesComplexes[i]], "s")
    
    # Découpage des motifs
    Motifs = [] # Format : Temps du 1er complexe, Temps du 2e complexe, Onde T, Onde P
    
    for i in range(len(indicesComplexes)-1): # 4 complexes détectés donnent 3 motifs complets entre ces complexes
        Motifs.append([Temps[indicesComplexes[i]], Temps[indicesComplexes[i+1]]])
    
    # Détections des ondes T et P
    for motif in Motifs:
        # Onde T
        quinzePourcent = motif[0]+0.15*(motif[1]-motif[0])
        trentePourcent = motif[0]+0.30*(motif[1]-motif[0])
        
        ondeT = quinzePourcent
        for i in indicesMax:
            if Temps[i] > quinzePourcent and Temps[i] <= trentePourcent and Temps[i] > ondeT:
                ondeT = Temps[i]
        motif.append(ondeT)
        print("Onde T identifiée au temps",ondeT,"s")
        
        # Onde P
        soixanteDixPourcent = motif[0]+0.70*motif[2]
        quatreVingtCinqPourcent = motif[0]+0.85*motif[2]
        
        ondeP = soixanteDixPourcent
        for i in indicesMax:
            if Temps[i] > soixanteDixPourcent and Temps[i] <= quatreVingtCinqPourcent and Temps[i] > ondeP:
                ondeP = Temps[i]
        motif.append(ondeP)
        print("Onde P identifiée au temps", ondeP, "s")
    
    print(Motifs)
    return(indicesComplexes, Motifs)

"""
Rythme cardiaque
"""
 
def rythmeCardiaque(Temps, Signal, indicesComplexes):
    TempsDesMax = [Temps[i] for i in indicesComplexes]
    Periodes = [TempsDesMax[i+1]-TempsDesMax[i] for i in range(len(TempsDesMax)-1)]
    
    print("\n")
    print(TempsDesMax)
    print((1/(sum(Periodes)/len(Periodes)))*60, "battements par minute")

"""
Raccourci de configuration
"""

a=1
b=1
c=1
d=1
e=1
f=1

Temps, Signal = lireCSV("enregistrements/lycée/EX8.csv")
if a==1: tracerECG(Temps, Signal)
if b==1: analyseFourier(Temps, Signal)
if c==1: SignalFiltre = lissage(Temps, Signal)
if d==1: derivee(Temps, Signal)
if e==1: indicesComplexes = decoupageSignal(Temps, SignalFiltre)[0]
if f==1: rythmeCardiaque(Temps, Signal, indicesComplexes)

"""
Autres
"""
# Affichage d'une ligne verticale ou horizontale :
# plt.axvline(x=1, color='black', linestyle='--')
# plt.axhline(y=1, color='black', linestyle='--')

# Centrage autour de y = 0
# maxabs = max(Signal+[-x for x in Signal])
# plt.ylim([-maxabs-1,maxabs+1])
