# -*- coding: utf-8 -*-

"""
Enregistrements

Nom du fichier : Enregistrements[i][0]
Temps : Enregistrements[i][1]
Signal : Enregistrements[i][2]
SignalFiltre : Enregistrements[i][3]
Parametres : Enregistrements[i][4]
Motifs : Enregistrements[i][5]
"""

Enregistrements = []

def data():
    print("La liste Enregistrements contient", len(Enregistrements), "enregistrement(s).")
    print("Le programme travaille actuellement avec l'enregistrement n°"+str(iECG)+".")
    
    for i in range(len(Enregistrements)):
        Temps = Enregistrements[i][1]
        Signal = Enregistrements[i][2]
        print()
        print("• Enregistrement n°"+str(i))
        print("""Fichier lu : "%s" """%(Enregistrements[i][0]))
        print("Nombre de points :", len(Signal))
        print("Période d'échantillonage :", Temps[1]-Temps[0])
        print("Temps d'acquisition total :", Temps[-1], "s")
        print("Rythme cardiaque :", round(rythmeCardiaque(i),2), "bpm")


"""
Lecture de fichier TXT
"""

import json

def lireTXT(cheminVersFichier):
    Temps, Signal = [], []
    
    with open(cheminVersFichier, 'r') as fichier:
        next(fichier)
        next(fichier)
        
        # Solution temporaire en attendant de reconstruire les fichiers avec json.loads()
        Infos = fichier.readline()
        Infos = Infos[1:-1].split(", ")
        Infos[1],Infos[3]=float(Infos[1]),int(Infos[3])
        
        for ligne in fichier:
            Signal.append(float(ligne))
    
    Temps = [Infos[1]*x for x in range(len(Signal))]
    
    print("\nFichier lu : %s \nNombre de points : %s \nPériode d'échantillonage : %s s \nTemps d'acquisition total : %s s"%(cheminVersFichier,len(Signal),Temps[1]-Temps[0],Temps[-1]))
    Enregistrements.append([cheminVersFichier, Temps, Signal, None, None, None])
    
"""
Détection des motifs

Temps du 1er complexe : Motifs[i][0]
Temps du 2e complexe : Motifs[i][1]
Temps de l'onde T : Motifs[i][2]
Temps de l'onde P : Motifs[i][3]
"""

def detectionMotifs(iECG):
    Temps = Enregistrements[iECG][1]
    Signal = Enregistrements[iECG][2]
    SignalFiltre = Enregistrements[iECG][3]
    Signal = list(SignalFiltre)
    
    print()
    
    # Détermination de la condition complexe
    SignalMax = max(Signal)
    conditionComplexe = 0.7*SignalMax
    print("La condition complexe choisie est à", conditionComplexe, "V.")

    # Détection de l'ensemble des maximums locaux
    indicesMax = []
    
    for i in range(1,len(Signal)-2):
        if (Signal[i]>Signal[i-1] and Signal[i]>Signal[i+1]):
            indicesMax.append(i)
    
    # Détection des complexes parmi les maximums locaux
    indicesComplexes = [i for i in indicesMax if Signal[i]>conditionComplexe]
    
    for i in range(len(indicesComplexes)):
        print("Complexe n°", i ,"identifié à", Signal[indicesComplexes[i]], "V au temps", Temps[indicesComplexes[i]], "s")
        
        
    # Découpage des motifs
    Motifs = [] # Format : 
    
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
    
    Enregistrements[iECG][4] = [conditionComplexe, indicesMax, indicesComplexes]
    Enregistrements[iECG][5] = Motifs


"""
Lissage du signal
"""

from scipy import signal

def lissage(iECG):
    Temps = Enregistrements[iECG][1]
    Signal = Enregistrements[iECG][2]
    PeriodeEchantillonnage = 0.003
    
    fe = 1/PeriodeEchantillonnage # Fréquence d'échantillonnage
    f_nyq = fe / 2 # Fréquence de Nyquist
    fc = 30 # Fréquence de coupure
    
    # Filtre de Butterworh
    b, a = signal.butter(4, fc/f_nyq, 'low', analog=False)
    SignalFiltre = signal.filtfilt(b, a, Signal)
    
    Enregistrements[iECG][3] = list(SignalFiltre)    


"""
Affichage du graphe
"""

import matplotlib.pyplot as plt
from matplotlib import tight_layout
import numpy as np

def tracerECG(iECG):
    Temps = Enregistrements[iECG][1]
    Signal = Enregistrements[iECG][2]
    SignalFiltre = Enregistrements[iECG][3]
    indicesComplexes = Enregistrements[iECG][4][2]
    conditionComplexe = Enregistrements[iECG][4][0]
    
    plt.figure("Graphe de l'ECG n°"+str(iECG))
    
    # Couleur de fond de la fenêtre d'affichage du graphe en blanc
    plt.rcParams["figure.facecolor"] = 'w'
    
    # Ligne horizontale à y = 0
    # plt.axhline(y=0, color='black', linestyle='--')
    # plt.axhline(y=sum(Signal)/len(Signal), color='black', linestyle='--')
    plt.axhline(y=conditionComplexe, color='black', linestyle='--')
    
    # Tracé de l'électrocardiogramme
    plt.title("Electrocardiogramme")
    plt.plot(Temps, Signal, color='silver', label="Signal (V)")
    plt.plot(Temps, SignalFiltre, color='red', label='Signal filtré (V)')
    plt.grid()
    plt.legend()
    plt.show()
    
    # Centrage autour de y = 0
    maxabs = max(Signal+[-x for x in Signal])
    plt.ylim([-maxabs-1,maxabs+1])
    
    # Affichage du graphe en plein écran directement
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    
    # Affichage des lignes horizontales aux abscisses des complexes
    for i in indicesComplexes:
        plt.axvline(x=Temps[i], color='black', linestyle='--')
    
    plt.show()


"""
Analyse de fourier
"""

import numpy as np
import scipy.fftpack

def analyseFourier(i):
    Temps = Enregistrements[iECG][1]
    Signal = Enregistrements[iECG][2]
    
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
Dérivée du signal
"""

import numpy as np

def derivee(iECG):
    Temps = Enregistrements[iECG][1]
    Signal = Enregistrements[iECG][2]
    
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


"""
Rythme cardiaque
"""
 
def rythmeCardiaque(iECG):
    Temps = Enregistrements[iECG][1]
    Signal = Enregistrements[iECG][2]
    indicesComplexes = Enregistrements[iECG][4][2]
    
    if indicesComplexes==None:
        return(0)
    
    TempsDesMax = [Temps[i] for i in indicesComplexes]
    Periodes = [TempsDesMax[i+1]-TempsDesMax[i] for i in range(len(TempsDesMax)-1)]
    
    return((1/(sum(Periodes)/len(Periodes)))*60)
        
"""
Configuration du programme
"""
lireTXT("enregistrements/lycée/EX8.txt")
lireTXT("enregistrements/ecg-id-database/01-Paul-1.txt")

Analyser = [0,1]
Afficher = [0,1]

for iECG in Analyser:
    lissage(iECG)
    detectionMotifs(iECG)

for iECG in Afficher:
    tracerECG(0)

"""
Autres
"""
# Affichage d'une ligne verticale ou horizontale :
# plt.axvline(x=1, color='black', linestyle='--')
# plt.axhline(y=1, color='black', linestyle='--')

# Centrage autour de y = 0
# maxabs = max(Signal+[-x for x in Signal])
# plt.ylim([-maxabs-1,maxabs+1])
