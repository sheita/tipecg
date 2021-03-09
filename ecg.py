# -*- coding: utf-8 -*-
# TIPE ECG 2021
# /!\ Voir la configuration du fichier à la fin du programme

"""
[Enregistrements]

Le données concernant chaque enregistrement sont accessibles par Enregistrements[i]

Meta : Enregistrements[i][0]
(nom du fichier, numéro de la personne, numéro de l'enregistrement, 
période d'échantillonnage, prénom, age, sexe, date de l'enregistrement)

Signaux : Enregistrements[i][1]
(signal, signal filtré)

ParametresAnalyse : Enregistrements[i][2]
(condition de complexe, indices max, indices complexes)

Motifs : Enregistrements[i][3]
(pour chaque motif : temps du 1er complexe, temps du 2e complexe, onde T, onde P)
"""

Enregistrements = []

def data():
    print("La liste Enregistrements contient", len(Enregistrements), "enregistrements.")
    
    for i in range(len(Enregistrements)):
        Meta = Enregistrements[i][0]
        PeriodeEchantillonnage = Meta[3]
        Signaux = Enregistrements[i][1]
        Signal = Signaux[0]
        Temps = [n*PeriodeEchantillonnage for n in range(len(Signal))]
    
        print()
        print("[Enregistrement n°"+str(i)+"]")
        print("""Fichier lu : "%s" """%(Meta[0]))
        print("Nombre de points :", len(Signal))
        print("Période d'échantillonnage :", PeriodeEchantillonnage)
        print("Temps d'acquisition total :", Temps[-1], "s")
        print("Rythme cardiaque :", round(rythmeCardiaque(i),2), "bpm")

"""
[Lecture de fichier TXT]
"""

import json

def lireTXT(cheminVersFichier):
    Signal = []
    
    with open(cheminVersFichier, 'r') as fichier:
        next(fichier)
        next(fichier)
        Meta = json.loads(fichier.readline())
        
        for ligne in fichier:
            Signal.append(float(ligne))
        
    Meta.insert(0, cheminVersFichier)
    Temps = [Meta[3]*x for x in range(len(Signal))]
    
    print("\nFichier lu : %s \nNombre de points : %s \nPériode d'échantillonage : %s s \nTemps d'acquisition total : %s s"%(Meta[0],len(Signal),Temps[1]-Temps[0],Temps[-1]))
    Enregistrements.append([Meta, [Signal, None], None, None])
    
"""
[Détection des motifs]
"""

def detectionMotifs(iECG):
    Meta = Enregistrements[iECG][0]
    Signaux = Enregistrements[iECG][1]
    
    Signal = Signaux[1]
    if Signal == None : Signal = Signaux[0]
    
    PeriodeEchantillonnage = Meta[3]
    Temps = [n*PeriodeEchantillonnage for n in range(len(Signal))]
    
    print()
    # Détection de l'ensemble des maximums locaux
    indicesMax = []
    
    for i in range(1,len(Signal)-2):
        if (Signal[i]>Signal[i-1] and Signal[i]>Signal[i+1]):
            indicesMax.append(i)
   
    # Méthode de la condition complexe
    SignalMax = max(Signal)
    conditionComplexe = round(0.7*SignalMax,2)
    indicesComplexes = [i for i in indicesMax if Signal[i]>conditionComplexe]
    print("[ECG n°"+str(iECG)+"] La condition complexe choisie est", conditionComplexe, "V")

    # Méthode de la dérivée seconde
    SignalDerive2nde = derivee(iECG, n=2, tracerDerivee=False)
    SignalDerive2ndeMin = min(SignalDerive2nde)
    indicesMin = []
    
    for i in range(1,len(SignalDerive2nde)-2):
        if (SignalDerive2nde[i]<SignalDerive2nde[i-1] and SignalDerive2nde[i]<SignalDerive2nde[i+1]):
            indicesMin.append(i)
    
    indicesComplexes = [i for i in indicesMin if SignalDerive2nde[i]<SignalDerive2ndeMin*0.65]
        
    # Découpage des motifs
    Motifs = []
    
    for i in range(len(indicesComplexes)-1): # 4 complexes détectés donnent 3 motifs complets entre ces complexes
        Motifs.append([Temps[indicesComplexes[i]], Temps[indicesComplexes[i+1]]])
    
    # Détections des ondes T et P
    for motif in Motifs:
        # Onde T
        quinzePourcent = motif[0]+0.15*(motif[1]-motif[0])
        trentePourcent = motif[0]+0.30*(motif[1]-motif[0])
        
        indiceOndeT = indicesMax[0]
        for i in indicesMax:
            if Temps[i] > quinzePourcent and Temps[i] <= trentePourcent and Temps[i] > Temps[indiceOndeT]:
                indiceOndeT = i
        motif.append([indiceOndeT, quinzePourcent, trentePourcent])
        
        # Onde P
        soixanteDixPourcent = motif[0]+0.70*(motif[1]-motif[0])
        quatreVingtCinqPourcent = motif[0]+0.85*(motif[1]-motif[0])
        
        indiceOndeP = indicesMax[0]
        for i in indicesMax:
            if Temps[i] > soixanteDixPourcent and Temps[i] <= quatreVingtCinqPourcent and Temps[i] > Temps[indiceOndeP]:
                indiceOndeP = i
        motif.append([indiceOndeP, soixanteDixPourcent, quatreVingtCinqPourcent])
    
    print("[ECG n°"+str(iECG)+"]", len(indicesComplexes), "complexes identifiés, soit", len(Motifs), "motifs")
    
    Enregistrements[iECG][2] = [conditionComplexe, indicesMax, indicesComplexes, indicesMin]
    Enregistrements[iECG][3] = Motifs


"""
[Lissage du signal]
"""

from scipy import signal

def lissage(iECG):
    # Fonction qui permet de filtrer l'enregistrement choisi
    Signaux = Enregistrements[iECG][1]
    Signal = Signaux[0]
    PeriodeEchantillonnage = Enregistrements[iECG][0][3]
    
    fe = 1/PeriodeEchantillonnage # Fréquence d'"chantillonnage
    f_nyq = fe / 2 # Fréquence de Nyquist
    fc = 30 # Fréquence de coupure
    
    # Filtre de Butterworth
    b, a = signal.butter(4, fc/f_nyq, 'low', analog=False)
    SignalFiltre = signal.filtfilt(b, a, Signal)
    
    Enregistrements[iECG][1][1] = list(SignalFiltre)    


"""
[Affichage du graphe]
"""

import numpy as np
import matplotlib.pyplot as plt
#from matplotlib import tight_layout

def tracerECG(iECG):
    Signaux = Enregistrements[iECG][1]
    ParametresAnalyse = Enregistrements[iECG][2]
    Motifs = Enregistrements[iECG][3]
    Signal = Signaux[0]
    SignalFiltre = Signaux[1]
    indicesMax, indicesMin = Enregistrements[iECG][2][1], Enregistrements[iECG][2][3]
    
    PeriodeEchantillonnage = Enregistrements[iECG][0][3]
    Temps = [n*PeriodeEchantillonnage for n in range(len(Signal))]
    
    plt.figure("Graphe de l'ECG n°"+str(iECG))
    plt.title("Electrocardiogramme de l'enregistrement n°"+str(iECG))
    
    # Couleur de fond de la fenêtre d'affichage du graphe en blanc
    plt.rcParams["figure.facecolor"] = 'w'
    
    # Axes en x et y
    plt.axhline(y=0, color='black', linestyle='-')
    plt.axvline(x=0, color='black', linestyle='-')
    
    estFiltré = SignalFiltre != None
    estAnalysé = Motifs != None
    
    # Tracé de l'électrocardiogramme
    if estAnalysé:
        conditionComplexe = ParametresAnalyse[0]
        plt.axhline(y=conditionComplexe, color='cornflowerblue', linestyle='-')
        plt.plot(Temps, Signal, color='silver', label="Signal (V)")
        plt.plot(Temps, SignalFiltre, color='red', label='Signal filtré (V)')
        for i in range(len(Motifs)):
            motif = Motifs[i]
            if i%2==0:
                plt.axvspan(motif[0], motif[1], facecolor='darkgrey', alpha=0.12)
            else:
                plt.axvspan(motif[0], motif[1], facecolor='darkgrey', alpha=0.01)
            
            plt.axvspan(motif[2][1], motif[2][2], facecolor='lightblue', alpha=0.2)
            plt.axvspan(motif[3][1], motif[3][2], facecolor='lightblue', alpha=0.2)
            plt.axvline(x=Temps[motif[2][0]], color='orange', linestyle='-')
            plt.axvline(x=Temps[motif[3][0]], color='orange', linestyle='-')
        
#        TempsMax = [Temps[i] for i in indicesMax]
#        TempsMin = [Temps[i] for i in indicesMin]
#        for i in TempsMax:
#            plt.axvline(x=i, color='red', linestyle='-')
#        for i in TempsMin:
#            plt.axvline(x=i, color='green', linestyle='-')
    
    elif estFiltré:
        plt.plot(Temps, Signal, color='silver', label="Signal (V)")
        plt.plot(Temps, SignalFiltre, color='red', label='Signal filtré (V)')
            
    else:
        Signal = Signaux[0]
        plt.plot(Temps, Signal, color='red', label="Signal (V)")
    
    # plt.grid()
    plt.ylabel("Signaux")
    plt.xlabel("Temps (s)")
    plt.legend()
    
    # Centrage autour de y = 0
    maxabs = max(Signal+[-x for x in Signal])
    plt.ylim([-maxabs-1,maxabs+1])
    
    # Affichage du graphe en plein écran directement
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    
    plt.show()


"""
[Pertinence des paramètres]
"""
def pertinenceParametres(iECG):
    Meta = Enregistrements[iECG][0]
    Signaux = Enregistrements[iECG][1]
    Signal = Signaux[0]
    Motifs = Enregistrements[iECG][3]
    PeriodeEchantillonnage = Meta[3]
    Temps = [n*PeriodeEchantillonnage for n in range(len(Signal))]
    
    SegmentsTP = []
    
    for motif in Motifs:
        SegmentsTP.append(Temps[motif[3][0]] - Temps[motif[2][0]])
    
    plt.ylabel("Longueur du segment TP")
    plt.xlabel("Motifs")
    plt.plot(range(len(SegmentsTP)),SegmentsTP)
    plt.ylim([0,max(SegmentsTP)*1.1])
    plt.show()


"""
[Analyse de Fourier]
"""

import numpy as np
import scipy.fftpack

def analyseFourier(iECG):
    Signaux = Enregistrements[iECG][1]
    Signal = Signaux[1]
    Meta = Enregistrements[iECG][0]
    PeriodeEchantillonnage = Meta[3]

    plt.figure("Analyse de Fourier")

    # L'analyse de Fourier utilise des tableaux numpy
    SignalNp = np.array(Signal)
    
    N = 500
    T = PeriodeEchantillonnage
    yf = scipy.fftpack.fft(SignalNp)
    xf = np.linspace(0, 1.0/(2.0*T), N/2)
    
    plt.stem(xf, 2.0/N * np.abs(yf[:N//2]), markerfmt=' ')
    plt.title("Analyse de Fourier du signal ECG")
    plt.xlabel("Fréquence (Hz)")
    plt.ylabel("Amplitude")
    plt.show()


"""
[Dérivée du signal]
"""

import numpy as np

def derivee(iECG, n=1, tracerDerivee=True):
    Signaux = Enregistrements[iECG][1]
    Signal = Signaux[1]
    if Signal == None : Signal = Signaux[0]
    PeriodeEchantillonnage = Enregistrements[iECG][0][3]
    Temps = [n*PeriodeEchantillonnage for n in range(len(Signal))]
    
    SignalDerive = list(np.diff([0]+Signal, n))
    TempsDerive = [n*PeriodeEchantillonnage for n in range(len(SignalDerive))]
    
    if tracerDerivee:
        plt.figure("Dérivée du signal d'ECG")
        
        # Ligne horizontale à  y = 0
        plt.axhline(y=0, color='black', linestyle='--')
        
        # Tracé de l'électrocardiogramme
        plt.title("Dérivée du signal de l'électrocardiogramme")
        plt.xlabel("Temps (s)")
        #plt.grid()
        plt.plot(Temps, Signal, color='silver', label='Signal (V)')
        
        # Centrage autour de y = 0
        maxabs = max(Signal+[-x for x in Signal])
        plt.ylim([-maxabs-1,maxabs+1])
        
        plt.plot(TempsDerive, SignalDerive, color='dodgerblue', label='Dérivée '+str(n)+['ère','nde','ème'][min([n-1,2])]+' du signal')
        plt.legend()
        
        plt.axhline(y=min(SignalDerive)*0.65, color='green', linestyle='-')
        
        plt.show()
        
    else:
        return(SignalDerive)

"""
[Rythme cardiaque]
"""
 
def rythmeCardiaque(iECG):
    Motifs = Enregistrements[iECG][3]
    
    tempsPremierComplexe = Motifs[0][0]
    tempsDernierComplexe = Motifs[-1][1]
    battementsParSeconde = (tempsDernierComplexe-tempsPremierComplexe)/len(Motifs)
    battementsParMinute = battementsParSeconde*60
    
    return(battementsParMinute)
        
"""
[Configuration du programme]
"""
lireTXT("enregistrements/lycée/EX8.txt")
for i in [1]:
    lireTXT("enregistrements/ecg-id-database/"+str(i)+".txt")

Lisser = [0]
Analyser = [0]
Afficher = [0]

for iECG in Lisser: lissage(iECG)
for iECG in Analyser: detectionMotifs(iECG)
for iECG in Afficher: tracerECG(iECG)
    
