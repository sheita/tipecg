# -*- coding: utf-8 -*-
# TIPE ECG 2021 - Candidat 190
# /!\ Voir la configuration du fichier à la fin du programme

import numpy as np
import matplotlib.pyplot as plt
import json # Module pour la lecture de fichiers
import random

# Modules scientifiques (dérivée, analyse de Fourier, filtrage)
from scipy import signal 
import scipy.fftpack 

# Amélioration de la lisibilité des graphes (grande police d'écriture, couleur de fond en blanc)
plt.rcParams["figure.facecolor"] = 'w'
plt.rcParams["font.size"]= 16
plt.rcParams["xtick.labelsize"]= 12
plt.rcParams["ytick.labelsize"]= 12
plt.rcParams["xtick.major.pad"]= 8
plt.rcParams["ytick.major.pad"]= 8
plt.rcParams["axes.labelpad"]= 10
plt.rcParams["savefig.dpi"]= 300
import warnings
warnings.simplefilter("ignore", category=RuntimeWarning)


"""
[Enregistrements]

Le données concernant chaque électrocardiogramme sont accessibles par Enregistrements[i]

Meta : Enregistrements[i][0]
(nom du fichier, numéro de la personne, numéro de l'enregistrement, 
période d'échantillonnage, prénom, age, sexe, date de l'enregistrement)

Signaux : Enregistrements[i][1]
(signal, signal filtré)

ParametresAnalyse : Enregistrements[i][2]
(condition de complexe, indices max, indices des complexes)

Motifs : Enregistrements[i][3]
(pour chaque motif : temps du 1er complexe, temps du 2e complexe, onde T, onde P)
"""

Enregistrements = []

def data():
    print("La liste Enregistrements contient", len(Enregistrements), "enregistrements.")
    print()
    
    for i in range(len(Enregistrements)):
        Meta = Enregistrements[i][0]
        PeriodeEchantillonnage = Meta[3]
        Signaux = Enregistrements[i][1]
        Signal = Signaux[0]
        Temps = [n*PeriodeEchantillonnage for n in range(len(Signal))]
        print("[%s] %s (%s%s)     durée %ss     période %ss"%(str(i), Meta[4], Meta[5], Meta[6],Temps[-1],PeriodeEchantillonnage))

"""
[Lecture de fichier TXT]

Chaque fichier texte contient les données d'un unique enregistrement ainsi que toutes les caractéristiques utiles à leur étude.

"""

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
    
    #print("\nFichier lu : %s \nNombre de points : %s \nPériode d'échantillonage : %s s \nTemps d'acquisition total : %s s"%(Meta[0],len(Signal),Temps[1]-Temps[0],Temps[-1]))
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
    
    # Détection de l'ensemble des maximums locaux
    indicesMax = []
    
    for i in range(1,len(Signal)-2):
        if (Signal[i]>Signal[i-1] and Signal[i]>Signal[i+1]):
            indicesMax.append(i)
   
    # Méthode de la condition complexe
    SignalMax = max(Signal)
    conditionComplexe = round(0.7*SignalMax,2)
    indicesComplexes = [i for i in indicesMax if Signal[i]>conditionComplexe]
    #print("[ECG n°"+str(iECG)+"] La condition complexe choisie est", conditionComplexe, "V")

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
        quinzePourcent = motif[0]+0.10*(motif[1]-motif[0])
        trentePourcent = motif[0]+0.40*(motif[1]-motif[0])
        
        indiceOndeT = indicesMax.index(min(indicesMax))
        for i in indicesMax:
            if Temps[i] > quinzePourcent and Temps[i] <= trentePourcent and Temps[i] > Temps[indiceOndeT]:
                indiceOndeT = i
        motif.append([indiceOndeT, quinzePourcent, trentePourcent])
        
        # Onde P
        soixanteDixPourcent = motif[0]+0.60*(motif[1]-motif[0])
        quatreVingtCinqPourcent = motif[0]+0.90*(motif[1]-motif[0])
        
        indiceOndeP = indicesMax.index(min(indicesMax))
        for i in indicesMax:
            if Temps[i] > soixanteDixPourcent and Temps[i] <= quatreVingtCinqPourcent and Temps[i] > Temps[indiceOndeP]:
                indiceOndeP = i
        motif.append([indiceOndeP, soixanteDixPourcent, quatreVingtCinqPourcent])
    
    #print("[ECG n°"+str(iECG)+"]", len(indicesComplexes), "complexes identifiés, soit", len(Motifs), "motifs")
    
    Enregistrements[iECG][2] = [conditionComplexe, indicesMax, indicesComplexes, indicesMin]
    Enregistrements[iECG][3] = Motifs

"""
[Lissage du signal]
"""

def lissage(iECG):
    # Fonction qui permet de filtrer l'enregistrement choisi
    Signaux = Enregistrements[iECG][1]
    Signal = Signaux[0]
    PeriodeEchantillonnage = Enregistrements[iECG][0][3]
    
    fe = 1/PeriodeEchantillonnage # Fréquence d'échantillonnage
    f_nyq = fe / 2 # Fréquence de Nyquist
    fc = 20 # Fréquence de coupure
    
    # Filtre de Butterworth
    b, a = signal.butter(4, fc/f_nyq, 'low', analog=False)
    SignalFiltre = signal.filtfilt(b, a, Signal)
    
    Enregistrements[iECG][1][1] = list(SignalFiltre)    


"""
[Affichage du graphe]
"""

def tracerECG(iECG):
    Meta = Enregistrements[iECG][0]
    Signaux = Enregistrements[iECG][1]
    ParametresAnalyse = Enregistrements[iECG][2]
    Motifs = Enregistrements[iECG][3]
    Signal = Signaux[0]
    SignalFiltre = Signaux[1]
    
    PeriodeEchantillonnage = Enregistrements[iECG][0][3]
    Temps = [n*PeriodeEchantillonnage for n in range(len(Signal))]
    
    plt.figure("Electrocardiogramme de "+Meta[4])
    
    # Axes en x et y
    plt.axhline(y=0, color='black', linestyle='-', alpha=0.5)
    plt.axvline(x=0, color='black', linestyle='-', alpha=0.5)
    
    estFiltré = SignalFiltre != None
    estAnalysé = Motifs != None
    
    # Tracé de l'électrocardiogramme
    if estAnalysé:
        indicesMax, indicesMin = ParametresAnalyse[1], ParametresAnalyse[3]
        conditionComplexe = ParametresAnalyse[0]
        plt.axhline(y=conditionComplexe, color='cornflowerblue', linestyle='-')
        plt.plot(Temps, Signal, color='silver', label="Signal")
        plt.plot(Temps, SignalFiltre, color='red', label="Signal filtré")
        for i in range(len(Motifs)):
            motif = Motifs[i]
            if motif[1]-motif[0]>=0.5 and motif[1]-motif[0]<=1: # La longueur d'un motif est d'environ 0.75
                plt.axvspan(motif[0], motif[1], facecolor='green', alpha=0.05)
            else:
                plt.axvspan(motif[0], motif[1], facecolor='red', alpha=0.12)
            

            # Complexes seuls
            plt.axvline(x=motif[0], color='darkgreen', linestyle='-', alpha=0.1)
            plt.axvline(x=motif[1], color='darkgreen', linestyle='-', alpha=0.1)

            plt.axvspan(motif[2][1], motif[2][2], facecolor='orange', alpha=0.05)
            plt.axvspan(motif[3][1], motif[3][2], facecolor='orange', alpha=0.05)
            plt.axvline(x=Temps[motif[2][0]], color='orange', linestyle='-')
            plt.axvline(x=Temps[motif[3][0]], color='orange', linestyle='-')
            
    
    elif estFiltré:
        plt.plot(Temps, Signal, color='silver', label="Signal")
        plt.plot(Temps, SignalFiltre, color='red', label="Signal filtré")
    
    else:
        Signal = Signaux[0]
        plt.plot(Temps, Signal, color='red', label="Signal")
    
    plt.ylabel("Potentiel (V)")
    plt.xlabel("Temps (s)")
    plt.legend()
    
    # Centrage autour de y = 0
    maxabs = max(Signal+[-x for x in Signal])
    plt.ylim([-maxabs-1,maxabs+1])
    
    plt.show()


"""
[Graphe de comparaison des paramètres]

Par exemple comparaisonParametres(range(0,311), saufLeGrapheNumero=12)

Le graphe qui n'est pas utilisé pour établir ce graphe sera utilisé comme
électrocardiogramme "anonyme" pour lequel on cherchera à qui il correspond
"""
def comparaisonParametres(iECGs, GrapheAnonyme=0, ParametresAnonymes=[]):
    SegmentsTP = [[] for i in range(91)]
    SegmentsTQRS = [[] for i in range(91)]
    SegmentsQRSP = [[] for i in range(91)]
    
    for iECG in iECGs:
        Meta = Enregistrements[iECG][0]
        Signaux = Enregistrements[iECG][1]
        Signal = Signaux[0]
        Motifs = Enregistrements[iECG][3]
        PeriodeEchantillonnage = Meta[3]
        Temps = [n*PeriodeEchantillonnage for n in range(len(Signal))]
        
        if not GrapheAnonyme==iECG:
            NumeroPersonne = Meta[1]
            
            for i in range(len(Motifs)):
                segmentTP = Temps[Motifs[i][3][0]]-Temps[Motifs[i][2][0]]
                SegmentsTP[NumeroPersonne].append(segmentTP)
                
                segmentTQRS = Motifs[i][1]-Temps[Motifs[i][2][0]]
                SegmentsTQRS[NumeroPersonne].append(segmentTQRS)
                
                segmentQRSP = Temps[Motifs[i][3][0]]-Motifs[i][0]
                SegmentsQRSP[NumeroPersonne].append(segmentQRSP)
    
    MoyennesTP = []
    MoyennesTQRS = []
    MoyennesQRSP = []
    
    EcartsTypesTP = []
    EcartsTypesTQRS = []
    EcartsTypesQRSP = []
    
    for i in range(len(SegmentsTP)):
        MoyennesTP.append(np.mean(SegmentsTP[i]))
        EcartsTypesTP.append(np.std(SegmentsTP[i]))
        
    for i in range(len(SegmentsTQRS)):
        MoyennesTQRS.append(np.mean(SegmentsTQRS[i]))
        EcartsTypesTQRS.append(np.std(SegmentsTQRS[i]))
    
    for i in range(len(SegmentsQRSP)):
        MoyennesQRSP.append(np.mean(SegmentsQRSP[i]))
        EcartsTypesQRSP.append(np.std(SegmentsQRSP[i]))
        
        
    plt.figure("Comparaison des longueurs TP")
        
    plt.axhline(y=0, color='black', linestyle='-', alpha=0.5)
    plt.axvline(x=0, color='black', linestyle='-', alpha=0.5)
    
    couleurs = [(random.random(), random.random(), random.random()) for i in range(311)]
    
    for i in range(len(MoyennesTP)):
        if EcartsTypesTP[i]<0.3:
            plt.scatter(i, MoyennesTP[i], color=couleurs[i])
            plt.errorbar(i, MoyennesTP[i], yerr = EcartsTypesTP[i], fmt = 'none', capsize = 5, color = couleurs[i], zorder = 3)
            
    if ParametresAnonymes:
        plt.axhline(y=ParametresAnonymes[0], color='red', zorder=0.5)
        
    plt.ylabel("Longueur du segment TP")
    plt.xlabel("N° de la personne")
    plt.xticks([10*x for x in range(10)])
    plt.show()
    
    
    plt.figure("Comparaison des longueurs TQRS")
        
    plt.axhline(y=0, color='black', linestyle='-', alpha=0.5)
    plt.axvline(x=0, color='black', linestyle='-', alpha=0.5)
    
    for i in range(len(MoyennesTQRS)):
        if EcartsTypesTQRS[i]<0.3:
            plt.scatter(i, MoyennesTQRS[i], color=couleurs[i])
            plt.errorbar(i, MoyennesTQRS[i], yerr = EcartsTypesTQRS[i], fmt = 'none', capsize = 5, color = couleurs[i], zorder = 3)
    
    if ParametresAnonymes:
        plt.axhline(y=ParametresAnonymes[1], color='red', zorder=0.5)
            
    plt.ylabel("Longueur du segment TQRS")
    plt.xlabel("N° de la personne")
    plt.xticks([10*x for x in range(10)])
    plt.show()
            
    
    plt.figure("Comparaison des longueurs QRSP")
        
    plt.axhline(y=0, color='black', linestyle='-', alpha=0.5)
    plt.axvline(x=0, color='black', linestyle='-', alpha=0.5)
    
    for i in range(len(MoyennesQRSP)):
        if EcartsTypesQRSP[i]<0.3:
            plt.scatter(i, MoyennesQRSP[i], color=couleurs[i])
            plt.errorbar(i, MoyennesQRSP[i], yerr = EcartsTypesQRSP[i], fmt = 'none', capsize = 5, color = couleurs[i], zorder = 3)
            
    if ParametresAnonymes:
        plt.axhline(y=ParametresAnonymes[2], color='red', zorder=0.5)

    plt.ylabel("Longueur du segment QRSP")
    plt.xlabel("N° de la personne")
    plt.xticks([10*x for x in range(10)])
    plt.show()

"""
[Authentification]
"""
def authentification(iECG):
    Meta = Enregistrements[iECG][0]
    Signaux = Enregistrements[iECG][1]
    Signal = Signaux[0]
    Motifs = Enregistrements[iECG][3]
    PeriodeEchantillonnage = Meta[3]
    Temps = [n*PeriodeEchantillonnage for n in range(len(Signal))]
    
    for i in range(len(Motifs)):
        segmentTP = Temps[Motifs[i][3][0]]-Temps[Motifs[i][2][0]]
        segmentTQRS = Motifs[i][1]-Temps[Motifs[i][2][0]]
        segmentQRSP = Temps[Motifs[i][3][0]]-Motifs[i][0]
    
    ParametresAnonymes = [segmentTP, segmentTQRS, segmentQRSP]
    
    comparaisonParametres(range(0,311), iECG, ParametresAnonymes)
    
    print("L'authentification est réalisée sur le graphe anonyme n°%s de %s (id %s) \
sur la base de l'ensemble des électrocardiogrammes différents de celui-ci"%(iECG, Meta[4], Meta[1]))


"""
[Analyse de Fourier et Dérivée du signal]
"""

def analyseFourier(iECG):
    Signaux = Enregistrements[iECG][1]
    Signal = Signaux[1]
    Meta = Enregistrements[iECG][0]
    PeriodeEchantillonnage = Meta[3]

    plt.figure("Analyse de Fourier de l'enregistrement n°"+str(iECG))

    # L'analyse de Fourier utilise des tableaux numpy
    SignalNp = np.array(Signal)
    
    N = 500
    T = PeriodeEchantillonnage
    yf = scipy.fftpack.fft(SignalNp)
    xf = np.linspace(0, 1.0/(2.0*T), N//2)
    
    plt.stem(xf, 2.0/N * np.abs(yf[:N//2]), markerfmt=' ')
    plt.xlabel("Fréquence (Hz)")
    plt.ylabel("Amplitude")
    plt.show()

def derivee(iECG, n=1, tracerDerivee=True):
    Signaux = Enregistrements[iECG][1]
    Signal = Signaux[1]
    if Signal == None : Signal = Signaux[0]
    PeriodeEchantillonnage = Enregistrements[iECG][0][3]
    Temps = [n*PeriodeEchantillonnage for n in range(len(Signal))]
    
    SignalDerive = list(np.diff([0]+Signal, n))
    TempsDerive = [n*PeriodeEchantillonnage for n in range(len(SignalDerive))]
    
    if tracerDerivee:
        plt.figure("Dérivée du signal de l'enregistrement n°"+str(iECG))
        
        # Ligne horizontale à  y = 0
        plt.axhline(y=0, color='black', linestyle='--')
        
        # Tracé de l'électrocardiogramme
        plt.xlabel("Temps (s)")
        plt.plot(Temps, Signal, color='silver', label='Potentiel (V)')
        
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

Utilisez la fonction lireTXT(emplacementFichier) pour ajouter des enregistrements au programme.
Les enregistrements prennent un indice iECG dans leur ordre d'arrivée.
"""

print("Chargement des fichiers...", end=' ')


# Enregistrement personnel (vu dans le diaporama)
lireTXT("enregistrements/lycée/EX8.txt")

# Enregistrements de la base de données Physionet
for i in range(1,311):
    lireTXT("enregistrements/ecg-id-database/"+str(i)+".txt")

Lisser =  range(0,311)
Analyser = range(0,311)
Afficher = []

for iECG in Lisser: lissage(iECG)
for iECG in Analyser: detectionMotifs(iECG)
for iECG in Afficher: tracerECG(iECG)


print("Terminé. \n\nUtilisez la fonction data() pour lister les enregistrements qui ont été chargés")
