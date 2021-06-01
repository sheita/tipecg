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
plt.rcParams["font.size"]=16
plt.rcParams["xtick.labelsize"]=12
plt.rcParams["ytick.labelsize"]=12
plt.rcParams["xtick.major.pad"]=8
plt.rcParams["ytick.major.pad"]=8
plt.rcParams["axes.labelpad"]=10
plt.rcParams["savefig.dpi"]=300
import warnings
warnings.simplefilter("ignore", category=RuntimeWarning)

"""
[Enregistrements]

Les données concernant chaque électrocardiogramme sont accessibles par Enregistrements[i]

Meta : Enregistrements[i][0]
(nom du fichier, numéro de la personne, numéro de l'enregistrement, 
période d'échantillonnage, prénom, age, sexe, date de l'enregistrement)

Signaux : Enregistrements[i][1]
(signal, signal filtré)

ParametresAnalyse : Enregistrements[i][2]
(condition de complexe, indices max, indices des complexes)

Motifs : Enregistrements[i][3]
(pour chaque motif : temps du 1er complexe, temps du 2e complexe, onde T, onde P)

ParametresAuthentification : Enregistrements[i][4]
(segment TP, segment TQRS, segment QRSP)
"""

Enregistrements = []

def data():
    print("La liste Enregistrements contient", len(Enregistrements), "enregistrements.")
    print()
    
    # Sur chaque ligne : affichage du numéro de l'enregistrement, du prénom, de l'âge, du sexe, de la durée de l'acquisition et de la période d'échantillonage
    for i in range(len(Enregistrements)):
        Meta = Enregistrements[i][0]
        PeriodeEchantillonnage = Meta[3]
        Signaux = Enregistrements[i][1]
        Signal = Signaux[0]
        Temps = [n*PeriodeEchantillonnage for n in range(len(Signal))]
        print("[%s] %s (%s%s)     durée %ss     période %ss"%(str(i), Meta[4], Meta[5], Meta[6],Temps[-1],PeriodeEchantillonnage))

"""
[Lecture de fichier TXT]
"""

def lireTXT(cheminVersFichier):
    Signal = []
    
    # Lecture du fichier texte : chaque fichier texte contient les données d'un unique enregistrement ainsi que toutes les caractéristiques utiles à son étude.
    with open(cheminVersFichier, 'r') as fichier:
        # On passe les deux premières lignes qui expliquent ce que contient le fichier
        next(fichier)
        next(fichier)
        # La liste Meta sur la 3eme ligne est récupérée comme telle par le programme
        Meta = json.loads(fichier.readline())
        
        # On récupère les valeurs du potentiel une à une
        for ligne in fichier:
            Signal.append(float(ligne))
    
    # On ajoute au début de la liste Meta l'emplacement du fichier
    Meta.insert(0, cheminVersFichier)

    # Les ECG sont ajoutés à la liste Enregistrements dans l'ordre où les fichiers ont été appelés par cette fonction
    Enregistrements.append([Meta, [Signal, None], None, None])
    
"""
[Détection des motifs]
"""

def detectionMotifs(iECG):
    Meta = Enregistrements[iECG][0]
    Signaux = Enregistrements[iECG][1]
    
    # La détection des motifs est réalisée sur le signal filtré, ou le signal non filtré le cas échéant
    Signal = Signaux[1]
    if Signal == None : Signal = Signaux[0]
    
    PeriodeEchantillonnage = Meta[3]
    Temps = [n*PeriodeEchantillonnage for n in range(len(Signal))]
    
    # Méthode de la condition complexe
    
    # Détection de l'ensemble des maximums locaux
    indicesMax = []
    
    for i in range(1,len(Signal)-2):
        if (Signal[i]>Signal[i-1] and Signal[i]>Signal[i+1]):
            indicesMax.append(i)
   
    # Choix de la valeur de la condition complexe à 70% de la hauteur de la valeur de potentiel maximale de l'ECG
    SignalMax = max(Signal)
    conditionComplexe = 0.7*SignalMax
    
    indicesComplexes = [i for i in indicesMax if Signal[i]>conditionComplexe]

    # Méthode de la dérivée seconde
    SignalDerive2nde = derivee(iECG, n=2, tracerDerivee=False)
    
    # Détection de l'ensemble des minimums locaux sur le graphe de la dérivée 2nde du signal
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
        # On recherche l'onde T entre 10% et 40% de chaque motif de l'ECG (entre chaque complexe)
        dixPourcent = motif[0]+0.10*(motif[1]-motif[0])
        quarantePourcent = motif[0]+0.40*(motif[1]-motif[0])
        
        # Initialisation de l'indice
        indiceOndeT = indicesMax.index(min(indicesMax))
        
        # Si un des indices max est plus haut que le précédent, il est remplacé
        for i in indicesMax:
            if Temps[i] > dixPourcent and Temps[i] <= quarantePourcent and Temps[i] > Temps[indiceOndeT]:
                indiceOndeT = i
        
        # On garde en mémoire l'indice de l'onde et la zone dans laquelle elle a été recherchée pour l'afficher sur le graphe
        motif.append([indiceOndeT, dixPourcent, quarantePourcent])
        
        # Même principe pour l'onde P entre 60% et 90%
        soixantePourcent = motif[0]+0.60*(motif[1]-motif[0])
        quatreVingtDixPourcent = motif[0]+0.90*(motif[1]-motif[0])
        indiceOndeP = indicesMax.index(min(indicesMax))
        for i in indicesMax:
            if Temps[i] > soixantePourcent and Temps[i] <= quatreVingtDixPourcent and Temps[i] > Temps[indiceOndeP]:
                indiceOndeP = i
        motif.append([indiceOndeP, soixantePourcent, quatreVingtDixPourcent])
    
    # On sauvegarde les paramètres qui ont servis à la détection ainsi que les motifs détectés dans la liste principale
    Enregistrements[iECG][2] = [conditionComplexe, indicesMax, indicesComplexes, indicesMin]
    Enregistrements[iECG][3] = Motifs
    
"""
[Lissage du signal]
"""

def lissage(iECG):
    # Application d'un filtre passe-bas de fréquence de coupure 20Hz au signal
    Meta = Enregistrements[iECG][0]
    Signaux = Enregistrements[iECG][1]
    Signal = Signaux[0]
    PeriodeEchantillonnage = Meta[3]
    
    fe = 1/PeriodeEchantillonnage # Fréquence d'échantillonnage
    f_nyq = fe / 2 # Fréquence de Nyquist
    fc = 20 # Fréquence de coupure
    
    # Filtre de Butterworth
    b, a = signal.butter(4, fc/f_nyq, 'low', analog=False)
    SignalFiltre = signal.filtfilt(b, a, Signal)
    
    # Le signal filtré est sauvegardé dans la liste principale, en plus du signal initial
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
    
    PeriodeEchantillonnage = Meta[3]
    Temps = [n*PeriodeEchantillonnage for n in range(len(Signal))]
    
    # Le prénom de la personne dont on trace l'enregistrement apparait dans le titre de la fenêtre
    plt.figure("Electrocardiogramme de "+Meta[4])
    
    # Affichage des axes en x et y
    plt.axhline(y=0, color='black', linestyle='-', alpha=0.5)
    plt.axvline(x=0, color='black', linestyle='-', alpha=0.5)
    
    # L'affichage du graphe est différent si le signal est filtré, filtré et analysé ou ni l'un ni l'autre
    estFiltré = SignalFiltre != None
    estAnalysé = Motifs != None
    
    # Tracé de l'électrocardiogramme
    if estAnalysé:
        # Cas où le signal a été filtré par la fonction lissage() et les motifs ont été détectés par detectionMotifs()
        for i in range(len(Motifs)):
            motif = Motifs[i]
            
            # Si la longueur d'un motif est aberrante, le fond du motif s'affiche en rouge afin d'identifier plus facilement les erreurs de détection
            if motif[1]-motif[0]>=0.5 and motif[1]-motif[0]<=1:
                plt.axvspan(motif[0], motif[1], facecolor='cornflowerblue', alpha=0.04)
            else:
                plt.axvspan(motif[0], motif[1], facecolor='red', alpha=0.12)

            # Affichage des complexes à chaque extremité du motif
            plt.axvline(x=motif[0], color='hotpink', linestyle='-', alpha=1)
            plt.axvline(x=motif[1], color='hotpink', linestyle='-', alpha=1)

            # Affichage de l'onde T dans la zone où il a été recherché
            plt.axvspan(motif[2][1], motif[2][2], facecolor='royalblue', alpha=0.1)
            plt.axvline(x=Temps[motif[2][0]], color='cornflowerblue', linestyle='-')
            
            # Affichage de l'onde P dans la zone où il a été recherché
            plt.axvspan(motif[3][1], motif[3][2], facecolor='darkorange', alpha=0.1)
            plt.axvline(x=Temps[motif[3][0]], color='coral', linestyle='-')
        
        # Tracé du signal (en gris) et du signal filtré (en rouge)
        plt.plot(Temps, Signal, color='silver', label="Signal", zorder=0.5)
        plt.plot(Temps, SignalFiltre, color='red', label="Signal filtré", zorder=10)
        
        # Affichage de la condition de complexe (pour illustrer le diapo)
        conditionComplexe = ParametresAnalyse[0]
        #plt.axhline(y=conditionComplexe, color='cornflowerblue', linestyle='-')
    
    elif estFiltré:
        # Cas où le signal a été filtré par la fonction lissage()
        plt.plot(Temps, Signal, color='silver', label="Signal")
        plt.plot(Temps, SignalFiltre, color='red', label="Signal filtré")
    
    else:
        # Cas où le signal est brut
        plt.plot(Temps, Signal, color='red', label="Signal")
    
    # Affichage des légendes du graphe
    plt.ylabel("Potentiel (V)")
    plt.xlabel("Temps (s)")
    plt.legend()
    
    # Centrage du graphe autour de y = 0 pour une meilleure visibilité
    maxabs = max(Signal+[-x for x in Signal])
    plt.ylim([-maxabs-1,maxabs+1])
    
    # Affichage du graphe dans une fenêtre
    plt.show()

"""
[Graphes de comparaison des paramètres]
"""
def comparaisonParametres(iECGs, GrapheAnonyme=0, ParametresAnonymes=[], AfficherGraphes=True):
    # L'ECG anonyme ne sera pas utilisé pour construire les graphes de comparaison des paramètres 
    
    # Chaque liste contiendra, pour chaque individu de la base de donnée, l'ensemble des mesures des segments pour chaque motif détectés
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
        
        # Pour tous les autres graphes que le graphe anonyme, on mesure les segments TP, T-QRS et QRS-P
        if not GrapheAnonyme==iECG:
            NumeroPersonne = Meta[1]
            
            for i in range(len(Motifs)):
                segmentTP = Temps[Motifs[i][3][0]]-Temps[Motifs[i][2][0]]
                SegmentsTP[NumeroPersonne].append(segmentTP)
                
                segmentTQRS = Motifs[i][1]-Temps[Motifs[i][2][0]]
                SegmentsTQRS[NumeroPersonne].append(segmentTQRS)
                
                segmentQRSP = Temps[Motifs[i][3][0]]-Motifs[i][0]
                SegmentsQRSP[NumeroPersonne].append(segmentQRSP)
    
    # On liste les moyennes et les écarts-types de la longueur des segments pour l'ensemble des individus
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
        
    # Si l'affichage des graphes est demandé à l'appel de la fonction, on trace les graphes de comparaison (oui par défaut)
    if AfficherGraphes:
        plt.figure("Comparaison des longueurs TP")
        
        # Affichage des axes en x et en y
        plt.axhline(y=0, color='black', linestyle='-', alpha=0.5)
        plt.axvline(x=0, color='black', linestyle='-', alpha=0.5)
        
        # Chaque individu obtient une couleur aléatoire sur le graphe
        couleurs = [(random.random(), random.random(), random.random()) for i in range(311)]
        
        # Pour chaque individu
        for i in range(len(MoyennesTP)):
            if EcartsTypesTP[i]<0.3:
                # Affichage de points correspondant à la moyenne de la longueur des segments TP pour chaque individu
                plt.scatter(i, MoyennesTP[i], color=couleurs[i])
                # Affichage d'une barre d'erreur qui représente l'écart-type correspondant à cette moyenne
                plt.errorbar(i, MoyennesTP[i], yerr = EcartsTypesTP[i], fmt = 'none', capsize = 5, color = couleurs[i], zorder = 3)
        
        # On affiche une barre rouge horizontale pour l'authentification lorsque les paramètres du graphe anonyme sont donnés en appel de la fonction
        if ParametresAnonymes:
            plt.axhline(y=ParametresAnonymes[0], color='red', zorder=0.5)
        
        # Affichage des légendes du graphe
        plt.ylabel("Longueur du segment TP")
        plt.xlabel("N° de la personne")
        plt.xticks([10*x for x in range(10)])
        
        # Affichage du graphe
        plt.show()
    
        
        # Même principe pour le segment T-QRS
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
                
        
        # Même principe pour le segment QRS-P
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
    
    # Si l'affichage des graphes n'est pas demandé, on retourne les données calculées
    else :
        return [[MoyennesTP, MoyennesTQRS, MoyennesQRSP],[EcartsTypesTP, EcartsTypesTQRS, EcartsTypesQRSP]]

"""
[Authentification]
"""

def authentification(iECGs, iGrapheAnonyme=0):
    Meta = Enregistrements[iGrapheAnonyme][0]
    Signaux = Enregistrements[iGrapheAnonyme][1]
    Signal = Signaux[0]
    Motifs = Enregistrements[iGrapheAnonyme][3]
    PeriodeEchantillonnage = Meta[3]
    Temps = [n*PeriodeEchantillonnage for n in range(len(Signal))]
    
    # Même principe que la fonction comparaisonParametres() mais pour un unique ECG
    SegmentsTP = []
    SegmentsTQRS = []
    SegmentsQRSP = []
    
    for i in range(len(Motifs)):
        segmentTP = SegmentsTP.append(Temps[Motifs[i][3][0]]-Temps[Motifs[i][2][0]])
        segmentTQRS = SegmentsTQRS.append(Motifs[i][1]-Temps[Motifs[i][2][0]])
        segmentQRSP = SegmentsQRSP.append(Temps[Motifs[i][3][0]]-Motifs[i][0])
    
    # La moyenne des segments TP, T-QRS et QRS-P est déterminée pour ce graphe uniquement
    moyenneTP = np.mean(SegmentsTP)
    moyenneTQRS = np.mean(SegmentsTQRS)
    moyenneQRSP = np.mean(SegmentsQRSP)
    
    ParametresAnonymes = [moyenneTP, moyenneTQRS, moyenneQRSP]
    
    # La fonction comparaisonParametres() est appelée sur l'ensemble des ECG à l'exception de l'ECG anonyme
    ParametresConnus = comparaisonParametres(iECGs, iGrapheAnonyme, ParametresAnonymes, AfficherGraphes=False)
    
    SegmentsConnus = ParametresConnus[0]
    EcartsTypesConnus = ParametresConnus[1]
    
    print("L'authentification est réalisée sur le graphe anonyme n°%s de %s (id %s) \
sur la base de l'ensemble des électrocardiogrammes différents de celui-ci"%(iECG, Meta[4], Meta[1]))

    # L'ensemble des individus de la base de données sont initialement des candidats potentiels à l'authentification du graphe anonyme
    PersonnesPossibles = [True for i in range(91)]
    
    # On élimine toutes les personnes dont les longueurs des segments ne sont pas compatibles avec celles du graphe anonyme
    for i in range(len(PersonnesPossibles)):
        if not (moyenneTP<SegmentsConnus[0][i]+EcartsTypesConnus[0][i] and moyenneTP>SegmentsConnus[0][i]-EcartsTypesConnus[0][i]\
        and moyenneTQRS<SegmentsConnus[1][i]+EcartsTypesConnus[1][i] and moyenneTQRS>SegmentsConnus[1][i]-EcartsTypesConnus[1][i]\
        and moyenneQRSP<SegmentsConnus[2][i]+EcartsTypesConnus[2][i] and moyenneQRSP>SegmentsConnus[2][i]-EcartsTypesConnus[2][i]):
            PersonnesPossibles[i]=False
    
    # On obtient une liste de personnes qui sont possiblement celles ayant enregistré l'ECG anonyme
    indicesPersonnesPossibles = [i for i in range(len(PersonnesPossibles)) if not PersonnesPossibles[i]]
    
    return (Meta[1] in indicesPersonnesPossibles, len(indicesPersonnesPossibles))

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
    
    T = PeriodeEchantillonnage
    yf = scipy.fftpack.fft(SignalNp)
    xf = np.linspace(0, 1.0/(2.0*T), 250)
    
    # Affichage de l'analyse de Fourier sous forme d'histogramme en amplitude
    plt.stem(xf, 2.0/N * np.abs(yf[:N//2]), markerfmt=' ')
    plt.xlabel("Fréquence (Hz)")
    plt.ylabel("Amplitude")
    plt.show()

def derivee(iECG, n=1, tracerDerivee=True):
    Signaux = Enregistrements[iECG][1]
    
    # La détection des motifs est réalisée sur le signal filtré, ou le signal non filtré le cas échéant
    Signal = Signaux[1]
    if Signal == None : Signal = Signaux[0]
    
    PeriodeEchantillonnage = Enregistrements[iECG][0][3]
    Temps = [n*PeriodeEchantillonnage for n in range(len(Signal))]
    
    # La dérivée n-ième d'un signal utilise des tableaux numpy
    SignalDerive = list(np.diff([0]+Signal, n))
    TempsDerive = [n*PeriodeEchantillonnage for n in range(len(SignalDerive))]
    
    # Si l'affichage du graphe est demandé à l'appel de la fonction, on trace la dérivée
    if tracerDerivee:
        plt.figure("Dérivée du signal de l'enregistrement n°"+str(iECG))
        
        # Affichage des axes en x et y
        plt.axhline(y=0, color='black', linestyle='-', alpha=0.5)
        plt.axvline(x=0, color='black', linestyle='-', alpha=0.5)
        
        # Tracé de l'électrocardiogramme
        plt.xlabel("Temps (s)")
        plt.plot(Temps, Signal, color='silver', label='Potentiel (V)')
        
        # Affichage du seuil à partir duquel la dérivée seconde correspond à un complexe
        plt.axhline(y=min(SignalDerive)*0.65, color='green', linestyle='-')
        
        # Centrage autour de y = 0
        maxabs = max(Signal+[-x for x in Signal])
        plt.ylim([-maxabs-1,maxabs+1])
        
        # Affichage de la dérinée n-ième du signal
        plt.plot(TempsDerive, SignalDerive, color='dodgerblue', label='Dérivée '+str(n)+['ère','nde','ème'][min([n-1,2])]+' du signal')
        plt.legend()
        plt.show()
        
    # Si l'affichage du graphe est demandé à l'appel de la fonction, on retourne les données calculées
    else:
        return(SignalDerive)

"""
[Rythme cardiaque]
"""
 
def rythmeCardiaque(iECG):
    # La fonction s'exécute seulement lorsque des motifs ont été détectés par la fonction detectionMotifs()
    Motifs = Enregistrements[iECG][3]
    if Motifs:
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
Afficher = [0]

for iECG in Lisser: lissage(iECG)
for iECG in Analyser: detectionMotifs(iECG)
for iECG in Afficher: tracerECG(iECG)

print("Terminé. \n\nUtilisez la fonction data() pour lister les enregistrements qui ont été chargés")
