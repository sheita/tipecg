# -*- coding: utf-8 -*-

import wfdb
import matplotlib.pyplot as plt

# record = wfdb.rdrecord('Person_01/rec_1') 
# wfdb.plot_wfdb(record=record)
# print(record.__dict__)
# Temps = [x[0] for x in list(record.__dict__["p_signal"])]
# Signal = [x[1] for x in list(record.__dict__["p_signal"])]
# print(len(Temps), len(Signal))

# plt.plot(Temps, Signal)
# plt.show()
# plt.plot(Temps, SignalFiltre)
# plt.show()

import os.path
from os import path
import json
import names

def dop():
    for i in range(1,91):
        if i<10:
            i = "0"+str(i)
        else:
            i = str(i)
        j = 1
        
        infos = wfdb.rdsamp("Person_"+str(i)+"/rec_"+str(j))[1]
        if infos["comments"][1][5:]=="male":
            sexe = 'male'
        else:
            sexe = 'female'
        nom = names.get_first_name(gender=sexe)
        
        
        
        while path.exists("Person_"+str(i)+"/rec_"+str(j)+".dat"):
            signaux, infos = wfdb.rdsamp("Person_"+str(i)+"/rec_"+str(j))
            
            personID = str(i)
            age = infos["comments"][0][5:]
            if infos["comments"][1][5:]=="male":
                sexe = 'H'
            else:
                sexe = 'F'
            TempsEchantillonnage = 1/infos["fs"]
            date = infos["comments"][2][10:]
            
            Signal = [x[0] for x in signaux]
            SignalFiltre = [x[1] for x in signaux]
            Temps = [i*TempsEchantillonnage for i in range(len(Signal))]
            

            with open('ecg-id-database/'+str(personID)+"-"+nom+"-"+str(j)+".txt", 'w') as fichier:
                fichier.write("# TIPE ECG 2021 - ECG-ID Database v1.0.0 par Physionet\n")
                fichier.write("# [Numéro de la personne, Période d'échantillonnage, Prénom, Age, Sexe, Date de l'enregistrement]\n")
                fichier.write("["+personID+", "+str(TempsEchantillonnage)+", "+nom+", "+age+", "+sexe+", "+date+"]\n")
                for y in Signal:
                    fichier.write(str(y)+"\n")
            with open('ecg-id-database/signaux-deja-filtres-par-physionet/'+str(personID)+"-"+nom+"-"+str(j)+"-F.txt", 'w') as fichier:
                fichier.write("# TIPE ECG 2021 - ECG-ID Database v1.0.0 par Physionet\n")
                fichier.write("# [Numéro de la personne, Période d'échantillonnage, Prénom, Age, Sexe, Date de l'enregistrement]\n")
                fichier.write("["+personID+", "+str(TempsEchantillonnage)+", "+nom+", "+age+", "+sexe+", "+date+"]\n")
                for y in SignalFiltre:
                    fichier.write(str(y)+"\n")
    
            j = j+1
