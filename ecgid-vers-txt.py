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
    
    PrenomsH = ['Lucas', 'Jules', 'Paul', 'Sacha', 'Tom', 'Martin', 'Enzo', 'Axel', 'Antoine', 'Valentin', 'Samuel', 'Maxence', 'Malo', 'Thomas', 'Oscar', 'William', 'Mike', 'Noah', 'Robert', 'Michel', 'Patrick', 'Alain', 'Didier', 'Hugo', 'Valentin', 'Alexandre', 'Nicolas', 'Benjamin', 'Matthieu', 'Samy', 'Marc', 'Hamza', 'Luc', 'Jacques', 'Philippe', 'Thibaut', 'Steven', 'Albert', 'Harry', 'Juan', 'Serge', 'Benoit', 'Zack', 'Edouard', 'Will', 'Tony', 'Damien', 'Henri', 'Stanislas', 'Yann']
    PrenomsF = ['Emma', 'Camille', 'Sarah', 'Alice', 'Eva', 'Clara', 'Manon', 'Jade', 'Lisa', 'Julie', 'Rose', 'Margot', 'Ambre', 'Claire', 'Lucile', 'Audrey', 'Coline', 'Candice', 'Laura', 'Anne', 'Jacqueline', 'Amy', 'Deborah', 'Enora', 'Flavie', 'Kate', 'Suzie', 'Marie', 'Agathe', 'Suzanne', 'Carla', 'Tess', 'Sophie', 'Fanny', 'Maud', 'Louane', 'Helena', 'Jeanne', 'Charlotte', 'Anais', 'Yasmine', 'Roxane', 'Nina', 'Billie', 'Margaux', 'Albane', 'Romane', 'Brigitte', 'Axelle', 'Ali']

    
    nbH = 0
    nbF = 0
    
    n=1
    for i in range(1,91):
        
        if i<10:
            i = "0"+str(i)
        else:
            i = str(i)
        j = 1
        
        infos = wfdb.rdsamp("ecg-id-database-1.0.0/Person_"+str(i)+"/rec_"+str(j))[1]
        
        
        if infos["comments"][1][5:]=="male":
            sexe = 'male'
            nom = PrenomsH.pop(0)
            nbH +=1
        else:
            sexe = 'female'
            nom = PrenomsF.pop(0)
            nbF +=1
        
        
        
        while path.exists("ecg-id-database-1.0.0/Person_"+str(i)+"/rec_"+str(j)+".dat"):
            signaux, infos = wfdb.rdsamp("ecg-id-database-1.0.0/Person_"+str(i)+"/rec_"+str(j))
            
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
            

            with open('enregistrements/ecg-id-database/'+str(n)+".txt", 'w') as fichier:
                fichier.write("# TIPE ECG 2021 - ECG-ID Database v1.0.0 par Physionet\n")
                fichier.write("# [N° de la personne, N° de l'enregistrement, Période d'échantillonnage, Prénom, Age, Sexe, Date de l'enregistrement]\n")
                fichier.write("["+str(int(i))+", "+str(int(j))+", "+str(TempsEchantillonnage)+", \""+nom+"\", "+age+", \""+sexe+"\", \""+date+"\"]\n")
                for y in Signal:
                    fichier.write(str(y)+"\n")
            with open('enregistrements/ecg-id-database/signaux-deja-filtres-par-physionet/'+str(n)+".txt", 'w') as fichier:
                fichier.write("# TIPE ECG 2021 - ECG-ID Database v1.0.0 par Physionet\n")
                fichier.write("# [N° de la personne, N° de l'enregistrement, Période d'échantillonnage, Prénom, Age, Sexe, Date de l'enregistrement]\n")
                fichier.write("["+str(int(i))+", "+str(int(j))+", "+str(TempsEchantillonnage)+", \""+nom+"\", "+age+", \""+sexe+"\", \""+date+"\"]\n")
                for y in SignalFiltre:
                    fichier.write(str(y)+"\n")
        
            j = j+1
            n = n+1
    print(nbH)
    print(nbF)
        
import prenoms

def selection():
    
    PrenomsH = []
    PrenomsF = []
    
    while len(PrenomsH)<50:
        nom = input()
        PrenomsH.append(str(nom))
        print("Reste",50-len(PrenomsH))
    
    while len(PrenomsF)<50:
        nom = input()
        PrenomsF.append(str(nom))
        print("Reste",50-len(PrenomsF))
    
    print(PrenomsH)
    print(PrenomsF)
