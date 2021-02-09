# -*- coding: utf-8 -*-

Signal = []
Temps = []

with open("enregistrements/lycée/EX8.csv", 'r') as fichier:
    next(fichier)
    
    for x in fichier:
        Signal.append(float(x.strip().replace(",",".").split(";")[1]))
        Temps.append(float(x.strip().replace(",",".").split(";")[0]))
    print(Signal)

with open("enregistrements/lycée/EX8.txt", 'w') as fichier:
    fichier.write("# TIPE ECG 2021 - Enregistrements réalisés au lycée\n")
    fichier.write("# [Numéro de la personne, Période d'échantillonnage, Prénom, Age, Sexe, Date de l'enregistrement]\n")
    fichier.write("""["00", 0.003, "Clément", 20, "H", "13.10.2020"]\n""")
    for y in Signal:
        fichier.write(str(y)+"\n")
