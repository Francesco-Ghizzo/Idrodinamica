# -*- coding: utf-8 -*-
#!/usr/bin/env python
# =================================================================
#
# ------------------
# Metodo di Engelund
# ------------------
#
# Schema per il calcolo della scala di deflusso di sezioni complesse
#
# Ia esercitazione di Idrodinamica
# Anno accademico 2019/2020
# Luca Adami
#
# ==================================================================


# =========================
# Import Pacchetti e Moduli
# =========================
from __future__ import division # Evita la divisione tra interi
import numpy as np # Libreria numerica
from matplotlib import pyplot as plt # Libreria grafica
import os

# ================
# Dati dell'utente
# ================
file_input = ''
file_output = ''
vpoints = 100 # N punti discretizzazione verticale
g = 9.81 # Accelerazione di gravita'
iF = 0.003 # Pendenza
NG = 2 # Numero di punti di Gauss
tol = 1e-03 # Tolleranza nel calcolo della profondita' critica
MAXITER = 100 # Numero massimo di iterazioni nel calcolo della profondità critica
fs = 16 # Font Size x grafici

# ========
# Funzioni
# ========

def MotoUniforme( iF, y, z, ks, Y, NG=2 ):
    '''
    Calcola i parametri di moto uniforme per assegnato tirante

    Argomenti
    ---------

    iF: float
       pendenza del canale
    y: numpy.ndarray
      coordinate trasversali dei punti della sezione
    z: numpy.ndarray
      coordinate verticali dei punti della sezione
    ks: numpy.ndarray
      coefficienti di scabrezza dei punti della sezione
    Y: float
      profondità alla quale calcolare i parametri di moto uniforme
    NG: int [default=2]
      numero di punti di Gauss

    Output
    ------
    Q: float
      portata alla quale si realizza la profondità Y di moto uniforme
    Omega: float
      area sezione bagnata alla profondita' Y
    b: float
      larghezza superficie libera alla profondita' Y
    alpha: float
      coefficiente di ragguaglio dell'energia alla profondita' Y
    beta: float
      coefficiente di ragguaglio della qdm alla profondita' Y
    '''

    # Inizializzo
    Omega = 0 # Area bagnata
    b = 0 # Larghezza superficie libera
    num_alpha = 0 # Numeratore di alpha
    num_beta = 0 # Numeratore di beta
    den = 0 # Base del denominatore di alpha e beta

    Yi = Y - (z-z.min())  # Distribuzione trasversale della profondita'
    N = Yi.size # Numero di punti sulla trasversale

    # N punti trasversali -> N-1 intervalli (trapezi)
    for i in range( N-1 ): # Per ogni trapezio

        if not (Yi[i] <= 0 and Yi[i+1] <= 0): # Non considero i tratti con tirante idraulico negativo (al di sopra della superficie dell'acqua)
        
        # calcolo i valori di y[i], y[i+1], Yi[i], Yi[i+1], dy e dz per i triangoli laterali

            if Yi[i] < 0:
                dy = y[i+1]-y[i]
                dy = dy*Yi[i+1]/(z[i]-z[i+1])
                Yi[i] = 0
                dz = Yi[i+1]
            elif Yi[i+1] < 0:
                dy = y[i+1]-y[i]
                dy = dy*Yi[i]/((z[i+1]-z[i]))
                Yi[i+1] = 0
                dz = Yi[i]
            else:

        # ... calcolare i valori per il singolo trapezio
        
                dy = y[i+1]-y[i]
                dz = z[i+1]-z[i]

            cos_phi = dy/np.sqrt(dy**2 + dz**2)
    
    # ... calcolare gli integrali
    
            b = b + dy
            Omega = Omega + (Yi[i]+Yi[i+1])*0.5*dy
            num_alpha = num_alpha + (cos_phi**2)*((ks[i]**3)*(Yi[i]**3) + ((ks[i+1]**3)*(Yi[i+1]**3)))*0.5*dy
            num_beta = num_beta + (cos_phi**(4/3))*((ks[i]**2)*(Yi[i]**(7/3)) + (ks[i+1]**2)*(Yi[i+1]**(7/3)))*0.5*dy
            den = den + (cos_phi**(2/3))*(ks[i]*Yi[i]**(5/3)+ks[i+1]*Yi[i+1]**(5/3))*0.5*dy

    Q = den*np.sqrt(iF)

    if den != 0: # Evito la divisione per 0  
        alpha = (Omega**2)*num_alpha/(den**3)
        beta = Omega*num_beta/(den**2)
    else:
        alpha = 0
        beta = 0

    return Q, Omega, b, alpha, beta



# =================
# Codice principale
# =================

# Carica File di Input
# --------------------
cs = np.loadtxt( file_input ) # Carica Sezione
y, z, ks = cs[:, [0]], cs[:, [1]], cs[:, [2]] #Carico i valori delle colonne del file di input .dat come vettori e assegno i valori a y, z e ks

# Calcolo della scala di deflusso
# -------------------------------
h = np.linspace( z.min(), z.max(), vpoints+1 )[1:] # Array dei livelli della superficie libera
Y = h - z.min() # Array dei tiranti

# Inizializzo gli array di output
# -------------------------------
Q = np.zeros( vpoints ) # Portata
Omega = np.zeros( vpoints ) # Area
b = np.zeros( vpoints ) # Larghezza superficie libera
alpha = np.zeros( vpoints ) # Coefficiente di ragguaglio dell'energia
beta = np.zeros( vpoints ) # Coefficiente di ragguaglio della qta' di moto
Yc = np.zeros( vpoints ) # Tirante critico


# Ciclo sui tiranti
# -----------------
# Per ciascun livello della superficie libera calcolo la portata defluente in condizioni
# di moto uniforme, i coefficienti di ragguaglio e la geometria della sezione bagnata 
for n in range( vpoints ):
    
    # Calcola i parametri di moto uniforme assegnato il tirante
    Q[n], Omega[n], b[n], alpha[n], beta[n] = MotoUniforme( iF, y, z, ks, Y[n], NG=NG )
    

# Salva File di Output
# --------------------
out_table = np.array([ Y, Q, Yc, b, Omega, alpha, beta ]).T
np.savetxt( file_output, out_table )


# Crea i grafici
# --------------

# Grafico 1: Profilo Sezione
plt.scatter(y, z)
plt.title("Profilo dell'alveo")
plt.xlabel("y[m]")
plt.ylabel("z[m]")
plt.plot(y, z)
plt.show()

# Grafico 2: Scala di Deflusso
plt.title("Scala di deflusso")
plt.xlabel("Q[m3/s]")
plt.ylabel("Y[m]")
plt.plot(Q, Y)
plt.show()

# Grafico 3: Coefficiente di Ragguaglio Alpha
plt.title("Coefficiente di ragguaglio α")
plt.xlabel("α[-]")
plt.ylabel("Y[m]")
plt.plot(alpha, Y)
plt.show()

# Grafico 4: Coefficiente di Ragguaglio Beta
plt.title("Coefficiente di ragguaglio ß")
plt.xlabel("ß[-]")
plt.ylabel("Y[m]")
plt.plot(beta, Y)
plt.show()
