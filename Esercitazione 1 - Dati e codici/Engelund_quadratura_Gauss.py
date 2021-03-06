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
def GaussPoints(NG):
    '''
    Funzione per il calcolo dei punti e dei pesi di Gauss
    
    Argomenti
    ---------
    NG: int
       numero di punti di Gauss

    Output
    ------
    p: numpy.ndarray
      array dei punti di Gauss
    w: numpy.ndarray
      array dei pesi
    '''
    p, w = None, None
    if NG==2:
        p = np.array([ -1/np.sqrt(3),
                       +1/np.sqrt(3) ])
        w = np.array([ 1,
                       1 ])
    elif NG==3:
        p = np.array([ -np.sqrt(3)/np.sqrt(5),
                       0,
                       +np.sqrt(3)/np.sqrt(5)
        ])
        w = np.array([ 5/9,
                       8/9,
                       5/9
        ])
    elif NG==4:
        p = np.array([ -np.sqrt(525+70*np.sqrt(30))/35,
                       -np.sqrt(525-70*np.sqrt(30))/35,
                       +np.sqrt(525-70*np.sqrt(30))/35,
                       +np.sqrt(525+70*np.sqrt(30))/35
        ])
        w = np.array([ (18-np.sqrt(30))/36,
                       (18+np.sqrt(30))/36,
                       (18+np.sqrt(30))/36,
                       (18-np.sqrt(30))/36
        ])

    return p, w


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
    # Punti e pesi di Gauss
    xj, wj = GaussPoints( NG ) # Calcola i punti e i pesi di Gauss

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
        
                dy = y[i+1]-y[i]
                dz = z[i+1]-z[i]

            cos_phi = dy/np.sqrt(dy**2 + dz**2)
    
    # ... calcolare gli integrali
    
            b = b + dy
            Omega = Omega + (Yi[i]+Yi[i+1])*0.5*dy
            ksj = (ks[i+1]-ks[i])*0.5*np.array(xj) + (ks[i]+ks[i+1])*0.5
            Yj = (Yi[i+1]-Yi[i])*0.5*np.array(xj) + (Yi[i]+Yi[i+1])*0.5
            num_alpha = num_alpha + dy*0.5*(cos_phi**2)*np.sum(wj*(ksj**3)*(Yj**3))
            num_beta = num_beta + dy*0.5*(cos_phi**(4/3))*np.sum(wj*(ksj**2)*(Yj**(7/3)))
            den = den + dy*0.5*(cos_phi**(2/3))*np.sum(wj*ksj*Yj**(5/3))


    Q = den*np.sqrt(iF)

    if den != 0: # Evito la divisione per 0  
        alpha = (Omega**2)*num_alpha/(den**3)
        beta = Omega*num_beta/(den**2)
    else:
        alpha = 0
        beta = 0

    return Q, Omega, b, alpha, beta


def Critica( iF, y, z, ks, Q, MAXITER=100, tol=1e-03, NG=2 ):
    '''
    Calcolo della profondita' critica ad assegnata portata

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
    Q: float
      portata
    MAXITER: int [default=100]
      numero massimo di iterazioni per l'algoritmo dicotomico
    tol: float
      tolleranza sull'eccesso di energia per l'algoritmo dicotomico
    NG: int [default=2]
      numero di punti di Gauss nel calcolo dei parametri

    Output
    ------
    Ym: float
       profondita' critica calcolata con il metodo dicotomico
    '''
    Ya = 1e-06
    Yb = z.max() - z.min()

    # Calcolo della profondita' critica con il metodo dicotomico
    # La funzione da annullare e' quella per l'eccesso di carico specifico sulla sezione

    for i in range(MAXITER):
        Ym = (Ya+Yb)/2
        fYm = Energia( iF, y, z, ks, Ym, Q, NG)
        fYa = Energia( iF, y, z, ks, Ya, Q, NG)
        if abs(fYm) < tol:
            break
        else:
            if fYm*fYa < 0:
                Yb = Ym
            else:
                Ya = Ym    
        
    return Ym


def Energia( iF, y, z, ks, Y, Q, NG=2 ):
    '''
    Eccesso di energia rispetto alla critica
    Funzione da annullare per trovare le condizioni critiche
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
      profondita'
    Q: float
      portata
    NG: int [default=2]
      numero di punti di Gauss nel calcolo dei parametri
    '''
    
    _, Omega, b, alpha, _ = MotoUniforme( iF, y, z, ks, Y, NG=NG )    
    Fr2 = (Q**2)*b/(9.81*Omega**3)
    return alpha*Fr2 - 1


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
    
    # Calcola il livello critico associato alla portata corrispondente
    # al livello di moto uniform corrente
    Yc[n] = Critica( iF, y, z, ks, Q[n], NG=NG, MAXITER=MAXITER, tol=tol ) 


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
plt.plot(Q, Yc)
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
