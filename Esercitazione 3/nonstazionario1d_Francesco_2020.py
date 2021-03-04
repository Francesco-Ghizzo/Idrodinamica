# -*- coding: utf-8 -*_
#!/usr/bin/env python

# ============================================================================== #
#                                                                                #
#                             Corso di Idrodinamica                              #
#                                                                                #
#                           Anno Accademico 2019/2020                            #
#                        III Esercitazione - Moto Vario                          #
#                         Autore: Luca Adami                                     #
#                                                                                #
#          Il codice risolve il sistema di equazioni di De Saint Venant          #
#                                                                                #
# Descrizione: Lo script risolve le equazioni unidimensionali del moto nel caso  #
#              semplificato di alveo rettangolare largo a pendenza costante.     #
#              I tre scenari possibili (diga, paratoia, piena) individuano tre   #
#              doverse combinazioni di condizioni iniziali e al contorno.        #
#              Alcuni parametri sono condivisi dai tre casi di studio, mentre    #
#              altri sono specifici per il singolo caso di studio.               #
#              Tutti i parametri sono esplicitati nel primo blocco di codice,    #
#              mentre il resto del codice è autonomo                             #
#                                                                                #
# ============================================================================== #


# DA FARE PER INIZIARE:
# 	- sostituire FORCE con LF nel principale

# ======================
# Importazione Pacchetti
# ======================
import numpy as np
import matplotlib.pyplot as plt
import os

# =====================
# Parametri dell'utente
# =====================
# Parametri Generali
# ------------------
problema = 'diga' # Opzioni: 'diga', 'paratoia', 'piena'
g = 9.81 # Accelerazione di gravita' [m/s^2]
iF = 0.001 # Pendenza costante del canale [-]
L = 100e3 # Lunghezza del canale [m] (Utilizzare una lunghezza di 300m per la diga, 100e3 per la piena)
IMAX = 100e3 # Numero di celle del dominio [-] (Utilizzare celle di 1m -> IMAX = L)
NMAX = 100 # Numero massimo di iterazioni temporali [-]
ks = 30 # Coefficiente di scabrezza del fondo [m^(1/3)/s]
TIMEOUT = ... # Tempo finale della simulazione [s]
CFL = ... # Coefficiente di Courant-Friedrichs-Levy # [-]
# Parametri Diga
# --------------
xdiga = ... # Posizione dello sbarramento [m]
YL = ... # Livello a sinistra dello sbarramento [m]
YR = ... # Livello a destra dello sbarramento [m]
# Parametri Paratoia
# ------------------
q0 = 1 # Portata unitaria di base [m^2/s]
t_chiusura = 5 # Tempo di chiusura della paratoia [s]
# Parametri Piena
# ---------------
file_idrogramma = 'idrogramma.txt' # File contenente due colonne x tempo e portata dell'idrogramma


# ========
# Funzioni
# ========

def CC( t, dt, U, problema ):
    '''Applicazione delle condizioni al contorno'''
    Y, q = U
    
    if problema == 'diga':
        Y[0], q[0] = Y[1], q[1] # Condizioni di monte
        Y[-1], q[-1] = Y[-2], q[-2] # Condizioni di valle
        
elif problema == 'paratoia':
        Y[0], q[0] = Y[1], q[1] # Condizioni di monte
        if t >= t_chiusura:
            Y[-1], q[-1] = Y[-2], q[-2] # Condizioni di valle
        else:
            Y[-1], q[-1] = Y[-2], -q[-2] # Condizioni di valle

    elif problema == 'piena':
        if t < t_hydro[0]: # Se la piena deve ancora arrivare...
            qt = ...
            Yt = UniFlow( qt, ks, iF )
        elif t > t_hydro[-1]: # Se la piena è già terminata...
            qt = ...
            Yt = ...
        else: # Sta transitando l'onda di piena
            qt = np.interp( t, t_hydro, q_hydro ) # vedere significato di np.interp!
            Fr = ... # Numero di Froude
            if Fr>=1: # Il controllo è a monte
                Yt = ...
            else: # Il controllo è a valle
                # Cerco il piede della caratteristica partendo dalla cella 0
                # ed andando a valle finchè non lo trovo.
                # L'iterazione non dovrebbe superare le prime 2-3 celle
                for i in range( x.size ):
                    ...
                    xR = ...
                    if x[i]-dx/2 <= xR <= x[i]+dx/2: break # Ho trovato xR!
                # Uso dei valori interpolati di Y e q
                YR = np.interp( xR, x, Y )
                qR = np.interp( xR, x, q )
                ...
                Yt = ...
        Y[0], q[0] = Yt, qt
        Y[-1], q[-1] = ...
        
    U[0,:] = Y
    U[1,:] = q
    
    return U

def UniFlow( q, ks, iF ):
    '''Calcola la profondita' di moto uniforme'''
    Y = (q/(ks*np.sqrt(iF)))**(3/5) #per alveo rettangolare largo
    return Y

def PhysFlux( U ):
    '''Calcolo dei flussi fisici'''
    Y, q = U
    Flux = np.array([q, q**2/Y + g*Y**2/2])
    return Flux

def LaxFriedrichs( U, dt, dx ):
    '''Flusso Numerico di Lax-Friedrichs'''
    NumFlux = 0.5*(PhysFlux()+PhysFlux())
    return NumFlux

def LaxWendroff( U, dt, dx ):
    '''Flusso Numerico di Lax-Wendroff'''
    ...
    return NumFlux

def FORCE( U, dt, dx ):
    '''Flusso Numerico FORCE di Toro'''
    ...
    return NumFlux

def Source( U ):
    '''Termine Sorgente'''
    ...
    return S

def RK2( S, U, dt ):
    '''Schema di Runge-Kutta del secondo ordine'''
    ...
    return ...


# =============================================
# Creazione della griglia e condizioni iniziali
# =============================================
# Mesh e quote del fondo
# ----------------------
dx = L / IMAX
x = np.linspace( -dx, L+dx, IMAX+2 ) # aggiungere due celle per condizioni al contorno
b = -iF*x # Quota del fondo (posta a 0 per x=0)

# Condizioni iniziali
# -------------------
if problema == 'diga':
    Y = np.where( ... ) # Profondità
    q = np.zeros( ... ) # Portata
elif problema == 'paratoia':
    Y = np.ones( IMAX+2 ) * UniFlow( q0, ks, iF )
    q = np.ones( IMAX+2 ) * q0
elif problema == 'piena':
    t_hydro, q_hydro = np.loadtxt( file_idrogramma ).T
    q0 = q_hydro[0]
    Y = ...
    q = ...
else:
    # Esci dal codice con un errore
    raise ValueError, '<problema> deve essere "diga", "paratoia" o "piena", non "%s"!' % problema

U = np.array([Y,q]) # Variabili conservate

time = 0 # Tempo corrente
times = [] # Lista dei tempi da salvare alla fine


# ==================================
# Predisposizione cartella di output
# ==================================
cartella_output = 'output_%s' % problema # Nome cartella di output
timesfile = '%s_times.txt' % problema # Nome file dove salvare i tempi
if not os.path.isdir( cartella_output ): os.mkdir( cartella_output ) # Se la cartella di output non esiste, creala
else: map( os.unlink, [os.path.join(cartella_output,f) for f in os.listdir(cartella_output)] ) # Se già esiste, rimuovi il contenuto prima (vedi utilizzo "map" e "list comprehension")


# ===============
# Ciclo temporale
# ===============
for n in range( NMAX ):

    print( 'Iteration: %8d, Time: %8.4f' % (n, time) )

    # Variabili di Comodo
    # -------------------
    Y, q = U
    u = q / Y

    # Output
    # ------
    filename = os.path.join( cartella_output, '%s_%010d.txt' % (problema,n) ) # Nome del file di output per l'iterazione n-esima
    if not n%10: np.savetxt( filename, (x[1:-1], b[1:-1], Y[1:-1], q[1:-1]) ) # Salva solo ogni 10 steps

    # Time control
    # ------------
    if time >= TIMEOUT: break
    ... # calcolare dt!
    if time+dt > TIMEOUT: dt = TIMEOUT-time

    # Condizioni al Contorno
    # ----------------------
    U = CC( time, dt, U, problema )

    # Flussi numerici
    # ---------------
    Flux = FORCE( U, dt, dx )

    # Aggoirna Parte Iperbolica
    # -------------------------
    U[:,1:-1] -= dt/dx * ( Flux[:,1:] - Flux[:,:-1] )

    # Soluzione Termine Sorgente
    U = RK2( Source, U, dt )

    # Update time
    time += dt
    times.append(time)

# Salva la lista dei tempi per avere il riferimento all'indice del ciclo
np.savetxt( timesfile, times )

# ==============
# Grafico Rapido
# ==============
file_list = sorted( os.listdir( cartella_output ) ) # Lista ordinata dei files
times = np.loadtxt( timesfile )
N = len(file_list)
cols = plt.cm.Spectral_r(np.linspace(0,1,N)) # Sequenza di colori
Nplots = 100
plt.figure()
for i in xrange(0, N, int(N/Nplots)):
    f = file_list[i]
    fname = os.path.join( cartella_output, f )
    data = np.loadtxt( fname )
    x, b, Y, q = data[0], data[1], data[2], data[3]
    plt.plot( x, Y, c=cols[i] )
plt.show()
