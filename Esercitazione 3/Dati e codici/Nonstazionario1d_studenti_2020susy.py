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
#              diverse combinazioni di condizioni iniziali e al contorno.        #
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
problema = 'diga' # Opzioni: 'diga', 'piena'
g = 9.81 # Accelerazione di gravita' [m/s^2]
iF = 0.001 # Pendenza costante del canale [-]
L = 100e3 # Lunghezza del canale [m]
IMAX =200 # Numero di celle del dominio [-]
NMAX =200# Numero massimo di iterazioni temporali [-]
ks = 30 # Coefficiente di scabrezza del fondo [m^(1/3)/s]
TIMEOUT = 3600e1 # Tempo finale della simulazione [s]
CFL = 0.9# Coefficiente di Courant-Friedrichs-Levy # [-]
q0 = 1 # Portata unitaria di base [m^2/s]



# Parametri Diga
# --------------
xdiga = 10 # Posizione dello sbarramento [m]
YL = 1 # Livello a sinistra dello sbarramento [m]
YR = 1e-03 # Livello a destra dello sbarramento [m]




# Parametri Piena
# ---------------
file_idrogramma = 'idrogramma.txt' # File contenente due colonne x tempo e portata dell'idrogramma


# ========
# Funzioni
# ========

    
def CC( t, dt, U, problema ):
    '''Applicazione delle condizioni al contorno'''
    Y, q = np.float128(U)
    
    if problema == 'diga':
        Y[0], q[0] = Y[1], q[1] # Condizioni di monte
        Y[-1], q[-1] = Y[-2], q[-2]  # Condizioni di valle
        
    elif problema == 'piena':
        if t < t_hydro[0]: # Se la piena deve ancora arrivare...
            qt = q0
            Yt = UniFlow( qt, ks, iF )
        elif t > t_hydro[-1]: # Se la piena è già terminata...
            qt = q0
            Yt = UniFlow( qt, ks, iF )
        else: # Sta transitando l'onda di piena
            qt = np.interp( t, t_hydro, q_hydro ) #interpolazione lineare piecewise
            Fr = qt / (UniFlow(qt, ks, iF) ** (3 / 2) * np.sqrt(g)) # Numero di Froude
            if Fr>=1: # Il controllo è a monte
                Yt = UniFlow(qt, ks, iF)
            else: # Il controllo è a valle
                # Cerco il piede della caratteristica partendo dalla cella 0
                # ed andando a valle finchè non lo trovo.
                # L'iterazione non dovrebbe superare le prime 2-3 celle
                for i in range( x.size ):
                    c= np.sqrt(g*UniFlow(qt, ks, iF))
                    u= qt/UniFlow(qt, ks, iF)
                    xR = x[0]- dt*(u-c)
                    if x[i]-dx/2 <= xR <= x[i]+dx/2: break # Ho trovato xR!
                # Uso dei valori interpolati di Y e q
                YR = np.interp( xR, x, Y )
                qR = np.interp( xR, x, q )
                j= qR ** 2 / (YR ** (10 / 3) * ks ** 2)
                Yt = (dt * (g * YR * (iF - j)) + qR - qt) / (-qR / YR - np.sqrt(g * YR)) + YR
        Y[0], q[0] = Yt, qt
        Y[-1], q[-1] = Y[-2], q[-2]  #trasmissiva
        
    U[0,:] = Y
    U[1,:] = q
    
    return U

def UniFlow( q, ks, iF ):
    '''Calcola la profondita' di moto uniforme'''
    Y = (q / (ks * np.sqrt(iF))) ** (3 / 5)
    return Y

def PhysFlux( U ):
    '''Calcolo dei flussi fisici'''
    Y, q = U
    Flux = np.array([q, q ** 2 / Y + g * Y ** 2 / 2])
    return Flux

def LaxFriedrichs( U, dt, dx ):
    '''Flusso Numerico di Lax-Friedrichs'''
    Y, q = U
    NumFlux = 0.5 * (PhysFlux(U[:, :-1]) + PhysFlux(U[:, 1:])) - 0.5 * dx / dt * (U[:, 1:] - U[:, :-1])
    return NumFlux

def LaxWendroff( U, dt, dx ):
    '''Flusso Numerico di Lax-Wendroff'''
    Y, q = U
    NumFlux = PhysFlux (0.5*(U[:,:-1]+U[:,1:])-0.5*dt/(dx*(PhysFlux(U[:,1:])-PhysFlux(U[:,:-1]))))
    return NumFlux

def FORCE( U, dt, dx ):
    '''Flusso Numerico FORCE di Toro'''
    NumFlux = 0.5*(LaxFriedrichs( U, dt, dx )+LaxWendroff( U, dt, dx ))
    return NumFlux

def Source( U ):
    '''Termine Sorgente'''
    S=np.array([np.zeros(IMAX+2), g*U[0,:]*(iF-U[1,:]**2/(U[0,:]**(10/3)*ks**2))])
    return S

def RK2( S, U, dt ):
    '''Schema di Runge-Kutta del secondo ordine'''
    K1= Source( U )
    K2= Source( U + dt*K1 )
    U1= U + dt*(K1+K2)/2
    return U1

#def digaAnalitica(x,time, YL, xdiga):
#    '''
#    Soluzione analitica dam-breack
#    
#    Parameters
#    ---------------
#    x: np.array 
#        distanza [m] 
#    t: np.array  
#       tempo[s]
#    YL: float
#        tirante sx diga[m]
#    xdiga: float
#        coordinata x della [m]
#        
#    Returns
#    --------------
#    Y: np.array
#        tiranti
#        
#   '''
#    
#    Y=np.zeros([len(time), len(x)])
#    for i in range(len(time)):
#        for j in range(len(x)):
#            if (x[j]-xdiga)/time[i] < -np.sqrt(g*YL):
#                Y[i][j]=YL
#            elif(x[j]-xdiga)/time[i]> -np.sqrt(g*YL) and (x[j]-xdiga)/time[i]< 2*np.sqrt(g*YL):
#                Y[i][j]=1/(9*g)*(2*np.sqrt(g*YL)-((x[j]-xdiga)/time[i]))**2
#            elif (x[j]-xdiga)/time[1]>2*np.sqrt(g*YL) :
#                Y[i][j]=0
#                
#    return Y


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
    Y = np.where( x<xdiga, YL,YR) # Profondità
    q = np.zeros( IMAX+2 ) # Portata
 
    
elif problema == 'piena':
    t_hydro, q_hydro = np.loadtxt( file_idrogramma ).T
    q0 = q_hydro[0]
    Y = np.ones(IMAX + 2) * UniFlow(q0, ks, iF)
    q = np.ones(IMAX + 2) * q0
else:
    # Esci dal codice con un errore
    raise ValueError #'<problema> deve essere "diga" o "piena", non "%s"!' 

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

    c= np.sqrt(g*Y[1:-1])
    lambda_neg= u[1:-1]-c
    lambda_pos= u[1:-1]+c
    dt= CFL*dx/max(max(abs(lambda_neg)),max(abs(lambda_pos))) # calcolare dt!
    
    
    if time+dt > TIMEOUT: dt = TIMEOUT-time

    # Condizioni al Contorno
    # ----------------------
    U = CC( time, dt, U, problema )

    # Flussi numerici
    # ---------------
    Flux = LaxFriedrichs( U, dt, dx )

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
if problema == 'diga':
    file_list = sorted( os.listdir( cartella_output ) ) # Lista ordinata dei files
    times = np.loadtxt( timesfile )
    N = len(file_list)
    cols = plt.cm.Spectral_r(np.linspace(0,1,N)) # Sequenza di colori
    Nplots = 10
    plt.figure()
    plt.title("Portata diga")
    plt.xlabel("x [m]")
    plt.ylabel("q [m^(1/3)/(s m)]")

    #x_line=np.linspace(0,int(L-1),int(L))
    #y_line=np.ones(int(L))*(4/9*YL)
    
    for i in range(0, N, int(N/Nplots)):
        f = file_list[i]
        fname = os.path.join( cartella_output, f )
        data = np.loadtxt( fname )
        x, b, Y, q = data[0], data[1], data[2], data[3]
        plt.plot( x, q, c=cols[i] )   #q, Y
    
    #plt.plot(x_line,y_line,'-k',label='$\\frac{4}{9}Y_0$')    
    #plt.legend()



if problema == 'piena':
    file_list = sorted( os.listdir( cartella_output ) ) # Lista ordinata dei files
    times = np.loadtxt( timesfile )
    N = len(file_list)
    cols = plt.cm.Spectral_r(np.linspace(0,1,N)) # Sequenza di colori
    Nplots = 100
    plt.figure()
    plt.title("Portata dell'onda di piena")
    plt.xlabel("x [m]")
    plt.ylabel("q [m^(1/3)/(s m)]")
    
    
    for i in range(0, N, int(N/Nplots)):
        f = file_list[i]
        fname = os.path.join( cartella_output, f )
        data = np.loadtxt( fname )
        x, b, Y, q = data[0], data[1], data[2], data[3]
        plt.plot( x, q, c=cols[i] )   #q, Y





    # Cappio di piena
    dist=[10001, 20002, 30003, 40004, 50005]  #[10km, 20km, 30km,40km, 50km] sezioni in cui valuto il cappio
    colors=['palevioletred', 'cornflowerblue', 'mediumseagreen', 'coral', 'firebrick']
    fig=plt.figure(figsize=(12,12))
    for i,val in enumerate(dist):
        cappioY=[] # lista dei tiranti
        cappioQ=[] # lista delle portate
        for j in range(0, N, int(N/Nplots)):
            f = file_list[j]
            fname = os.path.join( cartella_output, f )
            data = np.loadtxt( fname )
            x, b, Y, q = data[0], data[1], data[2], data[3]
            x_round=[round(num,1) for num in x]
            index=x_round.index(float(val))
            cappioY.append(Y[index])
            cappioQ.append(q[index])
        plt.plot(cappioQ, cappioY, color=colors[i], label="Cappio di piena %d Km" %((i+1)*10))
    plt.xlabel('q [m^3/s m^(-1)]')
    plt.ylabel('Y[m]')
    plt.legend()


plt.show()










