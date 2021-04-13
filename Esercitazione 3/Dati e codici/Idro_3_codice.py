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
#                                                                                #
# Studenti: Federico Panconi                                                     #
#           Riccardo Busti                                                       #
#           Ludovico Agostini                                                    #
#                                                                                #
# ============================================================================== #

#Parametri utilizzati per diga
# L = 100 # Lunghezza del canale [m] 
# IMAX = 10000 # Numero di celle del dominio [-]
# NMAX = 10000 # Numero massimo di iterazioni temporali [-]L = 100 # Lunghezza del canale [m] 
# TIMEOUT = 3600 # Tempo finale della simulazione [s]

#Parametri utilizzati per piena
# L = 100e3 # Lunghezza del canale [m] 
# IMAX = 10000 # Numero di celle del dominio [-]
# NMAX = 100000 # Numero massimo di iterazioni temporali [-]L = 100 # Lunghezza del canale [m] 
# TIMEOUT = 3600e1 # Tempo finale della simulazione [s]

# ======================
# Importazione Pacchetti
# ======================
import numpy as np
import matplotlib
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
L = 300 # Lunghezza del canale [m] 
IMAX = 1000# Numero di celle del dominio [-]
NMAX = 1000# Numero massimo di iterazioni temporali [-]
ks = 30 # Coefficiente di scabrezza del fondo [m^(1/2)/s]
TIMEOUT = 3600 # Tempo finale della simulazione [s]
CFL = 0.9 # Coefficiente di Courant-Friedrichs-Levy  [-]

# Parametri Diga 
# --------------
xdiga = 5 # Posizione dello sbarramento [m]
YL = 1 # Livello a destra dello sbarramento [m]
YR = 1e-03 # Livello a sinistra dello sbarramento [m]

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
        
    # elif problema == 'paratoia':
    #     Y[0], q[0] = Y[1], q[1] # Condizioni di monte
    #     if t >= t_chiusura:
    #         Y[-1], q[-1] = Y[-2], q[-2] # Condizioni di valle
    #     else:
    #         Y[-1], q[-1] = Y[-2], -q[-2] # Condizioni di valle

    elif problema == 'piena':
        if t < t_hydro[0]: # Se la piena deve ancora arrivare...
            qt = q0
            Yt = UniFlow( qt, ks, iF )
        elif t > t_hydro[-1]: # Se la piena è già terminata...
            qt = q0
            Yt = UniFlow( qt, ks, iF )
        else: # Sta transitando l'onda di piena
            qt = np.interp( t, t_hydro, q_hydro )
            Fr = qt/(np.sqrt(g)*UniFlow(qt,ks,iF)**(3/2)) # Numero di Froude
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
                j= qR**2/(ks**2*YR**(10/3))
                Yt = (dt * (g * YR * (iF - j)) + qR - qt) / (-qR / YR - np.sqrt(g * YR)) + YR
        Y[0], q[0] = Yt, qt
        Y[-1], q[-1] = Y[-2], q[-2]  
        
    U[0,:] = Y
    U[1,:] = q
    
    return U

def UniFlow( q, ks, iF ):
    '''Calcola la profondita' di moto uniforme'''
    Y = (q/(ks*np.sqrt(iF)))**(3/5)
    return Y

def PhysFlux( U ):
    '''Calcolo dei flussi fisici'''
    Y, q = U
    Flux = np.array( [q, q**(2)/Y+g*0.5*Y**2] )
    return Flux

def LaxFriedrichs( U, dt, dx ):
    '''Flusso Numerico di Lax-Friedrichs'''
    NumFlux = 0.5 * (PhysFlux(U[:, :-1]) + PhysFlux(U[:, 1:])) - 0.5 * dx / dt * (U[:, 1:] - U[:, :-1])
    return NumFlux

def LaxWendroff( U, dt, dx ):
    '''Flusso Numerico di Lax-Wendroff'''
    NumFlux = PhysFlux (0.5 * (U[:,:-1]+U[:,1:]) - 0.5 * dt/dx * ( (PhysFlux(U[:,1:])-PhysFlux(U[:,:-1])) ))
    return NumFlux 

def FORCE( U, dt, dx ):
    '''Flusso Numerico FORCE di Toro'''
    NumFlux = 0.5 * ( LaxFriedrichs( U, dt, dx ) + LaxWendroff( U, dt, dx ) )
    return NumFlux

def Source( U ):
    '''Termine Sorgente'''
    S=np.array([np.zeros(IMAX+2), g*U[0,:]*(iF-U[1,:]**2/(ks**2*U[0,:]**(10/3)))])
    return S

def RK2(U, dt ):
    '''Schema di Runge-Kutta del secondo ordine'''
    k1= Source(U)
    k2= Source(U + dt*k1)
    Unew= U + dt*(k1+k2)/2
    return Unew

def digaAnalitica(x,time,Y0,xdiga):
    '''
    Soluzione analitica dam-break
    
    Parameters
    ----------
    x : np.array
        distanza [m]
    t : np.array
        tempo [s]
    Y0 : float
        Tirante sx diga [m]
    xdiga : float
        cooridinata x della diga [m]

    Returns
    -------
    Y : np.array
        tiranti

    '''
    Y=np.zeros([len(time), len(x)])
    for i,val in enumerate (time):   
        for j in range(len(x)):
            if (x[j]-xdiga)/val <= -np.sqrt(g*Y0):
                Y[i][j]=Y0
            elif (x[j]-xdiga)/val > -np.sqrt(g*Y0) and (x[j]-xdiga)/val< 2*np.sqrt(g*Y0):
                Y[i][j]=1/(9*g)*(2*np.sqrt(g*Y0) - ((x[j]-xdiga)/val))**2
            elif (x[j]-xdiga)/val >= 2*np.sqrt(g*Y0) :
                Y[i][j]=0
    return Y

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
    Y = np.where( x<xdiga, YL,YR ) # Profondità
    q = np.zeros( IMAX+2 ) # Portata
elif problema == 'paratoia':
    Y = np.ones( IMAX+2 ) * UniFlow( q0, ks, iF )
    q = np.ones( IMAX+2 ) * q0
elif problema == 'piena':
    t_hydro, q_hydro = np.loadtxt( file_idrogramma ).T
    q0 = q_hydro[0]
    Y = np.ones(IMAX + 2) * UniFlow(q0, ks, iF)
    q = np.ones(IMAX + 2) * q0
else:
    # Esci dal codice con un errore
    raise ValueError( '<problema> deve essere "diga", "paratoia" o "piena", non "%s"!' %problema)

U = np.array([Y,q]) # Variabili conservate

time = 0 # Tempo corrente
times = [] # Lista dei tempi da salvare alla fine


# ==================================
# Predisposizione cartella di output
# ==================================
cartella_output = 'output_%s' % problema # Nome cartella di output
timesfile = '%s_times.txt' % problema # Nome file dove salvare i tempi
if not os.path.isdir( cartella_output ): os.mkdir( cartella_output ) # Se la cartella di output non esiste, creala
else: list(map( os.unlink, [os.path.join(cartella_output,f) for f in os.listdir(cartella_output)] )) # Se già esiste, rimuovi tutto il contenuto 


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
    if not n%10: 
        np.savetxt( filename, (x[1:-1], b[1:-1], Y[1:-1], q[1:-1]) ) # Salva solo ogni 10 steps
        times.append(time) #Salvo il tempo corrente
        
    # Time control
    # ------------
    if time >= TIMEOUT: break
    
    # Calcolo dt
    # ----------
    c = np.sqrt(g*Y[1:-1]) # Celerità di propagazione onde di piccola ampiezza
    positiveLambda = u[1:-1]+c # Autovalore positivo
    negativeLambda = u[1:-1]-c # Autovalore negativo
    dt = CFL*dx/max(max(abs(positiveLambda)),max(abs(negativeLambda))) #Passo temporale
    
    if time+dt > TIMEOUT: dt = TIMEOUT-time

    # Condizioni al Contorno
    # ----------------------
    U = CC( time, dt, U, problema )

    # Flussi numerici
    # ---------------
    Flux = FORCE( U, dt, dx )

    # Aggiorna Parte Iperbolica
    # -------------------------
    U[:,1:-1] -= dt/dx * ( Flux[:,1:] - Flux[:,:-1] )

    # Soluzione Termine Sorgente
    # --------------------------
    U = RK2( U, dt )

    # Update time
    # -----------
    time += dt
    
    
    

# Salva la lista dei tempi per avere il riferimento all'indice del ciclo
np.savetxt( timesfile, times )


# ==============
# Plots
# ==============

# Import files
# -------------
file_list = sorted( os.listdir( cartella_output ) ) # Lista ordinata dei files
times = np.loadtxt( timesfile ) # file dei tempi
N = len(file_list) # numero di file
Nplots = 100 # step dei plot N/Nplot nel cicli di stampa

# Font set matplotlib
# ------------------
nice_fonts = {
        'figure.figsize': (14,7 ),
        "text.usetex": True, # necessario installare LaTeX
        "font.family": "serif",
        "axes.labelsize": 20,
        "font.size": 20,
        "legend.fontsize": 16,
        "xtick.labelsize": 16,
        "ytick.labelsize": 16,
}

plt.rcParams.update(nice_fonts)

# Set colorbar 1
# --------------
colors=np.linspace(0,times[-1],N)
norm = matplotlib.colors.Normalize( vmin=np.min(colors), vmax=np.max(colors))
c_m = matplotlib.cm.Spectral
s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
s_m.set_array([])

# Set colorbar 2
# --------------
colors_2 = np.linspace(0,int(100e3),int(100e3+1))
norm_2 = matplotlib.colors.Normalize( vmin=np.min(colors_2), vmax=np.max(colors_2))
c_m_2 = matplotlib.cm.Spectral
s_m_2 = matplotlib.cm.ScalarMappable(cmap=c_m_2, norm=norm_2)
s_m_2.set_array([])


# Piena
# -----
if problema=='piena':   
    
    # Plot idrogramma di piena
    fig=plt.figure()
    plt.plot( t_hydro, q_hydro)
    plt.title('Idrogramma di piena')
    plt.xlabel('Tempo [s]')
    plt.ylabel('Portata specifica [$m^3/s \cdot m^{-1}$]')
    plt.show()
    fig.savefig('img/idrogrammaPiena.png')

    # Onda di piena
    fig, ax1 = plt.subplots()
    for i in range(0, N, int(N/Nplots)):
        f = file_list[i]
        fname = os.path.join( cartella_output, f )
        data = np.loadtxt( fname )
        x, b, Y, q = data[0], data[1], data[2], data[3]
        ax1.plot( x, q, color=s_m.to_rgba(colors[i]))
    ax1.set_xlabel('Distanza [$m$]')
    ax1.set_ylabel('Portata specifica [$m^3/s \cdot m^{-1}$]')
    colorbar=plt.colorbar(s_m)
    colorbar.set_label(label='Tempo [$s$]', rotation=270, labelpad=25)
    plt.show()
    fig.savefig('img/ondaPiena_portata.png')
    
    # Onda di piena_tirante
    fig, ax1 = plt.subplots()
    for i in range(0, N, int(N/Nplots)):
        f = file_list[i]
        fname = os.path.join( cartella_output, f )
        data = np.loadtxt( fname )
        x, b, Y, q = data[0], data[1], data[2], data[3]
        ax1.plot( x, Y, color=s_m.to_rgba(colors[i]))
    ax1.set_xlabel('Distanza [$m$]')
    ax1.set_ylabel('Tirante [$m$]')
    colorbar=plt.colorbar(s_m)
    colorbar.set_label(label='Tempo [$s$]', rotation=270, labelpad=25)
    plt.show()
    fig.savefig('img/ondaPiena_tirante.png')
        
    # Cappio di piena
    dist=[10001, 20002, 30003, 40004, 50005]  #[10km, 20km, 30km,40km, 50km] sezioni in cui valuto il cappio
    colors=['red', 'orange', 'green', 'blue', 'purple']
    scalaQ=np.linspace(1,10)
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
        plt.plot(cappioQ, cappioY, color=colors[i], label="Cappio di piena x=%d Km" %((i+1)*10))
    plt.plot(scalaQ,UniFlow(scalaQ,ks,iF),"--",color="black",markersize=8, label="Scala di Deflusso Moto Uniforme")
    plt.xlabel('Portata specifica [$m^3/s \cdot m^{-1}$]')
    plt.ylabel('Tirante [$m$]')
    plt.legend(loc='best')
    plt.show()
    fig.savefig('img/cappioPiena.png')
    
    # Idrogramma di portata - tempo
    dist=[10001, 20002, 30003, 40004, 50005]  #[10km, 20km, 30km,40km, 50km] sezioni in cui valuto l'idrogramma
    colors=['red', 'orange', 'green', 'blue', 'purple']
    fig=plt.figure()
    Qmax  =[] # Portate massime
    Qmax_time = [] # Tempi dei massimi
    for i,val in enumerate(dist):
        idroT =[] # lista dei tiranti
        idroQ =[] # lista delle portate
        for j in range(0, N, int(N/Nplots)):
            f = file_list[j]
            fname = os.path.join( cartella_output, f )
            data = np.loadtxt( fname )
            x, b, Y, q = data[0], data[1], data[2], data[3]
            x_round=[round(num,1) for num in x]
            index=x_round.index(float(val))
            idroT.append(times[j])
            idroQ.append(q[index])
        Qmax.append(max(idroQ))
        Q_max_index=idroQ.index(max(idroQ))
        Qmax_time.append(idroT[Q_max_index])
        plt.plot(idroT, idroQ, color=colors[i], label="Idrogramma  x=%d Km" %((i+1)*10))
    plt.plot(Qmax_time, Qmax, '--', color='black', label='Inviluppo')
    plt.xlabel('Tempo [$s$]')
    plt.ylabel('Portata specifica [$m^3/s \cdot m^{-1}$]')
    plt.legend(loc='best', fontsize=14)
    plt.show()
    fig.savefig('img/idrogrammaPiena_portata.png')
    
    # Idrogramma di tirante - tempo
    dist=[10001, 20002, 30003, 40004, 50005]  #[10km, 20km, 30km,40km, 50km] sezioni in cui valuto l'idrogramma
    colors=['red', 'orange', 'green', 'blue', 'purple']
    fig=plt.figure()
    Ymax  =[] # Portate massime
    Ymax_time = [] # Tempi dei massimi
    for i,val in enumerate(dist):
        idroT =[] # lista dei tiranti
        idroY =[] # lista delle portate
        for j in range(0, N, int(N/Nplots)):
            f = file_list[j]
            fname = os.path.join( cartella_output, f )
            data = np.loadtxt( fname )
            x, b, Y, q = data[0], data[1], data[2], data[3]
            x_round=[round(num,1) for num in x]
            index=x_round.index(float(val))
            idroT.append(times[j])
            idroY.append(Y[index])
        Ymax.append(max(idroY))
        Y_max_index=idroY.index(max(idroY))
        Ymax_time.append(idroT[Y_max_index])
        plt.plot(idroT, idroY, color=colors[i], label="Idrogramma  x=%d Km" %((i+1)*10))
    plt.plot(Ymax_time, Ymax, '--', color='black', label='Inviluppo')
    plt.xlabel('Tempo [$s$]')
    plt.ylabel('Tirante [$m$]')
    plt.legend(loc='best', fontsize=14)
    plt.show()
    fig.savefig('img/idrogrammaPiena_tirante.png')
    
   
# Diga
# ----
elif problema=='diga':
    
    #Profilo dell'onda
    fig, ax1 = plt.subplots()
    for i in range(0, N, int(N/Nplots)):
        f = file_list[i]
        fname = os.path.join( cartella_output, f )
        data = np.loadtxt( fname )
        x, b, Y, q = data[0], data[1], data[2], data[3]
        ax1.plot( x, Y, color=s_m.to_rgba(colors[i]))
    ax1.set_xlabel('Distanza [$m$]')
    ax1.set_ylabel('Tirante [$m$]')
    colorbar=plt.colorbar(s_m)
    colorbar.set_label(label='Tempo [$s$]', rotation=270, labelpad=25)
    plt.show()
    fig.savefig('img/profiloDiga_num.png')
   
    #Profilo dell'onda vs 4/9Y0
    fig, ax1 = plt.subplots()
    x_line=np.linspace(0,int(L-1),int(L))
    y_line=np.ones(int(L))*(4/9*YL)
    for i in range(0, N, int(N/Nplots)):
        f = file_list[i]
        fname = os.path.join( cartella_output, f )
        data = np.loadtxt( fname )
        x, b, Y, q = data[0], data[1], data[2], data[3]
        ax1.plot( x, Y, color=s_m.to_rgba(colors[i]))
    ax1.plot( x_line, y_line, '--k', label='$\\frac{4}{9}Y_0$')
    ax1.set_xlabel('Distanza [$m$]')
    ax1.set_ylabel('Tirante [$m$]')
    ax1.legend(loc='best')
    colorbar=plt.colorbar(s_m)
    colorbar.set_label(label='Tempo [$s$]', rotation=270, labelpad=25)
    plt.show()
    fig.savefig('img/profiloDiga_num_4_9_Y0.png')
    
    # Soluzione analitica
    fig, ax1 = plt.subplots()
    x_solution=np.linspace(0,int(L),10000)
    digaSolution=digaAnalitica(x_solution,times,YL,xdiga)
    x_line=np.linspace(0,int(L-1),int(L))
    y_line=np.ones(int(L))*(4/9*YL)
    for i in range(0, N, int(N/Nplots)):
        ax1.plot( x_solution,digaSolution[i,:], color=s_m.to_rgba(colors[i]))
    ax1.plot( x_line, y_line, '--k', label='$\\frac{4}{9}Y_0$')
    ax1.set_xlabel('Distanza [$m$]')
    ax1.set_ylabel('Tirante [$m$]')
    ax1.legend(loc='best')
    colorbar=plt.colorbar(s_m)
    colorbar.set_label(label='Tempo [$s$]', rotation=270, labelpad=25)
    plt.show()
    fig.savefig('img/profiloDiga_anal_4_9_Y0.png')
    
    # Confronto numerica-analitica
    fig, ax1 = plt.subplots()
    x_solution=np.linspace(0,int(L),10000)
    digaSolution=digaAnalitica(x_solution,times,YL,xdiga)
    x_line=np.linspace(0,int(L-1),int(L))
    y_line=np.ones(int(L))*(4/9*YL)
    colors=['red', 'orange', 'green', 'blue', 'purple']
    index=0
    for i in range(0, 100,20) :
        f = file_list[i]
        fname = os.path.join( cartella_output, f )
        data = np.loadtxt( fname )
        x, b, Y, q = data[0], data[1], data[2], data[3]
        ax1.plot(x, Y, label="Numerica t=%1.2fs" %times[i] , color=colors[index])
        ax1.plot( x_solution,digaSolution[i,:],'--', label="Analitica t=%1.2fs" %times[i], color=colors[index])
        index+=1
    ax1.plot( x_line, y_line, ':k', label='$\\frac{4}{9}Y_0$')
    ax1.set_xlim(0,10)
    ax1.set_xlabel('Distanza [$m$]')
    ax1.set_ylabel('Tirante [$m$]')
    #plt.legend(loc='upper right', bbox_to_anchor=(1.25,1), ncol=1)
    plt.legend(loc='best', fontsize=14)
    plt.show()
    fig.savefig('img/confronto_num_anal.png')
