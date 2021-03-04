def UniFlow( q, ks, iF ):
    '''Calcola la profondita' di moto uniforme'''
    Y = (q/(ks*np.sqrt(iF)))**(3/5)
    return Y

def PhysFlux( U ):
    '''Calcolo dei flussi fisici'''
    Y, q = U
    F = np.array([q, q**2/Y+0.5*g*Y**2])
    return F

def LaxFriedrichs( U, dt, dx ):
    '''Flusso Numerico di Lax-Friedrichs'''
    F_LF = 0.5*(PhysFlux(U[]))+PhysFlux(U[]))-0.5*(dx/dt)*(U[]-U[])
    return F_LF

def LaxWendroff( U, dt, dx ):
    '''Flusso Numerico di Lax-Wendroff'''
    ...
    return F_LW

def FORCE( U, dt, dx ):
    '''Flusso Numerico FORCE di Toro'''
    ...
    return F_TO

def Source( U ):
    '''Termine Sorgente'''
    ...
    return S

def RK2( S, U, dt ):
    '''Schema di Runge-Kutta del secondo ordine'''
    ...
    return ...
