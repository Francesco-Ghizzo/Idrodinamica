def UniFlow( q, ks, iF ):
    '''Calcola la profondita' di moto uniforme'''
    Y = (q/(ks*np.sqrt(iF)))**(3/5)
    return Y

def PhysFlux( U ):
    '''Calcolo dei flussi fisici'''
    Y, q = U
    Flux = ...
    return Flux

def LaxFriedrichs( U, dt, dx ):
    '''Flusso Numerico di Lax-Friedrichs'''
    ...
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
