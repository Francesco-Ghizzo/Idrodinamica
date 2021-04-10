def UniFlow( q, ks, iF ):
    '''Calcola la profondita' di moto uniforme'''
    Y = (q/(ks*np.sqrt(iF)))**(3/5)
    return Y

def PhysFlux( U ):
    '''Calcolo dei flussi fisici'''
    Y, q = U
    Flux = np.array([q, q**2/Y+0.5*g*Y**2])
    return Flux

def LaxFriedrichs( U, dt, dx ):
    '''Flusso Numerico di Lax-Friedrichs'''
    NumFlux = 0.5*(PhysFlux(U[:,:-1]))+PhysFlux(U[:,:])-0.5*(dx/dt)*(U[:,:]-U[:,:-1])
    return NumFlux

def LaxWendroff( U, dt, dx ):
    '''Flusso Numerico di Lax-Wendroff'''
    U = 0.5*(U[:,:]+U[:,+1:])+0.5*(dt/dx)*(PhysFlux(U[:,:])-PhysFlux(U[:,+1:]))
    NumFlux = PhysFlux(U)
    return NumFlux

def FORCE( U, dt, dx ):
    '''Flusso Numerico FORCE di Toro'''
    NumFlux = 0.5(LaxFriedrichs( U, dt, dx )+LaxWendroff( U, dt, dx ))
    return NumFlux

def Source( U ):
    '''Termine Sorgente'''
    S = np.array([0, g*U[0]*(iF-j)])
    return S

def RK2( S, U, dt ):
    '''Schema di Runge-Kutta del secondo ordine'''
    K1 = dt*S[:,:]
    K2 = dt*Source(U+K1)
    Unknown = U + 0.5*(K1+K2)
    return Unknown
