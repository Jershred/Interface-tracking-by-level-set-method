import numpy as np
import sys
import BC
import scheme
import output

#########################
# Calcul du pas de temps
#########################
def compute_dt(u, v, dx, dy, CFL, N):
    # Variables locales
    uc    = np.zeros(N)
    vc    = np.zeros(N)
    #
    dt = 1.0e15
    for j in range(N):
        uc[:] = 0.5*(u[0:N,j]+u[1:N+1,j])
        dt = min(dt,CFL*dx/np.max(abs(uc)))
    #
    for i in range(N):
        vc[:] = 0.5*(v[i,0:N]+v[i,1:N+1])
        dt = min(dt,CFL*dy/np.max(abs(vc)))
    #
    return dt

###############################
# Calcul de la boucle en temps
###############################
def compute_timeloop(phi, x, y, u, v, dx, dy, CFL, N, itmax, tmax, timeSchemeDef, spaceSchemeDef):
    # Recupere le pas de temps
    dt = compute_dt(u, v, dx, dy, CFL, N)
    #
    # initialiation de la boucle en temps
    continuer = True
    it = 1
    t  = dt
    while (continuer):
        #print("it=",it," t=",t,"s")
        if (spaceSchemeDef == 'UPWIND_O1'):
            compute_timeScheme(phi, x, y, u, v, dx, dy, dt, N, timeSchemeDef, scheme.upwind_O1)
        elif (spaceSchemeDef == 'LAX_WENDROFF'):
            compute_timeScheme(phi, x, y, u, v, dx, dy, dt, N, timeSchemeDef, scheme.lax_wendroff)
        elif (spaceSchemeDef == 'WENO5'):
            compute_timeScheme(phi, x, y, u, v, dx, dy, dt, N, timeSchemeDef, scheme.WENO5)
        else:
            print('Bad value for _spaceSchemeDef='+spaceSchemeDef)
            sys.exit(-1)
        #
        #print('-------------------')
        #
        it = it + 1
        t  = t  + dt
        #
        if (it > itmax):
            continuer = False
            #it_out = it #décommenter si t < tmax
        else:
            continuer = True
            if (t > tmax):
                it_out = it #on stock l'itération pour l'affichage
                it = itmax
                dt = tmax - t + dt
                t  = tmax
    #
    return (it_out, dt)


##################
# Schema en temps
##################
def compute_timeScheme(phi, x, y, u, v, dx, dy, dt, N, timeSchemeDef, spaceScheme):
    if (timeSchemeDef == 'EULER'):
        f1 = np.zeros((N,N))
        #
        spaceScheme(f1, phi, x, y, u, v, dx, dy, dt, N)
        phi[3:N+3,3:N+3] = phi[3:N+3,3:N+3] + f1*dt
        BC.BC_phi(phi, N)
        
    elif (timeSchemeDef == 'RK2'):
        f1   = np.zeros((N,N))
        f2   = np.zeros((N,N))
        phi1 = np.zeros((N+6,N+6))
        #
        spaceScheme(f1, phi, x, y, u, v, dx, dy, dt, N)
        phi1[3:N+3,3:N+3] = phi[3:N+3,3:N+3] + dt*f1
        BC.BC_phi(phi1, N)
        #
        spaceScheme(f2, phi1, x, y, u, v, dx, dy, dt, N)
        phi[3:N+3,3:N+3] =  phi[3:N+3,3:N+3] + 0.5*dt*(f1+f2)
        BC.BC_phi(phi, N)
    
    elif (timeSchemeDef == 'RK4'):
        f1 = np.zeros((N,N))
        f2 = np.zeros((N,N))
        f3 = np.zeros((N,N))
        f4 = np.zeros((N,N))
        phi1 = np.zeros((N+6,N+6))
        phi2 = np.zeros((N+6,N+6))
        phi3 = np.zeros((N+6,N+6))
        #
        spaceScheme(f1, phi, x, y, u, v, dx, dy, dt    , N)
        phi1[3:N+3,3:N+3] = phi[3:N+3,3:N+3] + 0.5*dt*f1
        BC.BC_phi(phi1, N)
        #
        spaceScheme(f2, phi1, x, y, u, v, dx, dy, 0.5*dt, N)
        phi2[3:N+3,3:N+3] = phi[3:N+3,3:N+3] + 0.5*dt*f2
        BC.BC_phi(phi2, N)
        #
        spaceScheme(f3, phi2, x, y, u, v, dx, dy, 0.5*dt, N)
        phi3[3:N+3,3:N+3] = phi[3:N+3,3:N+3] + dt*f3
        BC.BC_phi(phi3, N)
        #
        spaceScheme(f4, phi3, x, y, u, v, dx, dy, dt    , N)
        phi[3:N+3,3:N+3] = phi[3:N+3,3:N+3] + dt/6.0*(f1+2.0*f2+2.0*f3+f4)
        BC.BC_phi(phi, N)
    else:
        print('Bad value for _timeSchemeDef='+timeSchemeDef)
        sys.exit(-1)
    #
    return
