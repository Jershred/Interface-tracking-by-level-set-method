import numpy as np
import sys
import BC
#from numba import jit

###############################
# Initialisation des variables
###############################
def init_variable(u, v, phi, x, y, L, N, testCase):
    # Cas test du serpentin
    if (testCase == 'serpentin'):
       init_variable_serpentin(u, v, phi, x, y, L, N)
    #
    # Cas test de Zalesak
    elif (testCase == 'zalesak'):
        init_variable_zalesak(u, v, phi, x, y, L, N)
    #
    else:
        print('Bad value for _testCase='+testCase)
        sys.exit(-1)
    return



##########################################
# Initialisation du cas test du serpentin
##########################################
def init_variable_serpentin(u, v, phi, x, y, L, N):
    # Initialisation de la fonction couleur phi
    Xcenter = 0.5*L
    Ycenter = 0.75*L
    Radius  = 0.15*L
    #
    xc = np.zeros(N)
    yc = np.zeros(N)
    xc[:] = 0.5*(x[0:N] + x[1:N+1])
    yc[:] = 0.5*(y[0:N] + y[1:N+1])
    #
    for j in range(N):
        for i in range(N):
            phi[i+3,j+3] = np.sqrt((xc[i]-Xcenter)**2+(yc[j]-Ycenter)**2)-Radius
    BC.BC_phi(phi, N)
    #
    # Initialisation de la vitesse
    # Vitesse u
    for j in range(N):
        for i in range(N+1):
            xu = x[i]
            yu = 0.5*(y[j] + y[j+1])
            u[i,j]  = 2.0*np.sin(np.pi*xu)**2*np.sin(np.pi*yu)*np.cos(np.pi*yu)
    #
    # Vitesse v
    for j in range(N+1):
        for i in range(N):
            xv = 0.5*(x[i] + x[i+1])
            yv = y[j]
            v[i,j]  = -2.0*np.sin(np.pi*xv)*np.cos(np.pi*xv)*np.sin(np.pi*yv)**2
    #
    return



##########################################
# Initialisation du cas test de Zalesak
##########################################
def init_variable_zalesak(u, v, phi, x, y, L, N):
    nb_cer  = 10000
    nb_cote = 6000
    nb_haut = 4000
    #
    nb_pts = nb_cer + 2*nb_cote + nb_haut
    lignex = np.zeros(nb_pts)
    ligney = np.zeros(nb_pts)
    ang    = np.arcsin(2.5/15.0)
    #
    xc = np.zeros(N)
    yc = np.zeros(N)
    xc[:] = 0.5*(x[0:N] + x[1:N+1])
    yc[:] = 0.5*(y[0:N] + y[1:N+1])
    #
    it = 0
    for i in range(nb_cer):
        theta  = -0.5*np.pi + ang + 2.0*(np.pi-ang)*float(i)/float(nb_cer-1)
        lignex[it] = 50.0 + 15.0*np.cos(theta)
        ligney[it] = 75.0 + 15.0*np.sin(theta)
        it = it + 1
    #
    print("1")
    for i in range(nb_cote):
        lignex[it] = 47.5
        ligney[it] = 75.0 - 15.0*np.cos(ang) + (10.0+15.0*np.cos(ang))*float(i)/float(nb_cote-1)
        it = it + 1
    #
    for i in range(nb_haut):
        lignex[it] = 47.5 + 5.0*float(i)/float(nb_haut-1)
        ligney[it] = 85.0
        it = it + 1
    #
    for i in range(nb_cote):
        lignex[it] = 52.5
        ligney[it] = 75.0 - 15.0*np.cos(ang) + (10.0+15.0*np.cos(ang))*float(i)/float(nb_cote-1)
        it = it + 1
    
    it = it - 1
    print("2")
    find_signed_distance(phi, xc, yc, lignex, ligney, it, N)
    print("fin")
    BC.BC_phi(phi, N)
    #
    # Vitesse u
    for j in range(N):
        for i in range(N+1):
            yu = 0.5*(y[j] + y[j+1])
            u[i,j] = np.pi/314.0*(50.0-yu)       
    # Vitesse v
    for j in range(N+1):
        for i in range(N):
            xv = 0.5*(x[i] + x[i+1])
            v[i,j] = np.pi/314.0*(xv-50.0)
    #
    return


#@jit #permet de vectorisé les boucles for
################################################
# Calcul de la distance signee (cas de Zalesak)
################################################
def find_signed_distance(phi, xc, yc, lignex, ligney, it, N):
    phi[:,:] = 1.0e15
    for j in range(N):
        for i in range(N):
            for k in range(it):
                phi[i+3,j+3] = min(phi[i+3,j+3], ((xc[i]-lignex[k])**2 + (yc[j]-ligney[k])**2)**0.5)
            #
            if ( ((xc[i]-50.0)**2+(yc[j]-75.0)**2)**0.5  <=  15.0 ):
                if ( (yc[j]  <  85.0) and (xc[i]  >  47.5) and (xc[i]  <  52.5) ):
                    phi[i+3,j+3] = -phi[i+3,j+3]
            else:
                phi[i+3,j+3] = -phi[i+3,j+3]
