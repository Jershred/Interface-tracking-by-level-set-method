import sys
import numpy as np
#import tool
import init
import output
import timeloop
import time


#########################
# Principales constantes
#########################
_testCase = 'serpentin'                           # Cas test du serpentin
#_testCase = 'zalesak'                              # Cas test disque tournant
if (_testCase == 'serpentin'):
    _L = 1.0                                       # Longueur du domaine
elif (_testCase == 'zalesak'):
    _L = 100.0                                     # Longueur du domaine
else:
    print('Bad value for _testCase='+_testCase)
    sys.exit(-1)
#
_N = 512                                           # Nombre de cellules par direction
if (_testCase == 'serpentin'):
    _tmax = 3.0                                    # Temps de simulation max (serpentin)
elif (_testCase == 'zalesak'):
    _tmax = 628.0                                  # Temps de simulation max (Zalesak)
else:
    print('Bad value for _testCase='+_testCase)
    sys.exit(-1)
#
#_itmax = 100000                                    # Nombre d'iterations temporelles max
_itmax = 100000                               # Nombre d'iterations temporelles max
_CFL   = 0.3                                       # CFL number
#
# Schema temporel
#_timeSchemeDef = 'EULER'                          # Euler ordre1 
#_timeSchemeDef = 'RK2'                            # Runge-Kutta ordre 2
#_timeSchemeDef = 'RK4'                             # Runge-Kutta ordre 4
#
# Schema spatial
#_spaceSchemeDef = 'UPWIND_O1'               # Decentre amonte (upwind) ordre 1
#_spaceSchemeDef = 'LAX_WENDROFF'           # Schema de Lax-Wendroff
#_spaceSchemeDef = 'WENO5'                  # Decentre WENO ordre 5

Time = ['EULER','RK2','RK4']                    #Vecteur temporel pour automatisation
Space = ['LAX_WENDROFF','UPWIND_O1','WENO5']    #Vecteur saptial pour automatisation

######################
# Programme principal
######################
# Initialisation des variables
dx   = _L/_N                                # Pas d'espace suivant x
dy   = _L/_N                                # Pas d'espace suivant y
x    = np.arange(0.0, _L+dx, dx)            # Abscisse des noeuds du maillage (taille N+1)
y    = np.arange(0.0, _L+dy, dy)            # Ordonnee des noeuds du maillage (taille N+1)
phi0 = np.zeros((_N+6, _N+6))               # Fonction level-set (champ initial)
phi  = np.zeros((_N+6, _N+6))               # Fonction level-set
u    = np.zeros((_N+1, _N))                 # Vecteur vitesse (composante suivant x)
v    = np.zeros((_N  , _N+1))               # Vecteur vitesse (composante suivant y)

###############################
# Initialisation des variables
###############################

for _spaceSchemeDef in Space :
    for _timeSchemeDef in Time :
    
        print('-------------------')
        print('Running :' + '{}'.format(_timeSchemeDef +' / '+ _spaceSchemeDef))
        
        start = time.time()
        
        init.init_variable(u, v, phi, x, y, _L, _N, _testCase)
        phi0[:,:] = phi[:,:]
        
        ##################
        # Boucle en temps
        ##################
        it_out, dt = timeloop.compute_timeloop(phi, x, y, u, v, dx, dy, _CFL, _N, _itmax, _tmax, _timeSchemeDef, _spaceSchemeDef)
        print('error=', output.compute_errors(phi0, phi, _N))
        print('t=', dt*it_out)
        output.write_output(x, y, phi0, phi, u, v, _N, _timeSchemeDef, _spaceSchemeDef, it_out, _testCase)
        
        end = time.time()
        
        h = int((end-start)/3600) #heure
        m = int((end-start)/60)-h*60 #minute
        s = int((end-start)-h*3600-m*60) #heure
        print("Duration:",h,"h",m,"min",s,"s")
        
        print('End')