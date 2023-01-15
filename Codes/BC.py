import numpy as np

################################
# Calcul des conditions limites
################################
def BC_phi(phi, N):
    # i=0, i=1, i=2 pour tout j
    phi[2,3:N+3] = 2.0*phi[3,3:N+3]-phi[4,3:N+3]
    phi[1,3:N+3] = 2.0*phi[2,3:N+3]-phi[3,3:N+3]
    phi[0,3:N+3] = 2.0*phi[1,3:N+3]-phi[2,3:N+3]
    #
    # i=N+3, i=N+4, i=N+5 pour tout j
    phi[N+3,3:N+3] = 2.0*phi[N+2,3:N+3]-phi[N+1,3:N+3]
    phi[N+4,3:N+3] = 2.0*phi[N+3,3:N+3]-phi[N+2,3:N+3]
    phi[N+5,3:N+3] = 2.0*phi[N+4,3:N+3]-phi[N+3,3:N+3]
    #
    # j=0, j=1, j=2 pour tout i
    phi[3:N+3,2] = 2.0*phi[3:N+3,3] + phi[3:N+3,4]
    phi[3:N+3,1] = 2.0*phi[3:N+3,2] + phi[3:N+3,3]
    phi[3:N+3,0] = 2.0*phi[3:N+3,1] + phi[3:N+3,2]
    #
    # j=N+3, j=N+4, j=N+5 pour tout i
    phi[3:N+3,N+3] = 2.0*phi[3:N+3,N+2] + phi[3:N+3,N+1]
    phi[3:N+3,N+4] = 2.0*phi[3:N+3,N+3] + phi[3:N+3,N+2]
    phi[3:N+3,N+5] = 2.0*phi[3:N+3,N+4] + phi[3:N+3,N+3]
    #
    return
    
