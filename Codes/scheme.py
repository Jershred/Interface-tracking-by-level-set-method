import numpy as np
import sys

#######################
# Schema spatial WENO5
#######################
def WENO5(f, phi, x, y, u, v, dx, dy, dt, N):
    # Variables locales
    uhalf  = np.zeros(N)
    vhalf  = np.zeros(N)
    qLess1 = np.zeros(N)
    qLess2 = np.zeros(N)
    qLess3 = np.zeros(N)
    qLess4 = np.zeros(N)
    qLess5 = np.zeros(N)
    qPlus1 = np.zeros(N)
    qPlus2 = np.zeros(N)
    qPlus3 = np.zeros(N)
    qPlus4 = np.zeros(N)
    qPlus5 = np.zeros(N)
    DxLess0 = np.zeros(N)
    DxLess1 = np.zeros(N)
    DxLess2 = np.zeros(N)
    DxPlus0 = np.zeros(N)
    DxPlus1 = np.zeros(N)
    DxPlus2 = np.zeros(N)
    DyLess0 = np.zeros(N)
    DyLess1 = np.zeros(N)
    DyLess2 = np.zeros(N)
    DyPlus0 = np.zeros(N)
    DyPlus1 = np.zeros(N)
    DyPlus2 = np.zeros(N)
    ISLess0 = np.zeros(N)
    ISLess1 = np.zeros(N)
    ISLess2 = np.zeros(N)
    ISPlus0 = np.zeros(N)
    ISPlus1 = np.zeros(N)
    ISPlus2 = np.zeros(N)
    omegaLess0 = np.zeros(N)
    omegaLess1 = np.zeros(N)
    omegaLess2 = np.zeros(N)
    omegaPlus0 = np.zeros(N)
    omegaPlus1 = np.zeros(N)
    omegaPlus2 = np.zeros(N)
    DxLess = np.zeros(N)
    DxPlus = np.zeros(N)
    DyLess = np.zeros(N)
    DyPlus = np.zeros(N)
    #
    eps = 1.0e-6
    #
    f[:,:] = 0.0
    # Advection suivant i
    for j in range(N):
        uhalf[:]  = 0.5*(u[0:N,j] + u[1:N+1,j])
        #        
        qLess1[:] = (phi[1:N+1,j+3] - phi[0:N  ,j+3])/dx
        qLess2[:] = (phi[2:N+2,j+3] - phi[1:N+1,j+3])/dx
        qLess3[:] = (phi[3:N+3,j+3] - phi[2:N+2,j+3])/dx
        qLess4[:] = (phi[4:N+4,j+3] - phi[3:N+3,j+3])/dx
        qLess5[:] = (phi[5:N+5,j+3] - phi[4:N+4,j+3])/dx
        #
        qPlus1[:] = (phi[6:N+6,j+3] - phi[5:N+5,j+3])/dx
        qPlus2[:] = (phi[5:N+5,j+3] - phi[4:N+4,j+3])/dx
        qPlus3[:] = (phi[4:N+4,j+3] - phi[3:N+3,j+3])/dx
        qPlus4[:] = (phi[3:N+3,j+3] - phi[2:N+2,j+3])/dx
        qPlus5[:] = (phi[2:N+2,j+3] - phi[1:N+1,j+3])/dx
        #
        DxLess0[:] =  (1.0/3.0)*qLess1[:] - (7.0/6.0)*qLess2[:] + (11.0/6.0)*qLess3[:]
        DxLess1[:] = -(1.0/6.0)*qLess2[:] + (5.0/6.0)*qLess3[:] + (1.0/3.0)*qLess4[:]
        DxLess2[:] =  (1.0/3.0)*qLess3[:] + (5.0/6.0)*qLess4[:] - (1.0/6.0)*qLess5[:]
        #
        DxPlus0[:] =  (1.0/3.0)*qPlus1[:] - (7.0/6.0)*qPlus2[:] + (11.0/6.0)*qPlus3[:]
        DxPlus1[:] = -(1.0/6.0)*qPlus2[:] + (5.0/6.0)*qPlus3[:] + (1.0/3.0)*qPlus4[:]
        DxPlus2[:] =  (1.0/3.0)*qPlus3[:] + (5.0/6.0)*qPlus4[:] - (1.0/6.0)*qPlus5[:]
        #
        ISLess0[:] = 13.0*(qLess1[:] - 2.0*qLess2[:] + qLess3[:])**2 + 3.0*(qLess1[:] - 4.0*qLess2[:] + 3.0*qLess3[:])**2 + eps
        ISLess1[:] = 13.0*(qLess2[:] - 2.0*qLess3[:] + qLess4[:])**2 + 3.0*(qLess4[:] - qLess2[:])**2 + eps
        ISLess2[:] = 13.0*(qLess3[:] - 2.0*qLess4[:] + qLess5[:])**2 + 3.0*(3.0*qLess3[:] - 4.0*qLess4[:] + qLess5[:])**2 + eps
        #
        ISPlus0[:] = 13.0*(qPlus1[:] - 2.0*qPlus2[:] + qPlus3[:])**2 + 3.0*(qPlus1[:] - 4.0*qPlus2[:] + 3.0*qPlus3[:])**2 + eps
        ISPlus1[:] = 13.0*(qPlus2[:] - 2.0*qPlus3[:] + qPlus4[:])**2 + 3.0*(qPlus4[:] - qPlus2[:])**2 + eps
        ISPlus2[:] = 13.0*(qPlus3[:] - 2.0*qPlus4[:] + qPlus5[:])**2 + 3.0*(3.0*qPlus3[:] - 4.0*qPlus4[:] + qPlus5[:])**2 + eps
        #
        omegaLess0[:] = ISLess0[:]*ISLess0[:]
        omegaLess1[:] = (1.0/6.0)*(ISLess1[:]*ISLess1[:])
        omegaLess2[:] = (1.0/3.0)*(ISLess2[:]*ISLess2[:])
        #
        omegaPlus0[:] = ISPlus0[:]*ISPlus0[:]
        omegaPlus1[:] = (1.0/6.0)*(ISPlus1[:]*ISPlus1[:])
        omegaPlus2[:] = (1.0/3.0)*(ISPlus2[:]*ISPlus2[:])
        #
        DxLess[:] = (omegaLess1[:]*omegaLess2[:]*DxLess0[:] + omegaLess0[:]*omegaLess2[:]*DxLess1[:] + omegaLess0[:]*omegaLess1[:]*DxLess2[:])/ \
            (omegaLess1[:]*omegaLess2[:] + omegaLess0[:]*omegaLess2[:] + omegaLess0[:]*omegaLess1[:])
        DxPlus[:] = (omegaPlus1[:]*omegaPlus2[:]*DxPlus0[:] + omegaPlus0[:]*omegaPlus2[:]*DxPlus1[:] + omegaPlus0[:]*omegaPlus1[:]*DxPlus2[:])/ \
            (omegaPlus1[:]*omegaPlus2[:] + omegaPlus0[:]*omegaPlus2[:] + omegaPlus0[:]*omegaPlus1[:])
        #
        f[:,j] = f[:,j] - np.maximum(uhalf[:], 0.0)*DxLess[:] - np.minimum(uhalf[:], 0.0)*DxPlus[:]

    # Advection suivant j
    for i in range(N):
        vhalf[:]  = 0.5*(v[i,0:N] + v[i,1:N+1])
        #
        qLess1[:] = (phi[i+3,1:N+1] - phi[i+3,0:N  ])/dy
        qLess2[:] = (phi[i+3,2:N+2] - phi[i+3,1:N+1])/dy
        qLess3[:] = (phi[i+3,3:N+3] - phi[i+3,2:N+2])/dy
        qLess4[:] = (phi[i+3,4:N+4] - phi[i+3,3:N+3])/dy
        qLess5[:] = (phi[i+3,5:N+5] - phi[i+3,4:N+4])/dy
        #
        qPlus1[:] = (phi[i+3,6:N+6] - phi[i+3,5:N+5])/dy
        qPlus2[:] = (phi[i+3,5:N+5] - phi[i+3,4:N+4])/dy
        qPlus3[:] = (phi[i+3,4:N+4] - phi[i+3,3:N+3])/dy
        qPlus4[:] = (phi[i+3,3:N+3] - phi[i+3,2:N+2])/dy
        qPlus5[:] = (phi[i+3,2:N+2] - phi[i+3,1:N+1])/dy
        #
        DyLess0[:] =  (1.0/3.0)*qLess1[:] - (7.0/6.0)*qLess2[:] + (11.0/6.0)*qLess3[:]
        DyLess1[:] = -(1.0/6.0)*qLess2[:] + (5.0/6.0)*qLess3[:] + (1.0/3.0)*qLess4[:]
        DyLess2[:] =  (1.0/3.0)*qLess3[:] + (5.0/6.0)*qLess4[:] - (1.0/6.0)*qLess5[:]
        #
        DyPlus0[:] =  (1.0/3.0)*qPlus1[:] - (7.0/6.0)*qPlus2[:] + (11.0/6.0)*qPlus3[:]
        DyPlus1[:] = -(1.0/6.0)*qPlus2[:] + (5.0/6.0)*qPlus3[:] + (1.0/3.0)*qPlus4[:]
        DyPlus2[:] =  (1.0/3.0)*qPlus3[:] + (5.0/6.0)*qPlus4[:] - (1.0/6.0)*qPlus5[:]
        #
        ISLess0[:] = 13.0*(qLess1[:] - 2.0*qLess2[:] + qLess3[:])**2 + 3.0*(qLess1[:] - 4.0*qLess2[:] + 3.0*qLess3[:])**2 + eps
        ISLess1[:] = 13.0*(qLess2[:] - 2.0*qLess3[:] + qLess4[:])**2 + 3.0*(qLess4[:] - qLess2[:])**2 + eps
        ISLess2[:] = 13.0*(qLess3[:] - 2.0*qLess4[:] + qLess5[:])**2 + 3.0*(3.0*qLess3[:] - 4.0*qLess4[:] + qLess5[:])**2 + eps
        #
        ISPlus0[:] = 13.0*(qPlus1[:] - 2.0*qPlus2[:] + qPlus3[:])**2 + 3.0*(qPlus1[:] - 4.0*qPlus2[:] + 3.0*qPlus3[:])**2 + eps
        ISPlus1[:] = 13.0*(qPlus2[:] - 2.0*qPlus3[:] + qPlus4[:])**2 + 3.0*(qPlus4[:] - qPlus2[:])**2 + eps
        ISPlus2[:] = 13.0*(qPlus3[:] - 2.0*qPlus4[:] + qPlus5[:])**2 + 3.0*(3.0*qPlus3[:] - 4.0*qPlus4[:] + qPlus5[:])**2 + eps
        #
        omegaLess0[:] = ISLess0[:]*ISLess0[:]
        omegaLess1[:] = (1.0/6.0)*(ISLess1[:]*ISLess1[:])
        omegaLess2[:] = (1.0/3.0)*(ISLess2[:]*ISLess2[:])
        #
        omegaPlus0[:] = ISPlus0[:]*ISPlus0[:]
        omegaPlus1[:] = (1.0/6.0)*(ISPlus1[:]*ISPlus1[:])
        omegaPlus2[:] = (1.0/3.0)*(ISPlus2[:]*ISPlus2[:])
        #
        DyLess[:] = (omegaLess1[:]*omegaLess2[:]*DyLess0[:] + omegaLess0[:]*omegaLess2[:]*DyLess1[:] + omegaLess0[:]*omegaLess1[:]*DyLess2[:])/ \
            (omegaLess1[:]*omegaLess2[:] + omegaLess0[:]*omegaLess2[:] + omegaLess0[:]*omegaLess1[:])
        DyPlus[:] = (omegaPlus1[:]*omegaPlus2[:]*DyPlus0[:] + omegaPlus0[:]*omegaPlus2[:]*DyPlus1[:] + omegaPlus0[:]*omegaPlus1[:]*DyPlus2[:])/ \
            (omegaPlus1[:]*omegaPlus2[:] + omegaPlus0[:]*omegaPlus2[:] + omegaPlus0[:]*omegaPlus1[:])
        #
        f[i,:] = f[i,:] - np.maximum(vhalf[:], 0.0)*DyLess[:] - np.minimum(vhalf[:], 0.0)*DyPlus[:]
    #
    return



################################
# Schema spatial upwind ordre 1
################################
def upwind_O1(f, phi, x, y, u, v, dx, dy, dt, N):
    # Variables locales
    DxLess = np.zeros(N)
    DxPlus = np.zeros(N)
    DyLess = np.zeros(N)
    DyPlus = np.zeros(N)
    uhalf  = np.zeros(N)
    vhalf  = np.zeros(N)
    #
    f[:,:] = 0.0
    # Advection suivant i
    for j in range(N):
        uhalf[:]   = 0.5*(u[0:N,j]+u[1:N+1,j])
        DxLess[:]  = (phi[3:N+3,j+3]-phi[2:N+2,j+3])/dx
        DxPlus[:]  = (phi[4:N+4,j+3]-phi[3:N+3,j+3])/dx
        f[:,j] = f[:,j] - np.maximum(uhalf[:],0.0)*DxLess[:] - np.minimum(uhalf[:],0.0)*DxPlus[:]
    #
    # Advection suivant j
    for i in range(N):
        vhalf[:]   = 0.5*(v[i,0:N]+v[i,1:N+1])
        DyLess[:]  = (phi[i+3,3:N+3]-phi[i+3,2:N+2])/dy
        DyPlus[:]  = (phi[i+3,4:N+4]-phi[i+3,3:N+3])/dy
        f[i,:] = f[i,:] - np.maximum(vhalf[:],0.0)*DyLess[:] - np.minimum(vhalf[:],0.0)*DyPlus[:]
    #
    return




##########################################################
# Schema centre + diffusion numerique (type Lax-Wendroff)
##########################################################
def lax_wendroff(f, phi, x, y, u, v, dx, dy, dt, N):
    # Variables locales
    convx = np.zeros(N)
    convy = np.zeros(N)
    diffx = np.zeros(N)
    diffy = np.zeros(N)
    uhalf = np.zeros(N)
    vhalf = np.zeros(N)
    #
    f[:,:] = 0.0
    # Advection suivant i
    for j in range(N):
        uhalf[:]   = 0.5*(u[1:N+1,j]+u[0:N,j])
        convx[:]   = uhalf[:]*(phi[4:N+4,j+3]-phi[2:N+2,j+3])/(2*dx)
        diffx[:]   = uhalf[:]**2*dt/2.0*(phi[4:N+4,j+3]+phi[2:N+2,j+3]-2.0*phi[3:N+3,j+3])/dx**2
        f[:,j] = f[:,j] - (convx[:] - diffx[:])
    #
    # Advection suivant j
    for i in range(N):
        vhalf[:]   = 0.5*(v[i,1:N+1]+v[i,0:N])
        convy[:]   = vhalf[:]*(phi[i+3,4:N+4]-phi[i+3,2:N+2])/(2*dy)
        diffy[:]   = vhalf[:]**2*dt/2.0*(phi[i+3,4:N+4]+phi[i+3,2:N+2]-2.0*phi[i+3,3:N+3])/dy**2
        f[i,:] = f[i,:] - (convy[:] - diffy[:])
    #
    return
