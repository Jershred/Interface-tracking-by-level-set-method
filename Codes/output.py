import numpy as np
import matplotlib.pyplot as plt

########################################################
# Affiche la superposition des 2 solutions phi0 et phi1
########################################################
def write_output(x, y, phi0, phi1, u, v, N, timeSchemeDef, spaceSchemeDef, it_out, testCase):
    xout     = np.zeros(N)
    yout     = np.zeros(N)
    uout     = np.zeros((N, N))
    vout     = np.zeros((N, N))
    #
    xout[:]  = 0.5*(x[0:N] + x[1:N+1])
    yout[:]  = 0.5*(y[0:N] + y[1:N+1])
    #
    for j in range(N):
        uout[:,j] = 0.5*(u[0:N,j]+u[1:N+1,j])
    #
    for i in range(N):
        vout[i,:] = 0.5*(v[i,0:N]+v[i,1:N+1])
    #    #
    error = compute_errors(phi0, phi1, N)
    #
    fig2, ax2 = plt.subplots()
    #
    CS1 = ax2.contour(xout, yout, np.transpose(phi0[3:N+3,3:N+3]), levels=[0], colors='blue')
    CS2 = ax2.contour(xout, yout, np.transpose(phi1[3:N+3,3:N+3]), levels=[0], colors='red')
    #
    # Obtenir les coordonnées du polygone
    contour1 = CS1.collections[0]
    contour2 = CS2.collections[0]
    
    vs1 = contour1.get_paths()[0].vertices
    vs2 = contour2.get_paths()[0].vertices
    
    # Calculer l'aire du polygone grâce au sommet
    a1 = area(vs1)
    a2 = area(vs2)
    
    print('Initial area :', a1)
    print('Final area :', a2)
    
    print('Area difference :', abs(a2-a1))
    print('Area realtive error :', abs(a2-a1)/a1*100) #erreur relative en pourcent
    #
    ax2.grid('true')
    #
    ax2.set(xlabel='x')
    ax2.set(ylabel='y') 
    ax2.text(0.02, 0.02, 'error= ' + str(format(error,'.2e')), fontsize=10, transform = ax2.transAxes) 
    ax2.set_title('{}'.format(timeSchemeDef +' / '+ spaceSchemeDef)) 
    h1,_ = CS1.legend_elements()
    h2,_ = CS2.legend_elements()
    ax2.legend([h1[0], h2[0]], ['it=0', 'it={}'.format(it_out)])
    #
    plt.savefig('Images/{}.png'.format(testCase+'_'+timeSchemeDef +'_'+spaceSchemeDef)) #Sauvegarde les graphs dans un dossier
    plt.show()
    
    return

##########################################################################
# Use Green's theorem to compute the area enclosed by the given contour
##########################################################################

def area(vs):
    a = 0
    x0,y0 = vs[0]
    for [x1,y1] in vs[1:]:
        dx = x1-x0
        dy = y1-y0
        a += 0.5*(y0*dx - x0*dy)
        x0 = x1
        y0 = y1
    return a


#####################################
# Ecrit une solution au format ASCII
#####################################
def write_output_ASCII(x, y, phi, t, N):
    filename = 'output_'+str(t)+'.dat'
    xc = np.zeros(N)
    yc = np.zeros(N)
    xc[:] = 0.5*(x[0:N] + x[1:N+1])
    yc[:] = 0.5*(y[0:N] + y[1:N+1])
    with open(filename, 'w') as f:
        for j in range(N):
            for i in range(N):
                f.write('%.10f %.10f %.10f\n' % (xc[i], yc[j], phi[i,j]))
    return


########################################
# Calcul l'erreur L2 entre phi0 et phi1
########################################
def compute_errors(phi0, phi1, N):
    return ((1.0/(N*N))*np.sum((phi1[3:N+3,3:N+3]-phi0[3:N+3,3:N+3])**2))**0.5
