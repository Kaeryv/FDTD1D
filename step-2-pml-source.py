from matplotlib.pyplot import hist2d, plot, show, draw, title, figure, pause, axis
from numpy import zeros
from math import exp, pow

# Initialization
# - grille
# - paramètres

def gaussian(x, x0, sigma):
    return exp(  -pow( (x - x0) / sigma, 2.)   /2.)

N = 300 # Nombre de cellules de Yee

# on considère uniquement de l'air e=1

# Le coefficient C0dt/dz = 0.5
m = 0.5

TMAX = 200 # Nombre d'itérations de durée dt

Ex = zeros((N))
Hy = zeros((N))


#for i in range(100, 200):
#    Ex[i] = gaussian(i, 150, 10)

#fig = figure()
hEx = plot(Ex, 'r-')
hHy = plot(Hy, 'b-')
axis([0, TMAX, -1, 1])
E1, E2, E3 = 0, 0, 0
H1, H2, H3 = 0, 0, 0
for t in range(TMAX):

    # Ne pas oublier les conditions initiales! ou frontières
    # Calculer Hy en fonction de Ex
    for k in range(N-1):
        Hy[k] = Hy[k] - m * (Ex[k+1]-Ex[k])

    # pour k = N, condition PEC
    Hy[N-1] = Hy[N-1] - m * (E3 - Ex[N-1])
    Hy[9] += gaussian(t-1.5, 50, 10)
    H3 = H2
    H2 = H1
    H1 = Hy[0]
    # Calculer Ex en fonction de Hy

    # Pour k = 0, condition PEC
    Ex[0] = Ex[0] - m * (Hy[0] - H3)
    for k in range(1, N):
        Ex[k] = Ex[k] - m * (Hy[k]-Hy[k-1])
    E3 = E2
    E2 = E1
    E1 = Ex[-1]
    Ex[10] += gaussian(t, 50, 10)

    if not (t % 5):
        title(t)
        hEx[0].set_ydata(Ex)
        hHy[0].set_ydata(Hy)
        # En matlab: set(hEx,'Ydata', Ex)
        pause(0.00000001)
    
    '''
    En matlab:
    h = plot(rand(10,1))
    for n = 1:20
        set(h,'Ydata',rand(10,1))
        pause(.0000001)
    end
    '''

show()