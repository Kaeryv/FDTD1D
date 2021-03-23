from matplotlib.pyplot import hist2d, plot, show, draw, title, figure, pause, axis
from numpy import zeros, ones
from math import exp, pow



# Initialization
# - grille
# - paramètres

def gaussian(x, x0, sigma):
    return exp(  -pow( (x - x0) / sigma, 2.)   /2.)

N = 300 # Nombre de cellules de Yee

# on considère uniquement de l'air e=1

# Le coefficient C0 dt/dz = 0.5   
m = 0.5

TMAX = 1000 # Nombre d'itérations de durée dt

Ex = zeros((N))
Hy = zeros((N))
eps = ones((N))

eps[100:150] = 3.
eps[180:200] = 3.
# Ex[149] = 1.0
# Ex[150] = 0.5
# Ex[151] = 1.0
# for i in range(100, 200):
#     Ex[i] = gaussian(i, 150, 10)

#fig = figure()
hEx = plot(Ex, 'r-')
hHy = plot(Hy, 'b-')
heps = plot((eps-1)/max(eps-1), 'b-')

E1, E2, E3 = 0, 0, 0
H1, H2, H3 = 0, 0, 0

source_z = 10

axis([0, N, -2.2, 2.2])
for t in range(TMAX):

    # Ne pas oublier les conditions initiales! ou frontières
    # Calculer Hy en fonction de Ex
    for k in range(N-1):
        Hy[k] = Hy[k] - m * (Ex[k+1]-Ex[k])

    # pour k = N, condition PEC
    Hy[N-1] = Hy[N-1] - m * (E3 - Ex[N-1])

    # On propage ici vers la gauche la valeur en bout de grille à gauche. L'ordre des opérations est important.
    H3 = H2
    H2 = H1
    H1 = Hy[0] # Hy(1) en matlab, premier élément de la grille
    # La gaussienne est décalée de 1.0 unités de temps plus 0.5 à cause d'un décalage dans l'espace et le temps.
    Hy[source_z-1] += gaussian(t-0.5-1.0, 50, 10)
    # L'explication du source_z-1 provient de la décomposition de la densité de courant.
    # Calculer Ex en fonction de Hy

    # Pour k = 0, condition PEC
    Ex[0] = Ex[0] - m/eps[k] * (Hy[0] - H3)
    for k in range(1, N):
        Ex[k] = Ex[k] - m/eps[k] * (Hy[k]-Hy[k-1])

    E3 = E2
    E2 = E1
    E1 = Ex[N-1] # Ex[N] en matlab, dernier élément de la grille Ex

    Ex[source_z] += gaussian(t, 50, 10)
    # t % 10 modulo
    if not (t % 2):
        title(t)
        hEx[0].set_ydata(Ex)
        hHy[0].set_ydata(Hy)
        # En matlab: set(hEx,'Ydata', Ex)
        pause(0.00001)
    
    '''
    En matlab:
    h = plot(rand(10,1))
    for n = 1:20
        set(h,'Ydata',rand(10,1))
        pause(.0000001)
    end
    '''

show()