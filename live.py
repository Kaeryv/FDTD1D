from matplotlib.pyplot import hist2d, legend, plot, show, draw, title, figure, pause, axis, xlabel
from numpy import zeros, ones, zeros_like, complex64, exp
from math import pow, pi

from numpy.core.function_base import linspace



# Initialization
# - grille
# - paramètres

def gaussian(x, x0, sigma):
    return exp(  -pow( (x - x0) / sigma, 2.)   /2.)

N = 150 # Nombre de cellules de Yee
NF = 500 # nbr de valeurs de fréquences
TMAX = 1000 # Nombre d'itérations de durée dt

# on considère uniquement de l'air e=1

# Le coefficient c dt/dz = 0.5   
# Nombre de Courant  C0 = v dt / dz = 0.5 pour que ça marche en 1D
m = 0.5


Ex = zeros((N))
Hy = zeros((N))
eps = ones((N))
f = linspace(1/200, 1/60., NF)
K = exp(-2*pi*1j*f) # Partie de la TF indépendante de l'itération.
refl = zeros_like(f, dtype=complex64) # Juste un vecteur comme f, mais avec des zéros.
trsf = zeros_like(f, dtype=complex64) # Pas besoin de l'histoire de dtype en matlab
src  = zeros_like(f, dtype=complex64)


eps[50:100] = 3. # Un slab de 50 dz
#eps[180:200] = 3.
# Ex[149] = 1.0
# Ex[150] = 0.5
# Ex[151] = 1.0


E1, E2, E3 = 0, 0, 0
H1, H2, H3 = 0, 0, 0

source_z = 10

display = False

if display:
    hEx = plot(Ex, 'r-')
    hHy = plot(Hy, 'b-')
    heps = plot((eps-1)/max(eps-1), 'k-')
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

    ''' Calcul des spectres '''
    trsf += Ex[-1] * K**t
    refl += Ex[1] * K**t
    src += gaussian(t, 50, 10) * K**t

    ''' DISPLAY '''
    # t % 10 modulo
    if not (t % 2) and display:
        title(t)
        hEx[0].set_ydata(Ex)
        hHy[0].set_ydata(Hy)
        pause(0.00001)


dz = 10 * 1e-9 # m soit 10 nm => slab de 50 * 10 nm, soit 500nm de slab
# Nombre de Courant
# Courant = 0.5 = c dt / dz = m
c = 300000 * 1e3 # m/s
dt = dz * 0.5 / c # D'après la formule de Courant, en s
f = f / dt # en Hz
lamb = c / f # en m
if True:
    trsf = abs(trsf/2)**2
    refl = abs(refl/2)**2
    src = abs(src)**2
    plot(lamb/1e-9, trsf/src, 'r-', label='Trans')
    # hold on en matlab
    plot(lamb/1e-9, refl/src, 'b-', label='Refl')
    plot(lamb/1e-9, src/max(src), 'k-')
    xlabel("Fréquence en Hz")
    xlabel("Longueur d'onde en nm")
    legend()

show()