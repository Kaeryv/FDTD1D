from matplotlib.pyplot import hist2d, plot, show, draw, title, figure, pause, axis, grid, xlabel, axvline, axvspan, legend
from numpy import zeros, ones, max, linspace, zeros_like, exp, complex64, absolute, pi, asarray
from math import pow


# Initialization
# - grille
# - paramètres

def gaussian(x, x0, sigma):
    return exp(  -pow( (x - x0) / sigma, 2.)   /2.)

N = 300 # Nombre de cellules de Yee
eps = list() #ones((N))
n1= 1.5
n2 = 2.
d1 = 16
d2 = 12

eps.extend([ 1 for i in range(50) ])
for i in range(10):
    eps.extend([ n2**2 for i in range(d2) ])
    eps.extend([ n1**2 for i in range(d1) ])

eps.extend([ 1 for i in range(50) ])

N = len(eps)
eps = asarray(eps)
# on considère uniquement de l'air e=1

# Le coefficient C0dt/dz = 0.5
m = 0.5

TMAX = 10000 # Nombre d'itérations de durée dt

Ex = zeros((N))
Hy = zeros((N))


E1, E2, E3 = 0, 0, 0
H1, H2, H3 = 0, 0, 0
f = linspace(0,0.017,2000)
K = exp(-1j * 2 * pi * f)

spectum_source = zeros_like(f, dtype=complex64)
spectum_reflec = zeros_like(f, dtype=complex64)
spectum_transm = zeros_like(f, dtype=complex64)
display=False


if display:
    hEx = plot(Ex, 'r-')
    hHy = plot(Hy, 'b-')
    heps = plot((eps-1)/max(eps), 'k-')
    axis([0, N, -1, 1])



for t in range(TMAX):

    # Ne pas oublier les conditions initiales! ou frontières
    # Calculer Hy en fonction de Ex
    for k in range(N-1):
        Hy[k] = Hy[k] - m * (Ex[k+1]-Ex[k])

    # pour k = N, condition PML
    Hy[N-1] = Hy[N-1] - m * (E3 - Ex[N-1])
    H3 = H2
    H2 = H1
    H1 = Hy[0]
    Hy[9] += gaussian(t-1.5, 50, 10)
    # Calculer Ex en fonction de Hy

    # Pour k = 0, condition PML
    Ex[0] = Ex[0] - m/eps[k] * (Hy[0] - H3)
    for k in range(1, N):
        Ex[k] = Ex[k] - m/eps[k] * (Hy[k]-Hy[k-1])
    E3 = E2
    E2 = E1
    E1 = Ex[-1]

    Ex[10] += gaussian(t, 50, 10)
    spectum_reflec += K**t * Ex[3]
    spectum_source += K**t * gaussian(t, 50, 10)
    spectum_transm += K**t * Ex[-2]
    if not (t % 20) and display:
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
figure()
src = absolute(spectum_source)**2
refl = absolute(spectum_reflec)**2
tran = absolute(spectum_transm)**2

dz = 10 *1e-9 # m
c0 = 3e+8 # m/s
dt = dz/2/c0 #   dt = dz / 2c0
f = f/dt

l = c0 / f * 1e+9
plot(l, src/max(src), 'k-', label='profil source')
plot(l, tran/src/4, 'r-', label='transmission')
plot(l, refl/src/4, 'b-', label='reflexion')
axis([700, 1200, 0, 1])
legend()
#axvspan(2.39, 2.41, color='k', alpha=0.5)

title("Spectres de réflexion et transmission pour le bouclier, 10 bicouches.")
xlabel("Longueur d'onde [nm]")
grid()
show()
