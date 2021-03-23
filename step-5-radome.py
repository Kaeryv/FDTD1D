from matplotlib.pyplot import hist2d, plot, show, draw, title, figure, pause, axis, grid, xlabel, axvline, axvspan, legend
from numpy import zeros, ones, max, linspace, zeros_like, exp, complex64, absolute, pi, asarray
from math import pow


# Initialization
# - grille
# - paramètres

def gaussian(x, x0, sigma):
    return exp(  -pow( (x - x0) / sigma, 2.)   /2.)

eps = list()
eps.extend([ 1 for i in range(100) ])
#eps.extend([ 3.46 for i in range(8) ])
eps.extend([ 12 for i in range(75) ])
#eps.extend([ 3.46 for i in range(8) ])
eps.extend([ 1 for i in range(100) ])

N = len(eps)
eps = asarray(eps)

# Le coefficient C0dt/dz = 0.5
m = 0.5

TMAX = 5000 # Nombre d'itérations de durée dt

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

    # pour k = N, condition PEC
    Hy[N-1] = Hy[N-1] - m * (E3 - Ex[N-1])
    H3 = H2
    H2 = H1
    H1 = Hy[0]
    Hy[9] += gaussian(t-1.5, 50, 10)
    # Calculer Ex en fonction de Hy

    # Pour k = 0, condition PEC
    Ex[0] = Ex[0] - m/eps[k] * (Hy[0] - H3)
    for k in range(1, N):
        Ex[k] = Ex[k] - m/eps[k] * (Hy[k]-Hy[k-1])
    E3 = E2
    E2 = E1
    E1 = Ex[-1]

    Ex[10] += gaussian(t, 50, 10)

    spectum_source += K**t * gaussian(t, 50, 10)
    if t > 100:
        spectum_reflec += K**t * Ex[1]
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

dz = 2e-3 # m
c0 = 3e+8 # m/s
dt = dz/2/c0 #   dt = dz / 2c0
plot(f/dt, src/max(src), 'k-', label='profil source')
plot(f/dt, tran/src/4, 'r-', label='transmission')
plot(f/dt, refl/src/4, 'b-', label='reflexion')
legend()
#axvspan(2.39, 2.41, color='k', alpha=0.5)

title("Spectres de réflexion et transmission pour le radome")
xlabel("Fréquence [GHz]")
grid()
show()
