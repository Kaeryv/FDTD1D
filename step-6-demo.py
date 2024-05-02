from matplotlib.pyplot import hist2d, plot, show, draw, title, figure, pause, axis, grid, xlabel, axvline, axvspan, subplots
from numpy import zeros, ones, max, linspace, zeros_like, exp, complex64, absolute, pi, asarray
from math import pow
import numpy as np


# Initialization
# - grille
# - paramètres

def gaussian(x, x0, sigma):
    return exp(-pow( (x - x0) / sigma, 2.)   /2.)

TMAX = 2*5000 # Nombre d'itérations de durée dt
NF = 400
TP = 1e-9
eps = [] 
dbuffer=100
fmax = 1 # Normalisée
a = 0.5e-6
na = 50
ffield=0.365 # Fréquence track mode
dz = a / na
#dz =  10 * 1e-9 # [m]
#na = a / dz
c0 = 3e+8 # [m/s]
# Struture 1: pour diminuer réflexion dans un gamme spectrale (tête de fusée)
# Attention, étendre la gamme en fréquence vers f=1.6 sinon on ne voit pas la zone
if False:
    eps.extend([ 1 for i in range(100) ])
    eps.extend([ 3.46 for i in range(8) ])
    eps.extend([ 12 for i in range(60) ])
    eps.extend([ 3.46 for i in range(8) ])
    eps.extend([ 1 for i in range(100) ])

# Structure 2 Bragg avec défaut
if True:
    d1 = int(0.8*na)
    d2 = int(0.2*na)
    eps.extend([ 1 for i in range(dbuffer) ])
    for i in range(3):
        eps.extend([ 5 for i in range(d2) ])
        eps.extend([ 1 for i in range(d1) ])
    eps.extend([1 for i in range(d1) ])
    for i in range(3):
        eps.extend([ 5 for i in range(d2) ])
        eps.extend([ 1 for i in range(d1) ])
    eps.extend([ 1 for i in range(dbuffer) ])

# Structure 3: Miroir de Bragg du Joanopoulos
if False:
    d1 = int(0.8 * na)
    d2 = int(0.2 * na)
    eps.extend([1] * dbuffer)
    for i in range(6):
        eps.extend([1] *  d1) # 0.8 * a
        eps.extend([13] * d2) # 0.2 * a
    
    eps.extend([1] * dbuffer)

N = len(eps)
eps = asarray(eps)

# Le coefficient C0dt/dz = 0.5 (nombre de Courant)
m = 0.5


Ex = zeros((N))
Hy = zeros((N))
R = zeros((NF))
T = zeros((NF))


dt = dz / c0 / 2.0 # Formule de Courant
f = linspace(0, fmax*c0/a*dt, NF) # La fréquence maximale dépend de la résolution temporelle.

'''
    Set-up des plots
'''
x = np.arange(len(Ex)) * dz * 1e6
fig, (ax1, ax2) = subplots(2,1)
hEx = ax1.plot(x, Ex, 'r-')
hHy = ax1.plot(x, Hy, 'b-')
heps = ax1.plot(x, (eps-1)/max(eps), 'k-')
axis([0, N, -1, 1])

hTran, = ax2.plot(R, 'r-')
hRefl, = ax2.plot(T, 'b-')
ax2.set_xlim(0, f[-1]/dt/c0*a)
ax2.set_ylim(0, 1)

# Conditions PML
E1, E2, E3 = 0, 0, 0
H1, H2, H3 = 0, 0, 0

# Kernel pour la DiscrFourierTransform
K = exp(-1j * 2 * pi * f)
# Kernel pour la carte de champs fréq.
fc = ffield * c0 / a * dt 
Kc = exp(-1j * 2 * pi * fc)
Eomega = np.zeros_like(Ex, dtype='complex')
spectum_source = zeros_like(f, dtype=complex64)
spectum_reflec = zeros_like(f, dtype=complex64)
spectum_transm = zeros_like(f, dtype=complex64)

display=True
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

    Eomega += Kc**t * Ex
    if not (t % 20) and display:
        title(t)
        hEx[0].set_ydata(Ex)
        hHy[0].set_ydata(Hy)

        src = absolute(spectum_source)**2
        refl = absolute(spectum_reflec)**2 
        tran = absolute(spectum_transm)**2
        hRefl.set_xdata(f/dt/c0*a)
        hRefl.set_ydata(refl/src)
        hTran.set_xdata(f/dt/c0*a)
        hTran.set_ydata(tran/src)
        # En matlab: set(hEx,'Ydata', Ex)
        pause(TP)
    
    '''
    En matlab:
    h = plot(rand(10,1))
    for n = 1:20
        set(h,'Ydata',rand(10,1))
        pause(.0000001)
    end
    '''
show()
src = absolute(spectum_source)**2
refl = absolute(spectum_reflec)**2 
tran = absolute(spectum_transm)**2
fig, ax = subplots()
ax.plot(f/dt, refl/src, 'b-')
ax.plot(f/dt, tran/src, 'r-')
#axvspan(2.4, color='k', alpha=0.5)

title("Spectres de réflexion et transmission pour le randome")
ax.set_xlabel("Fréquence normalisée")
grid()
show()


fig, ax = subplots()
ax.plot(x, Eomega, 'b-')
ax1 = ax.twinx()
ax1.plot(x, eps, 'k-')
title(f"Champ Ex f={fc}")
ax.set_xlabel("Espace [microns]")
grid()
show()
