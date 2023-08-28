from numpy import pi
import matplotlib.pyplot as plt

dt = 0.0005  # time discretization. Keep this number low
N = 540  # inverse space discretization. Keep this number high!

# model constants
beta = 1.0  # thermodynamic beta: 1/kT
m1 = m2 = 1.0  # masses of Fo and F1

# model-specific parameters
gamma1 = gamma2 = 1000.0  # drag coefficients of Fo and F1

E0 = 2.0 # energy scale of Fo
Ecouple = 30.0 # energy scale of coupling between Fo and F1
E1 = 2.0 # energy scale of F1
mu_Hp = 4.0 #  mu_{H+}: energy INTO (positive) Fo by F1
mu_atp = 2.0 # mu_{ATP}: energy INTO (positive) F1 by Fo

n1 = 3.0  # number of minima in the potential of Fo
n2 = 30.0  # number of minima in the potential of F1
phase = 0.0  # how much sub-systems are offset from one another

dx = (2*pi)/N  # space discretization: total distance / number of points

time_check = dx/((0.5*(Ecouple+E0*n1)-mu_Hp)/(m1*gamma1) + (0.5*(Ecouple+E1*n2)-mu_atp)/(m2*gamma2))

if dt > time_check:
    # bail if user is stupid
    print("!!!TIME UNSTABLE!!! No use in going on. Aborting...\n")
    exit(1)
else:
    print("(All's good pardner) __^__^__")
    print("                  v  =('. ')=  |")
    print("                      |`[]``\  /")
    print("                      |______|/")

