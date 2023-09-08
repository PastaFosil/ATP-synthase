from numpy import linspace, loadtxt, append, pi, empty, sqrt, zeros, asarray, trapz
from numpy import array as nparray
import math
import os
import matplotlib.pyplot as plt
import matplotlib.font_manager
from matplotlib import rc

from itertools import cycle

rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
rc('text', usetex=True)
colors = cycle(['#c1272d','#0000a7','#eecc16','#008176'])

path = os.getcwd()

N = 360  # N x N grid is used for Fokker-Planck simulations
dx = 2 * math.pi / N  # spacing between gridpoints
positions = linspace(0, 2 * math.pi - dx, N)  # gridpoints
timescale = 1.5 * 10**4  # conversion factor between simulation and experimental timescale

E0 = 2.0  # barrier height Fo
E1 = 2.0  # barrier height F1
Ecouple_arr = [1.0, 10.0, 20.0, 30.0]  # coupling strengths
mu0 = 4.0  # chemical driving force on Fo
mu1 = -2.0  # chemical driving force on F1
num_minima0 = 3.0  # number of barriers in Fo's landscape
num_minima1 = linspace(3.0,30.0,28)  # number of barriers in F1's landscape

phase = 0.0

def plot_power_efficiency_Ecouple(path,timescale, num_minima0, num_minima1, E0, E1, Ecouple_arr, mu0, psi_2):  # plot power and efficiency vs number of barriers n2

    f, axarr = plt.subplots(2, 1, sharex='all', sharey='none', figsize=(6, 8))

    # power plot
    axarr[0].axhline(0, color='black', linewidth=1)  # x-axis
    optimal_n2 = num_minima0/(psi_2/mu0 + sqrt(1+(psi_2/mu0)*(psi_2/mu0)))
    axarr[0].axvline(optimal_n2, color='black', linestyle='--', linewidth=1)  # lining up features in the two plots
    axarr[0].fill_between([1, 250], 0, 31, facecolor='grey', alpha=0.2)  # shading power output
    
    tight_power_array = empty(0) # curve at tight coupling limit
    for ii, minima1 in enumerate(num_minima1):
        tight_power = -(psi_2 + (num_minima0*mu0/minima1-psi_2)/(1+(num_minima0/minima1)*(num_minima0/minima1)))*psi_2
        tight_power_array = append(tight_power_array, tight_power)
    axarr[0].plot(num_minima1, tight_power_array, 'o', color='black', label='Tight coupling limit', markersize=8)

    axarr[0].yaxis.offsetText.set_fontsize(14)
    axarr[0].tick_params(axis='y', labelsize=14)
    axarr[0].set_ylabel(r'$\beta \mathcal{P}_{\rm ATP} (\rm s^{-1}) $', fontsize=20)
    axarr[0].spines['right'].set_visible(False)
    axarr[0].spines['top'].set_visible(False)
    axarr[0].spines['bottom'].set_visible(False)
    axarr[0].set_xlim((1.7, 135))
    axarr[0].set_ylim((-5, 5))
    axarr[0].set_yticks([-50, -25, 0, 25])

    leg = axarr[0].legend(title=r'$\beta E_{\rm o} = \beta E_1$', fontsize=14, loc='lower right', frameon=False)
    leg_title = leg.get_title()
    leg_title.set_fontsize(14)

    #####################################################
    # efficiency plot
    axarr[1].axhline(0, color='black', linewidth=1)  # x axis
    axarr[1].axvline(optimal_n2, color='black', linestyle='--', linewidth=1)  # lining up features in the two plots
    axarr[1].axhline(1, color='black', linestyle=':', linewidth=1)  # max efficiency
    axarr[1].fill_between([1, 250], 0, 1, facecolor='grey', alpha=0.2)  # shading power output
    
    tight_eff_array = empty(0)
    for ii, minima1 in enumerate(num_minima1):
        tight_eff = -num_minima0*psi_2/(minima1*mu0)
        tight_eff_array = append(tight_eff_array, tight_eff)

    axarr[1].plot(num_minima1, tight_eff_array, 'o', color='black', markersize=8)

    axarr[1].set_xlabel(r'$n_0$', fontsize=20)
    axarr[1].set_ylabel(r'$\eta$', fontsize=20)
    axarr[1].set_xscale('log')
    axarr[1].set_xlim((1.7, 135))
    axarr[1].set_ylim((-0.5, 1.05))
    axarr[1].spines['right'].set_visible(False)
    axarr[1].spines['top'].set_visible(False)
    axarr[1].spines['bottom'].set_visible(False)
    axarr[1].set_yticks([-0.5, 0, 0.5, 1.0])
    axarr[1].tick_params(axis='both', labelsize=14)

    f.text(0.05, 0.95, r'$\mathbf{a)}$', ha='center', fontsize=20)
    f.text(0.05, 0.48, r'$\mathbf{b)}$', ha='center', fontsize=20)
    f.subplots_adjust(hspace=0.01)

    f.tight_layout()
    plt.show()
    #f.savefig("tight_coupling.pdf")

if __name__ == "__main__":

    path = os.getcwd()
    plot_power_efficiency_Ecouple(path,timescale, num_minima0, num_minima1, E0, E1, Ecouple_arr, mu0, mu1)
