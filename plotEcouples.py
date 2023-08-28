from numpy import linspace, loadtxt, append, pi, empty, sqrt, zeros, asarray, trapz
from numpy import array as nparray
import math
import os
import matplotlib.pyplot as plt
import matplotlib.font_manager
from matplotlib import rc
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
rc('text', usetex=True)

path = os.getcwd()

if not os.path.exists(path+'/plots'):  # create a diretory to output plots
    os.mkdir(path+'/plots')

N = 360  # N x N grid is used for Fokker-Planck simulations
dx = 2 * math.pi / N  # spacing between gridpoints
positions = linspace(0, 2 * math.pi - dx, N)  # gridpoints
timescale = 1.5 * 10**4  # conversion factor between simulation and experimental timescale

E0 = 2.0  # barrier height Fo
E1 = 2.0  # barrier height F1
Ecouple_arr = [1.0, 10.0, 20.0, 30.0]  # coupling strengths
psi_1 = 4.0  # chemical driving force on Fo
psi_2 = -2.0  # chemical driving force on F1
num_minima1 = 3.0  # number of barriers in Fo's landscape
num_minima2 = linspace(3.0,30.0,28)  # number of barriers in F1's landscape

phase = 0.0

def plot_power_efficiency_Ecouple(path,timescale, num_minima1, num_minima2, E0, E1, Ecouple_arr, psi_1, psi_2):  # plot power and efficiency vs number of barriers n2

    output_file_name = (
            path + "P_ATP_eff_" + "E0_{0}_E1_{1}_psi1_{2}_psi2_{3}_n1_{4}" + "_.pdf")
    f, axarr = plt.subplots(2, 1, sharex='all', sharey='none', figsize=(6, 8))

    # power plot
    axarr[0].axhline(0, color='black', linewidth=1)  # x-axis
    optimal_n2 = num_minima1/(psi_2/psi_1 + sqrt(1+(psi_2/psi_1)*(psi_2/psi_1)))
    axarr[0].axvline(optimal_n2, color='black', linestyle='--', linewidth=1)  # lining up features in the two plots
    axarr[0].fill_between([1, 250], 0, 31, facecolor='grey', alpha=0.2)  # shading power output
    
    for Ecouple in Ecouple_arr:
    # Fokker-Planck results (barriers)
        i = 0  # only use phase=0 data
        power_y_array = []
        for ii, minima2 in enumerate(num_minima2):
            input_file_name = (path + "/Ecouple" + str(int(Ecouple)) + "/plots/" + "flux_power_efficiency_"
                            + "E0_{0}_E1_{1}_psi1_{2}_psi2_{3}_n1_{4}_n2_{5}_Ecouple_{6}" + "_outfile.dat")
            try:
                data_array = loadtxt(
                    input_file_name.format(E0, E1, psi_1, psi_2, num_minima1, minima2, Ecouple),
                    usecols=(0, 4))
                
                power_y = nparray(data_array[1])
                power_y_array = append(power_y_array, power_y)
            except OSError:
                print('Missing file flux')

        axarr[0].plot(num_minima2, -2.0*pi*power_y_array*timescale, 'o', color='C1', label='Fokker-Planck', markersize=8)
        
    
    tight_power_array = empty(0) # curve at tight coupling limit
    for ii, minima2 in enumerate(num_minima2):
        tight_power = -(psi_2 + (num_minima1*psi_1/minima2-psi_2)/(1+(num_minima1/minima2)*(num_minima1/minima2)))*psi_2
        tight_power_array = append(tight_power_array, tight_power)
    axarr[0].plot(num_minima2, tight_power_array, 'o', color='black', label='Tight coupling limit', markersize=8)

    axarr[0].yaxis.offsetText.set_fontsize(14)
    axarr[0].tick_params(axis='y', labelsize=14)
    axarr[0].set_ylabel(r'$\beta \mathcal{P}_{\rm ATP} (\rm s^{-1}) $', fontsize=20)
    axarr[0].spines['right'].set_visible(False)
    axarr[0].spines['top'].set_visible(False)
    axarr[0].spines['bottom'].set_visible(False)
    axarr[0].set_xlim((1.7, 135))
    axarr[0].set_ylim((-60, 31))
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
    
    # Fokker-Planck results (barriers)
    for Ecouple in Ecouple_arr:
        eff_array = []
        for ii, minima2 in enumerate(num_minima2):
            input_file_name = (
                    path + "/plots/flux_power_efficiency_"
                    + "E0_{0}_E1_{1}_psi1_{2}_psi2_{3}_n1_{4}_n2_{5}_Ecouple_{6}" + "_outfile.dat")
            try:
                data_array = loadtxt(
                    input_file_name.format(E0, E1, psi_1, psi_2, num_minima1, minima2, Ecouple), usecols=5)

                eff_array = append(eff_array, data_array)
            except OSError:
                print('Missing file efficiency')
        axarr[1].plot(num_minima2, eff_array, 'o', color='C1', markersize=8)

    tight_eff_array = empty(0)
    for ii, minima2 in enumerate(num_minima2):
        input_file_name = (
                path + "/plots/flux_power_efficiency_"
                + "E0_{0}_E1_{1}_psi1_{2}_psi2_{3}_n1_{4}_n2_{5}_Ecouple_{6}" + "_outfile.dat")
        
        tight_eff = -num_minima1*psi_2/(minima2*psi_1)
        tight_eff_array = append(tight_eff_array, tight_eff)

    axarr[1].plot(num_minima2, tight_eff_array, 'o', color='black', markersize=8)

    axarr[1].set_xlabel(r'$\beta E_{\rm couple}$', fontsize=20)
    axarr[1].set_ylabel(r'$\eta / \eta^{\rm max}$', fontsize=20)
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
    f.savefig(output_file_name.format(E0, E1, psi_1, psi_2, num_minima1))

if __name__ == "__main__":

    path = os.getcwd()
    plot_power_efficiency_Ecouple(path,timescale, num_minima1, num_minima2, E0, E1, Ecouple_arr, psi_1, psi_2)
