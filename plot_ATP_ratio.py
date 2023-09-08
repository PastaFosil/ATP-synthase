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
Ecouple = 10.0  # coupling strengths
mu0 = 4.0  # chemical driving force on Fo
mu1 = -2.0  # chemical driving force on F1
num_minima0 = 3.0  # number of barriers in Fo's landscape
num_minima1 = linspace(3.0,30.0,28)  # number of barriers in F1's landscape

phase = 0.0

def calc_flux(p_now, drift_at_pos, diffusion_at_pos, flux_array, N):
    # explicit update of the corners
    # first component
    flux_array[0, 0, 0] = (
        (drift_at_pos[0, 0, 0]*p_now[0, 0])
        -(diffusion_at_pos[0, 1, 0]*p_now[1, 0]-diffusion_at_pos[0, N-1, 0]*p_now[N-1, 0])/(2.0*dx)
        -(diffusion_at_pos[1, 0, 1]*p_now[0, 1]-diffusion_at_pos[1, 0, N-1]*p_now[0, N-1])/(2.0*dx)
        )
    flux_array[0, 0, N-1] = (
        (drift_at_pos[0, 0, N-1]*p_now[0, N-1])
        -(diffusion_at_pos[0, 1, N-1]*p_now[1, N-1]-diffusion_at_pos[0, N-1, N-1]*p_now[N-1, N-1])/(2.0*dx)
        -(diffusion_at_pos[1, 0, 0]*p_now[0, 0]-diffusion_at_pos[1, 0, N-2]*p_now[0, N-2])/(2.0*dx)
        )
    flux_array[0, N-1, 0] = (
        (drift_at_pos[0, N-1, 0]*p_now[N-1, 0])
        -(diffusion_at_pos[0, 0, 0]*p_now[0, 0]-diffusion_at_pos[0, N-2, 0]*p_now[N-2, 0])/(2.0*dx)
        -(diffusion_at_pos[1, N-1, 1]*p_now[N-1, 1]-diffusion_at_pos[1, N-1, N-1]*p_now[N-1, N-1])/(2.0*dx)
        )
    flux_array[0, N-1, N-1] = (
        (drift_at_pos[0, N-1, N-1]*p_now[N-1, N-1])
        -(diffusion_at_pos[0, 0, N-1]*p_now[0, N-1]-diffusion_at_pos[0, N-2, N-1]*p_now[N-2, N-1])/(2.0*dx)
        -(diffusion_at_pos[1, N-1, 0]*p_now[N-1, 0]-diffusion_at_pos[1, N-1, N-2]*p_now[N-1, N-2])/(2.0*dx)
        )

    # second component
    flux_array[1, 0, 0] = (
        (drift_at_pos[1, 0, 0]*p_now[0, 0])
        -(diffusion_at_pos[2, 1, 0]*p_now[1, 0]-diffusion_at_pos[2, N-1, 0]*p_now[N-1, 0])/(2.0*dx)
        -(diffusion_at_pos[3, 0, 1]*p_now[0, 1]-diffusion_at_pos[3, 0, N-1]*p_now[0, N-1])/(2.0*dx)
        )
    flux_array[1, 0, N-1] = (
        (drift_at_pos[1, 0, N-1]*p_now[0, N-1])
        -(diffusion_at_pos[2, 1, N-1]*p_now[1, N-1]-diffusion_at_pos[2, N-1, N-1]*p_now[N-1, N-1])/(2.0*dx)
        -(diffusion_at_pos[3, 0, 0]*p_now[0, 0]-diffusion_at_pos[3, 0, N-2]*p_now[0, N-2])/(2.0*dx)
        )
    flux_array[1, N-1, 0] = (
        (drift_at_pos[1, N-1, 0]*p_now[N-1, 0])
        -(diffusion_at_pos[2, 0, 0]*p_now[0, 0]-diffusion_at_pos[2, N-2, 0]*p_now[N-2, 0])/(2.0*dx)
        -(diffusion_at_pos[3, N-1, 1]*p_now[N-1, 1]-diffusion_at_pos[3, N-1, N-1]*p_now[N-1, N-1])/(2.0*dx)
        )
    flux_array[1, N-1, N-1] = (
        (drift_at_pos[1, N-1, N-1]*p_now[N-1, N-1])
        -(diffusion_at_pos[2, 0, N-1]*p_now[0, N-1]-diffusion_at_pos[2, N-2, N-1]*p_now[N-2, N-1])/(2.0*dx)
        -(diffusion_at_pos[3, N-1, 0]*p_now[N-1, 0]-diffusion_at_pos[3, N-1, N-2]*p_now[N-1, N-2])/(2.0*dx)
        )

    for i in range(1, N-1):
        # explicitly update for edges not corners
        # first component
        flux_array[0, 0, i] = (
            (drift_at_pos[0, 0, i]*p_now[0, i])
            -(diffusion_at_pos[0, 1, i]*p_now[1, i]-diffusion_at_pos[0, N-1, i]*p_now[N-1, i])/(2.0*dx)
            -(diffusion_at_pos[1, 0, i+1]*p_now[0, i+1]-diffusion_at_pos[1, 0, i-1]*p_now[0, i-1])/(2.0*dx)
            )
        flux_array[0, i, 0] = (
            (drift_at_pos[0, i, 0]*p_now[i, 0])
            -(diffusion_at_pos[0, i+1, 0]*p_now[i+1, 0]-diffusion_at_pos[0, i-1, 0]*p_now[i-1, 0])/(2.0*dx)
            -(diffusion_at_pos[1, i, 1]*p_now[i, 1]-diffusion_at_pos[1, i, N-1]*p_now[i, N-1])/(2.0*dx)
            )
        flux_array[0, N-1, i] = (
            (drift_at_pos[0, N-1, i]*p_now[N-1, i])
            -(diffusion_at_pos[0, 0, i]*p_now[0, i]-diffusion_at_pos[0, N-2, i]*p_now[N-2, i])/(2.0*dx)
            -(diffusion_at_pos[1, N-1, i+1]*p_now[N-1, i+1]-diffusion_at_pos[1, N-1, i-1]*p_now[N-1, i-1])/(2.0*dx)
            )
        flux_array[0, i, N-1] = (
            (drift_at_pos[0, i, N-1]*p_now[i, N-1])
            -(diffusion_at_pos[0, i+1, N-1]*p_now[i+1, N-1]-diffusion_at_pos[0, i-1, N-1]*p_now[i-1, N-1])/(2.0*dx)
            -(diffusion_at_pos[1, i, 0]*p_now[i, 0]-diffusion_at_pos[1, i, N-2]*p_now[i, N-2])/(2.0*dx)
            )

        # second component
        flux_array[1, 0, i] = (
            (drift_at_pos[1, 0, i]*p_now[0, i])
            -(diffusion_at_pos[2, 1, i]*p_now[1, i]-diffusion_at_pos[2, N-1, i]*p_now[N-1, i])/(2.0*dx)
            -(diffusion_at_pos[3, 0, i+1]*p_now[0, i+1]-diffusion_at_pos[3, 0, i-1]*p_now[0, i-1])/(2.0*dx)
            )
        flux_array[1, i, 0] = (
            (drift_at_pos[1, i, 0]*p_now[i, 0])
            -(diffusion_at_pos[2, i+1, 0]*p_now[i+1, 0]-diffusion_at_pos[2, i-1, 0]*p_now[i-1, 0])/(2.0*dx)
            -(diffusion_at_pos[3, i, 1]*p_now[i, 1]-diffusion_at_pos[3, i, N-1]*p_now[i, N-1])/(2.0*dx)
            )
        flux_array[1, N-1, i] = (
            (drift_at_pos[1, N-1, i]*p_now[N-1, i])
            -(diffusion_at_pos[2, 0, i]*p_now[0, i]-diffusion_at_pos[2, N-2, i]*p_now[N-2, i])/(2.0*dx)
            -(diffusion_at_pos[3, N-1, i+1]*p_now[N-1, i+1]-diffusion_at_pos[3, N-1, i-1]*p_now[N-1, i-1])/(2.0*dx)
            )
        flux_array[1, i, N-1] = (
            (drift_at_pos[1, i, N-1]*p_now[i, N-1])
            -(diffusion_at_pos[2, i+1, N-1]*p_now[i+1, N-1]-diffusion_at_pos[2, i-1, N-1]*p_now[i-1, N-1])/(2.0*dx)
            -(diffusion_at_pos[3, i, 0]*p_now[i, 0]-diffusion_at_pos[3, i, N-2]*p_now[i, N-2])/(2.0*dx)
            )

        # for points with well defined neighbours
        for j in range(1, N-1):
            # first component
            flux_array[0, i, j] = (
                (drift_at_pos[0, i, j]*p_now[i, j])
                -(diffusion_at_pos[0, i+1, j]*p_now[i+1, j]-diffusion_at_pos[0, i-1, j]*p_now[i-1, j])/(2.0*dx)
                -(diffusion_at_pos[1, i, j+1]*p_now[i, j+1]-diffusion_at_pos[1, i, j-1]*p_now[i, j-1])/(2.0*dx)
                )
            # second component
            flux_array[1, i, j] = (
                (drift_at_pos[1, i, j]*p_now[i, j])
                -(diffusion_at_pos[2, i+1, j]*p_now[i+1, j]-diffusion_at_pos[2, i-1, j]*p_now[i-1, j])/(2.0*dx)
                -(diffusion_at_pos[3, i, j+1]*p_now[i, j+1]-diffusion_at_pos[3, i, j-1]*p_now[i, j-1])/(2.0*dx)
                )

def flux_power_efficiency(path, 
                          dx,N,num_minima0,num_minima1,phase,E0,E1,Ecouple,mu0,mu1): # processing of raw data
    
    phase_array = nparray([phase])
    psi1_array = nparray([mu0])
    psi2_array = nparray([mu1])

    for mu0 in psi1_array:
        for mu1 in psi2_array:
            integrate_flux_X = empty(phase_array.size)
            integrate_flux_Y = empty(phase_array.size)
            integrate_power_X = empty(phase_array.size)
            integrate_power_Y = empty(phase_array.size)
            efficiency_ratio = empty(phase_array.size)

            for minima1 in num_minima1:
                for ii, phase_shift in enumerate(phase_array):
                    input_file_name = (path +"/results/"+
                                       "reference_E0_{0}_Ecouple_{1}_E1_{2}_psi1_{3}_psi2_{4}_n1_{5}_n2_{6}_phase_{7}" +
                                       "_outfile.dat")

                    output_file_name = (path + "/plots/" + "flux_power_efficiency_" +
                                        "E0_{0}_E1_{1}_psi1_{2}_psi2_{3}_n1_{4}_n2_{5}_Ecouple_{6}" + "_outfile.dat")

                    print("Calculating flux for " + f"mu0 = {mu0}, mu1 = {mu1}, " +
                          f"Ecouple = {Ecouple}, num_minima0 = {num_minima0}, num_minima1 = {minima1}")
                    
                    try:
                        data_array = loadtxt(input_file_name.format(E0, Ecouple, E1, mu0, mu1, num_minima0,
                                                                    minima1, phase_shift),
                                             usecols=(0, 3, 4, 5, 6, 7, 8))
                        N = int(sqrt(len(data_array)))  # check grid size
                        print('Grid size: ', N)

                        prob_ss_array = data_array[:, 0].reshape((N, N))
                        drift_at_pos = data_array[:, 1:3].T.reshape((2, N, N))
                        diffusion_at_pos = data_array[:, 3:].T.reshape((4, N, N))

                        flux_array = zeros((2, N, N))
                        calc_flux(prob_ss_array, drift_at_pos, diffusion_at_pos, flux_array, N)
                        flux_array = asarray(flux_array)/(dx*dx)

                        # Note that the factor of 2 pi actually needs to be removed to get the right units.
                        # Currently, all the powers being plotted in this script are multiplied by 2 pi
                        # to make up for this factor
                        integrate_flux_X[ii] = (1/(2*pi))*trapz(trapz(flux_array[0, ...], dx=dx, axis=1), dx=dx)
                        integrate_flux_Y[ii] = (1/(2*pi))*trapz(trapz(flux_array[1, ...], dx=dx, axis=0), dx=dx)

                        integrate_power_X[ii] = integrate_flux_X[ii]*mu0
                        integrate_power_Y[ii] = integrate_flux_Y[ii]*mu1
                    except OSError:
                        print('Missing file')
                        print(input_file_name.format(E0, Ecouple, E1, mu0, mu1, num_minima0, minima1,
                                                     phase_shift))
                if abs(mu0) <= abs(mu1):
                    efficiency_ratio = -(integrate_power_X/integrate_power_Y)
                else:
                    efficiency_ratio = -(integrate_power_Y/integrate_power_X)

                with open(output_file_name.format(E0, E1, mu0, mu1, num_minima0, minima1, Ecouple), "w") as \
                        ofile:
                    for ii, phase_shift in enumerate(phase_array):
                        ofile.write(
                            f"{phase_shift:.15e}" + "\t"
                            + f"{integrate_flux_X[ii]:.15e}" + "\t"
                            + f"{integrate_flux_Y[ii]:.15e}" + "\t"
                            + f"{integrate_power_X[ii]:.15e}" + "\t"
                            + f"{integrate_power_Y[ii]:.15e}" + "\t"
                            + f"{efficiency_ratio[ii]:.15e}" + "\n")
                    ofile.flush()

def plot_power_efficiency_Ecouple(path,timescale, num_minima0, num_minima1, E0, E1, Ecouple, mu0, mu1):  # plot power and efficiency vs number of barriers n2

    output_file_name = (
            path + "/plots/" + "P_ATP_eff_Ecouple_" + "E0_{0}_Ecouple_{1}_E1_{2}_psi1_{3}_psi2_{4}_n1_{5}" + "_.pdf")
    f, axarr = plt.subplots(2, 1, sharex='all', sharey='none', figsize=(6, 8))

    # power plot
    axarr[0].axhline(0, color='black', linewidth=1)  # x-axis
    optimal_n2 = num_minima0/(mu1/mu0 + sqrt(1+(mu1/mu0)*(mu1/mu0)))
    axarr[0].axvline(optimal_n2, color='black', linestyle='--', linewidth=1)  # lining up features in the two plots
    axarr[0].fill_between([1, 250], 0, 31, facecolor='grey', alpha=0.2)  # shading power output
    '''
    # zero-barrier results
    input_file_name = (path + "/plotting_data/"
                       + "flux_zerobarrier_psi1_{0}_psi2_{1}_outfile.dat")
    data_array = loadtxt(input_file_name.format(mu0, mu1))
    num_minima1_array2 = array(data_array[:, 0])
    flux_y_array = array(data_array[:, 2])
    power_y = -flux_y_array * mu1
    axarr[0].plot(num_minima1_array2, 2*pi*power_y*timescale, '-', color='C0', label='$0$', linewidth=2)
    '''
    # Fokker-Planck results (barriers)
    i = 0  # only use phase=0 data
    power_y_array = []
    tight_power_array = empty(0) # curve at tight coupling limit
    for ii, minima1 in enumerate(num_minima1):
        input_file_name = (path + "/plots/" + "flux_power_efficiency_"
                           + "E0_{0}_E1_{1}_psi1_{2}_psi2_{3}_n1_{4}_n2_{5}_Ecouple_{6}" + "_outfile.dat")
        try:
            data_array = loadtxt(
                input_file_name.format(E0, E1, mu0, mu1, num_minima0, minima1, Ecouple),
                usecols=(0, 4))
            
            power_y = nparray(data_array[1])
            power_y_array = append(power_y_array, power_y)

            tight_power = -(mu1 + (num_minima0*mu0/minima1-mu1)/(1+(num_minima0/minima1)*(num_minima0/minima1)))*mu1
            tight_power_array = append(tight_power_array, tight_power)
        except OSError:
            print('Missing file flux')

    axarr[0].plot(num_minima1, -2.0*pi*power_y_array*timescale, 'o', color='C1', label=r'$E_{couple}=$'+str(Ecouple), markersize=8)
    axarr[0].plot(num_minima1, tight_power_array, 'o', color='black', label='Tight coupling limit', markersize=8)
    
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
    '''
    # zero-barrier curve
    input_file_name = (
            path + "/plotting_data/"
            + "flux_zerobarrier_psi1_{0}_psi2_{1}_outfile.dat")
   
    try:
        data_array = loadtxt(input_file_name.format(mu0, mu1))
        Ecouple_array2 = array(data_array[1:, 0])
        Ecouple_array2 = append(Ecouple_array2, 128.0)  # add point to end up with curves of equal length
        flux_x_array = array(data_array[1:, 1])
        flux_y_array = array(data_array[1:, 2])  # skip the point at zero, which is problematic on a log scale
        flux_x_array = append(flux_x_array, flux_x_array[-1])  # copy last point to add one
        flux_y_array = append(flux_y_array, flux_y_array[-1])
        axarr[1].plot(num_minima1_array2, flux_y_array / (flux_x_array), '-', color='C0', linewidth=2)
    except:
        print('Missing data efficiency')
    '''
    # Fokker-Planck results (barriers)
    eff_array = []
    tight_eff_array = empty(0)
    for ii, minima1 in enumerate(num_minima1):
        input_file_name = (
                path + "/plots/flux_power_efficiency_"
                + "E0_{0}_E1_{1}_psi1_{2}_psi2_{3}_n1_{4}_n2_{5}_Ecouple_{6}" + "_outfile.dat")
        try:
            data_array = loadtxt(
                input_file_name.format(E0, E1, mu0, mu1, num_minima0, minima1, Ecouple), usecols=5)

            eff_array = append(eff_array, data_array)
            tight_eff = -num_minima0*mu1/(minima1*mu0)
            tight_eff_array = append(tight_eff_array, tight_eff)
        except OSError:
            print('Missing file efficiency')
    axarr[1].plot(num_minima1, eff_array, 'o', color='C1', markersize=8)
    axarr[1].plot(num_minima1, tight_eff_array, 'o', color='black', markersize=8)

    axarr[1].set_xlabel(r'$n_2$', fontsize=20)
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
    f.savefig(output_file_name.format(E0, Ecouple, E1, mu0, mu1, num_minima0))

if __name__ == "__main__":

    path = os.getcwd()
    flux_power_efficiency(path,dx,N,num_minima0,num_minima1,phase,E0,E1,Ecouple,mu0,mu1)
    plot_power_efficiency_Ecouple(path,timescale, num_minima0, num_minima1, E0, E1, Ecouple, mu0, mu1)
