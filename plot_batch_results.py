### Import Python libraries
import numpy as np
np.set_printoptions(precision = 4)

import sys


### Import utilities for plotting
from plotting import*
from settings_util import*
from toimport import*

# Plot_styles:
# 0: wealth, alpha_k over time
# 1: power/FDR over pi1


def plot_results(plot_style, whichrun, FDRrange, pirange, hyprange, mu_gap, sigma, NUMHYP, num_runs, markov_lag_list, mod_choice, NUMDRAWS = 1):

    plot_dirname = './plots'
    numrun = 100000

    #%%%%%%%%%%%%%%%%%%%%  PLOTS vs. Hyp (time)  %%%%%%%%%%%%%%%%%%%%%%

    if plot_style == 0:

        numFDR = len(FDRrange)
        
        pi = pirange[0]

        # ----------- LOAD DATA --------
        FDR_mat = [None]*len(FDRrange)
        wealth_mat = [None]*len(FDRrange)
        TDR_mat = [None]*len(FDRrange)
        alpha_mat = [None]*len(FDRrange)

        for FDR_j, FDR in enumerate(FDRrange):

            for markov_lag in markov_lag_list:
                filename_pre = 'MG%.1f_Si%.1f_FDR%d_NH%d_ND%d_L%d_PM%.2f_NR%d_L%d_MOD%d' % (mu_gap, sigma, FDR, NUMHYP, NUMDRAWS, markov_lag, pi, num_runs, mod_choice)
                all_filenames = [filename for filename in os.listdir('./dat') if filename.startswith(filename_pre)]

                if all_filenames == []:
                    print("No files found!")
                    print(filename_pre)
                    sys.exit()

                # Load results
                result_mat = np.loadtxt('./dat/%s' % all_filenames[0])
                FDR_vec = result_mat[0:NUMHYP, whichrun]
                rej_vec = result_mat[1*NUMHYP:2*NUMHYP, whichrun]
                falrej_vec = result_mat[2*NUMHYP:3*NUMHYP, whichrun]
                wealth_vec = result_mat[3*NUMHYP:4*NUMHYP, whichrun]
                alpha_vec = result_mat[5*NUMHYP:6*NUMHYP, whichrun]

                # Get true Hypo vector
                Hypo = get_hyp(pi, NUMHYP)
                Hypo = Hypo.astype(int)

                # Save to matrix
                FDR_mat[FDR_j] = FDR_vec
                wealth_mat[FDR_j] = wealth_vec
                # TDR_mat[FDR_j] = TDR_vec
                alpha_mat[FDR_j] = alpha_vec
        
        # -------- PLOT ---------------
        # Set x axis
        if len(hyprange) == 1:
            xs = range(NUMHYP)
            hyplen = NUMHYP
        else:
            # Cut the matrices
            xs = hyprange
            hyplen = len(hyprange)
            FDR_mat = np.array(FDR_mat)[:,0:len(hyprange)]
            wealth_mat = np.array(wealth_mat)[:,0:len(hyprange)]
            alpha_mat = np.array(alpha_mat)[:,0:len(hyprange)]
        
        legends_list = np.array(proc_list).take(FDRrange)[0:numFDR]


        leg_col = 1

        #### Wealth vs HYP ####
        filename = 'WealthvsHP_Plot%d_MG%.1f_Si%.1f_NH%d_ND%d_PM%.2f_HR%d_R%d_L%d' %  (plot_style, mu_gap, sigma, NUMHYP, NUMDRAWS, pi, hyplen, whichrun, markov_lag)
        plot_curves_mat(xs, wealth_mat, legends_list, plot_dirname, filename,  'Hypothesis index', 'Wealth($J$)', 0, leg_col = leg_col)

        #### alpha vs. HYP ####
        filename = 'alphavsHP_MG%.1f_Si%.1f_NH%d_ND%d_PM%.2f_HR%d_R%d_L%d' %  (mu_gap, sigma, NUMHYP, NUMDRAWS, pi, hyplen, whichrun, markov_lag)
        plot_curves_mat(xs, alpha_mat, legends_list, plot_dirname, filename,  'Hypothesis index', '$alpha(J)$', 0, leg_col = leg_col)

        #%%%%%%%%%%%%%%%%%%%  PLOTS vs. pi1 %%%%%%%%%%%%%%%%%%%%%%%%%%
        
    elif plot_style == 1:

        TDR_av = []
        TDR_std = []
        FDR_av = []
        FDR_std = []
        ind = 0
        
        for index, FDR in enumerate(FDRrange):

            for markov_lag in markov_lag_list:
                filename_pre = 'MG%.1f_Si%.1f_FDR%d_NH%d_ND%d_L%d' % (mu_gap, sigma, FDR, NUMHYP, NUMDRAWS, markov_lag)
                all_filenames = [filename for filename in os.listdir('./dat') if filename.startswith(filename_pre)]

                if all_filenames == []:
                    print("No file found!")
                    print(filename_pre)
                    sys.exit()

                # Get different pis
                pos_PM_start = [all_filenames[i].index('PM') for i in range(len(all_filenames))]
                pos_PM_end = [all_filenames[i].index('_NR') for i in range(len(all_filenames))]
                PM_vec = [float(all_filenames[i][pos_PM_start[i] + 2:pos_PM_end[i]]) for i in range(len(all_filenames))]

                order = np.argsort(PM_vec)
                PM_list = sorted(set(np.array(PM_vec)[order]))

                # Initialize result matrices
                TDR_av.append(np.zeros([1, len(PM_list)]))
                TDR_std.append(np.zeros([1, len(PM_list)]))
                FDR_av.append(np.zeros([1, len(PM_list)]))
                FDR_std.append(np.zeros([1, len(PM_list)]))
                TDR_vec = np.zeros(len(PM_list))
                FDR_vec = np.zeros(len(PM_list))
                TDR_vec_std = np.zeros(len(PM_list))
                FDR_vec_std = np.zeros(len(PM_list))

                # Merge everything with the same NA and NH
                for k, PM in enumerate(PM_list):
                    indices = np.where(np.array(PM_vec) == PM)[0]
                    result_mat = []
                    # Load resultmats and append
                    for j, idx in enumerate(indices):
                        result_mat_cache = np.loadtxt('./dat/%s' % all_filenames[idx])
                        result_mat_cache = result_mat_cache[6*NUMHYP:6*NUMHYP+2,0:200]
                        if (j == 0):
                            result_mat = result_mat_cache
                        else:
                            result_mat = np.c_[result_mat, result_mat_cache]

                    numrun = len(result_mat[0])
                    # Get first vector for TDR
                    TDR_vec[k] = np.average(result_mat[0])
                    TDR_vec_std[k] = np.true_divide(np.std(result_mat[0]),np.sqrt(numrun))
                    # FDR
                    FDR_vec[k] = np.average(result_mat[1])
                    FDR_vec_std[k] = np.true_divide(np.std(result_mat[1]), np.sqrt(numrun))
                TDR_av[ind] = [TDR_vec[k] for k in range(len(PM_list))]
                TDR_std[ind] = [TDR_vec_std[k] for k in range(len(PM_list))]
                FDR_av[ind] = [FDR_vec[k] for k in range(len(PM_list))]
                FDR_std[ind] = [FDR_vec_std[k] for k in range(len(PM_list))]

                ind = ind + 1



        # -------- PLOT ---------------
        xs = PM_list
        x_label = '$\pi_1$'

        # Create legend
        legends_list = np.array(proc_list).take(FDRrange)
        # legends_list = ['L=0','L=1','L=2','L=3']

        ##### FDR vs pi #####

        filename = 'FDRvsPI_MG%.1f_Si%.1f_NH%d_ND%d_L%d_MOD%d' %  (mu_gap, sigma, NUMHYP, NUMDRAWS, markov_lag_list[-1], mod_choice)
        plot_errors_mat(xs, FDR_av, FDR_std, legends_list, plot_dirname, filename, x_label, 'FDR')

        ##### TDR vs pi ####
        filename = 'PowervsPI_MG%.1f_Si%.1f_NH%d_ND%d_L%d_MOD%d' %  (mu_gap, sigma, NUMHYP, NUMDRAWS, markov_lag_list[-1], mod_choice)
        plot_errors_mat(xs, TDR_av, TDR_std, legends_list, plot_dirname, filename, x_label, 'Power')


