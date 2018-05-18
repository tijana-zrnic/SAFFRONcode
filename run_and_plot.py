import logging, argparse
import numpy as np
from exp_FDR_batch_new import*
from plot_batch_results import*
from toimport import *


def main():

    if not os.path.exists('./dat'):
        os.makedirs('./dat')

    #########%%%%%%  SET PARAMETERS FOR RUNNING EXPERIMENT %%%%%%%##########

    FDRrange = str2list(args.FDRrange)
    pirange = str2list(args.pirange, 'float')
    mu_gap = args.mu_gap
    markov_lag_list = str2list(args.markov_lag, 'int')
    hyprange = [0]

    ########%%%%%%%%%%%%%%%%% RUN EXPERIMENT %%%%%%%%########################
    
    for pi in pirange:
        # Run single FDR
        for FDR in FDRrange:
            for markov_lag in markov_lag_list:
                # Prevent from running if data already exists
                filename_pre = 'MG%.1f_Si1.0_FDR%d_NH%d_ND%d_L%d_PM%.2f_NR%d_MOD%d' % (args.mu_gap, FDR, args.num_hyp, 1, markov_lag, pi, args.num_runs, args.mod_choice)
                all_filenames = [filename for filename in os.listdir('./dat') if filename.startswith(filename_pre)]
                # all_filenames = []

                # Run experiment if data doesn't exist yet
                if all_filenames == []:
                    print("Running experiment for FDR procedure %s and pi %.1f with lag %d" % (proc_list[FDR], pi, markov_lag))
                    run_single(args.num_runs, args.num_hyp, 1, mu_gap, pi, args.alpha0, markov_lag, args.mod_choice, FDR, sigma = 1, verbose = False)
                else:
                    print("Experiments for FDR procedure %s and pi %.1f with lag %d are already run" % (proc_list[FDR], pi, markov_lag))

    # Plot different measures over hypotheses for different FDR
    print("Now plotting ... ")
    plot_results(args.plot_style, 0, FDRrange, pirange, hyprange, mu_gap, 1, args.num_hyp, args.num_runs, markov_lag_list, args.mod_choice)
            


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--FDRrange', type=str, default = '5') # choice of algorithms and parameters (listed in exp_FDR_batch_new)
    parser.add_argument('--num-runs', type=int, default = 200) # number of independent trials
    parser.add_argument('--num-hyp', type=int, default = 1000) # number of hypotheses
    parser.add_argument('--plot-style', type = int, default = 1) # 0 for plots vs hyp, 1 for plots vs pi1
    parser.add_argument('--alpha0', type=float, default = 0.05) # test level
    parser.add_argument('--mu-gap', type=float, default = 3) # mu_c for gaussian tests
    parser.add_argument('--mod-choice', type=int, default = 1) # 1 for gaussian tests, 2 for beta alternatives (*set mu-gap=0 when not doing gaussian tests*)
    parser.add_argument('--pirange', type=str, default = '0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9') # range of pi1
    parser.add_argument('--markov-lag', type=str, default='0,1,2,3')
    args = parser.parse_args()
    logging.info(args)
    main()
