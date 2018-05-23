--------------------------- README --------------------

Plots saved as .pdf files in the folder "plots", data saved as .dat files in the folder "dat"

Note that the plots may look different than the ones in the paper because the observations are randomly generated

--------------------------------------------------------

The main file is run_and_plot.py. The experiments vary depending on the following passed arguments:
--FDRrange - integers encoding the choice of algorithms and parameters (listed in exp_FDR_batch_new)
--num-runs - number of independent trials
--num-hyp - number of hypotheses
--plot-style - 0 for plots vs hyp, 1 for plots vs pi1
--alpha0 - test level
--mu-gap - used for gaussian tests as mu_c, where observations under the alternative are N(Z,1), Z~N(mu_c,1)
--mod-choice - 1 for gaussian tests, 2 for beta alternatives (*set mu-gap=0 when not doing gaussian tests*)
--pirange - range of pi1
--markov-lag - parameter of local dependence between p-values

If you spot any issues or bugs, please contact me at tijana(at)eecs(dot)berkeley(dot)edu

This code borrowed substantial parts from Fanny Yang's code available at: https://github.com/fanny-yang/OnlineFDRCode


