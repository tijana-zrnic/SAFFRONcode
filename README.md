--------------------------- README --------------------

Plots saved as .pdf files in the folder "plots", data saved as .dat files in the folder "dat"

Note that the plots may look different than the ones in the paper because the observations are randomly generated

To reproduce the plots in the paper:

1. FDR and Power of SAFFRON and LORD vs pi1 for sequences of varying aggressiveness, Gaussian tests

for SAFFRON: python run_and_plot.py --FDRrange 2,4,5,6 --mu-gap 2 OR python run_and_plot.py --FDRrange 2,4,5,6 --mu-gap 3
for LORD: python run_and_plot.py --FDRrange 7,9,10,11 --mu-gap 2 OR python run_and_plot.py --FDRrange 7,10,11,12 --mu-gap 3

2. FDR and Power of SAFFRON, LORD and Alpha-investing vs pi1, Gaussian tests

python run_and_plot.py --FDRrange 2,9,13 --mu-gap 2 OR python run_and_plot.py --FDRrange 2,12,13 --mu-gap 3

3. FDR and Power of SAFFRON and LORD vs pi1 for sequences of varying aggressiveness, beta alternatives (don't forget to set mu-gap to 0)

for SAFFRON: python run_and_plot.py --FDRrange 1,3,5,6 --mu-gap 0 --mod-choice 2
for LORD: python run_and_plot.py --FDRrange 8,10,11,12 --mu-gap 0 --mod-choice 2

4. FDR and Power of SAFFRON, LORD and Alpha-investing vs pi1, beta alternatives (don't forget to set mu-gap to 0)

python run_and_plot.py --FDRrange 3,10,13 --mu-gap 0 --mod-choice 2

5. FDR and Power of alpha-investing and SAFFRON alpha-investing vs pi1, Gaussian tests

python run_and_plot.py --FDRrange 13,14 --mu-gap 3

6. Test level alpha_k of alpha-investing and SAFFRON alpha-investing vs hypothesis index, Gaussian tests

python run_and_plot.py --FDRrange 13,14 --mu-gap 3 --plot-style 0


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

If you spot any issues or bugs, please contact me at tijana(at)eecs(dot)berkeley(dot)edu

This code borrowed substantial parts from Fanny Yang's code available at: https://github.com/fanny-yang/OnlineFDRCode


