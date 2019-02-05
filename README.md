# SACE_for_PSEM
Code implementing simulations and application in the manuscript "Post-randomization Biomarker Effect Modification in an HIV Vaccine Clinical Trial" https://arxiv.org/abs/1811.03930

### Reproducing main results from the paper

Simulations:

* Run main_NETE.R (or main_NHETE.R)

* Rdata files containing the results will be put in the working directory and named according to the number of sims specified in "nsims <- __"

* Run make_plots_NETE.R (or make_plots_NHETE.R), changing the names of the loaded Rdata files to match the previous step

Application:

* Run main_app.R, making sure that a .csv file of the data is in the working directory