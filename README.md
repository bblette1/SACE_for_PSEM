# SACE_for_PSEM
Code implementing simulations and application in the manuscript "Post-randomization Biomarker Effect Modification in an HIV Vaccine Clinical Trial" https://arxiv.org/abs/1811.03930

### Reproducing main results from the paper

Simulations:

* Run main_NEE.R (or main_NEH.R)

* Rdata files containing the results will be put in the working directory and named according to the number of sims specified in "nsims <- __"

* Run make_plots_NEE.R (or make_plots_NEH.R), changing the names of the loaded Rdata files to match the previous step

Application:

* Run main_app.R, making sure that a .csv file of the data is in the working directory

### Using methods on other data

* Think about which assumptions your data meet and run analyze_NEE.R, analyze_NEH.R, or analyze_NEB_app.R to obtain point estimates and ignorance intervals

* To obtain EUIs, use output from above and run code chunks labeled as variance estimation in main_NEE.R, main_NEH.R, or (for NEB) main_app.R