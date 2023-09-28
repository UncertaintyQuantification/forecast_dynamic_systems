# forecast_dynamical_systems
Code for "Probabilistic forecast of nonlinear dynamical systems with uncertainty quantification"

Mengyang Gu, Yizi Lin, Victor Chang Lee, and Diana Y. Qiu (2023)


This folder contains data and code reproducing the results from the paper.

data: The folder contains the simulated data of quantum many-body systems far from equilib-
rium, which is associated with Example 2: Time-dependent Green's function.

functions: The folder contains the code to implement the Dynamic Mode Decomposition (DMD) algorithm.

eg1_lorenz96_fig1_tab2.R: The code for Fig.1 (forecast of the Lorenz 96 system for 900 steps by AR(1), DMD, HODMD, and PPGP) and Table 2 (forecast accuracy and uncertainty assessment on the held-out data).

eg1_lorenz96_fig2.R: The code for generating the forecast of PP-GP in Lorenz 96 example (Fig.2).

eg1_lorenz96_fig3.R: The code for generating the predictive standard deviation and cumulative mean absolute error from PP-GP in Lorenz 96 example (Fig.3).


eg2_time_dependent_green_function_tab3_fig4.R: The code for comparing the forecast and uncertainty assessment by AR(1), DMD, HODMD, and PP-GP (Table 3, Fig.4). 2500-time steps are used as training data.

eg2_time_dependent_green_function_tab4_fig5.R: The code for comparing the forecast and uncertainty assessment by AR(1), DMD, HODMD, and PP-GP (Table 4, Fig.5). 3500-time steps are used as training data.
