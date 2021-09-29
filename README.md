# CLSNA-2Party-Polarization
This repository contains code and data supporting the article "Disentangling positive and negative partisanship in affective polarization using a coevolving latent space network with attractors model": 

[https://arxiv.org/abs/2109.13129](https://arxiv.org/abs/2109.13129)

## Data

Data used for analysis are 1) Twitter congressional hashtag networks; 2) Reddit comment networks. 

They are stored as R binary files (.RData). Data will be publicly available later, upon publication, in `Data/`. 

## Code
The corresponding code to conduct both simulation and analyses on real data sets are provided in `Code/`.

* `Code/model/`: code for implementations of model and inference with no change-point.

* `Code/model_change_point/`: code for implementations of model and inference with multiple change-points.
* The main files for simulation study with synthetic data are `Code/model/main_simulation.R` and `Code/model_change_point/main_simulation_mc.R`.

* **Figure 1,2**: `Code/app_plot_density.R`.
* **Figure 4**:`Code/sim_nets_traces.R`.
* **Figure 5, 7**, **Table 1,2**: `Code/model/main_application.R`, `Code/model/posterior_plot_app.R`.
* **Figure 6, 8**: `Code/model_change_point/main_application_mc.R`, `Code/model_change_point/plot_sub_window_1cp.R`.

If there are any bugs in the code, please contact Xiaojing Zhu at xiaojzhu@bu.edu
