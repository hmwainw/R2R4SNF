This is the reactor-to-repository framework for evaluating spent nuclear fuel (SNF) from small modular reactors and other reactors.
It consists of includes (a) reactor physics models and post-discharge decay calculations, (b) repository footprint models, and (c) SNF repository performance assessment (PA) models with explicit radionuclide release/transport simulations.

![image](https://github.com/hmwainw/R2R4SNF/assets/110697247/6ef9d4b5-6382-4773-9d70-168b1de5759d)

- [OpenMC](openmc.org) first simulates reactor physics, including neutronics, depletion analysis, radionuclide generation, and decay-chain calculations during the reactor operation and post-discharge periods. It provides commonly used SNF metrics, including SNF mass (heavy metal equivalent) and volume as well as total activity, radiotoxicity, and decay heat.

- [NWPY](https://escholarship.org/uc/item/4n9157tz) then quantifies the repository area per package and per GWe.y given the thermal constraints.
  
![image](https://github.com/hmwainw/R2R4SNF/assets/110697247/9440e667-d0b2-4cf8-b4d4-dfd1ebefd15f)

- [PFLOTRAN](https://www.pflotran.org/) is then used for the generic repository perfrmance assessment model to compute the release of radionuclides from the waste forms, their migration in the geosphere, and the peak dose rates over one million years. The conceptural model is from [Stein et al. (2018)](https://www.osti.gov/servlets/purl/1513634)
  
![image](https://github.com/hmwainw/R2R4SNF/assets/110697247/775380db-e53e-4909-86a5-6a1726824866)


This repo includes:

0. OpenMC [input files](https://github.com/hmwainw/R2R4SNF/tree/master/OpenMC_input)
1. Post-processing of OpenMC to compute the SNF metrics (mass, volume, radiotoxicity, activity and decay heat): [1_Openmc_post.ipynb](https://github.com/hmwainw/R2R4SNF/blob/master/1_Openmc_post_processing.ipynb)
2. Repository footprint analysis with NWPY: [2_Repository_footprint.ipynb](https://github.com/hmwainw/R2R4SNF/blob/master/2_Repository_footprint_simulation.ipynb)
3. Repository PA with PFLOTRAN: [4_I129 release_transport.ipynb](https://github.com/hmwainw/R2R4SNF/blob/master/3_I129_release_transport.ipynb). Note that the simulations are run in [PFLOTRAN/](https://github.com/hmwainw/R2R4SNF/tree/master/Results/I129_release_transport). run_pf.py can be used to run multiple simulations in parallel.
4. Comparison of the overall statistics: [5_Overall_comparison.ipynb](https://github.com/hmwainw/R2R4SNF/blob/master/5_Overall_comparison.ipynb)

Note that these codes are not optimized for the speed/efficiency and kept simple for educational purposes. 




