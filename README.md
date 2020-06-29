[![DOI](https://zenodo.org/badge/223398945.svg)](https://zenodo.org/badge/latestdoi/223398945)

# ecosystem-transpiration

Code and examples of how to estimate transpiration from eddy covariance data.

This is a companion repository to Nelson et al (2020) which compares three evapotranspiration partitioning methods to estimate transpiration from eddy covariance data. Each method is fully described in the following reference manuscripts:

- **Pérez-Priego**: Perez-Priego et al (2018). Partitioning Eddy Covariance Water Flux Components Using Physiological and Micrometeorological Approaches. Journal of Geophysical Research: Biogeosciences. https://doi.org/10.1029/2018JG004637

- **Nelson**: Nelson et al (2018). Coupling Water and Carbon Fluxes to Constrain Estimates of Transpiration: The TEA Algorithm. Journal of Geophysical Research: Biogeosciences. https://doi.org/10.1029/2018JG004727

- **uWUE**: Zhou et al (2016). Partitioning evapotranspiration based on the concept of underlying water use efficiency: ET PARTITIONING. Water Resources Research, 52(2), 1160–1175. https://doi.org/10.1002/2015WR017766

## Code

The code for each method can be found:

- the Pérez-Priego method is implemented in R [here](https://github.com/oscarperezpriego/ETpartitioning)
- the TEA metod (Nelson 2018) is implemented as a python package [here](https://github.com/jnelson18/TranspirationEstimationAlgorithm)
- the uWUE method is implemented in python in this repository as [zhou.py](zhou.py)

## Tutorials

Tutorials for the TEA, uWUE, and Pérez-Priego methods are contianed in this repository as jupyter notebooks. Each can be run within a browser without any installation via Binder:

- TEA [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/jnelson18/ecosystem-transpiration/master?filepath=TEA_tutorial.ipynb)
- Zhou [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/jnelson18/ecosystem-transpiration/master?filepath=Zhou_tutorial.ipynb)
- Pérez-Priego [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/jnelson18/ecosystem-transpiration/master?filepath=Perez-Priego_tutorial.ipynb)


All code required to run the partitionin methods, as well as the tutorials, can be installed using the aAnaconda environment defined in [environment.yml](environment.yml). Instructions for installing Anaconda can be found [here](https://docs.anaconda.com/anaconda/install/).

To install the Anaconda environment, which includes software required to run the three pratitioning methods, first download the repository. If using Anaconda Navigator, follow the [Importing an environment instructions](https://docs.anaconda.com/anaconda/navigator/tutorials/manage-environments/#importing-an-environment). If using the command line, follow [Creating an environment from an environment.yml file](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file).

Once installed, activate the environment. The tutorials can then be run locally, for example:

```
jupyter notebook TEA_tutorial.ipynb
```

In order to run the Pérez-Priego method, the R packages found in [install.R](install.R) much be installed before the tutroial will run. This can be installed from the command line with:

```
Rscript install.R
```

