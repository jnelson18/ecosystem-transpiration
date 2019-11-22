# ecosystem-transpiration
Code and examples of how to estimate transpiration from eddy covariance data.

This is a companion repository for a soon to be published manuscript on evapotranspiration partitioning of eddy covariance data. As the manuscript is yet to be published, this page will be updated in the coming months. The paper utilized three partitioning methods:

- Perez-Priego et al (2018). Partitioning Eddy Covariance Water Flux Components Using Physiological and Micrometeorological Approaches. Journal of Geophysical Research: Biogeosciences. https://doi.org/10.1029/2018JG004637

- Nelson et al (2018). Coupling Water and Carbon Fluxes to Constrain Estimates of Transpiration: The TEA Algorithm. Journal of Geophysical Research: Biogeosciences. https://doi.org/10.1029/2018JG004727

- Zhou et al (2016). Partitioning evapotranspiration based on the concept of underlying water use efficiency: ET PARTITIONING. Water Resources Research, 52(2), 1160–1175. https://doi.org/10.1002/2015WR017766

To find the code for each method:

- the Perez-Priego method is implemented in R [here](https://github.com/oscarperezpriego/ETpartitioning)
- the TEA metod (Nelson 2018) is implemented as a python package [here](https://github.com/jnelson18/TranspirationEstimationAlgorithm)
- the Zhou method is implemented in python in this repository as [zhou.py](zhou.py)

## Tutorials

Tutorials for the TEA and Zhou methods are contianed in this repository. Each can be run within a browser without any installation:

- TEA [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/jnelson18/ecosystem-transpiration/master?filepath=TEA_tutorial.ipynb)
- Zhou [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/jnelson18/ecosystem-transpiration/master?filepath=Zhou_tutorial.ipynb)


The tutorials can also be run directly using the anaconda environment defined in [environment.yml](environment.yml)
