README

This repository includes different python algorithms that allow the user to create and analyse ocean temperature related time series from loaded data files. 
This project was part of a BSc thesis under the supervision of Spanish Ocean Institute (Instituto Español de Oceanografía, IEO). The aim of this thesis was 
to expand the work of Somavilla et al. (2013), where the authors created and analyzed time series from IEO buoy data for a 15 year span in the Bay of Biscay, 
subsequently studying the dynamics of heat and salinity balance processes and the behaviour of the Mixed Layer. Some of the cells in the main Jupyter file 
are improved python implementions of Matlab modules that were originally used in Somavilla et al. (2013); which in turn belong to the Matlab air-sea toolbox 2.0
(https://github.com/sea-mat/air-sea/blob/master/hfbulktc.m). In this main jupyter file, all the modules were self developed, even if other python implementations
of the same functions exist in different repositories. 

As for the other modules in this repository:

- atmosphere.py, windstress.py and constants.py are part of the Python version of Matlab air-sea toolbox 2.0 (Filipe P. A. Fernandes):
    - atmosphere.py : https://github.com/pyoceans/python-airsea/blob/master/airsea/atmosphere.py
    - windstress.py : https://github.com/pyoceans/python-airsea/blob/master/airsea/windstress.py
    - constants.py : https://github.com/pyoceans/python-airsea/blob/master/airsea/constants.py
      
- functions_for_fluxes.py : contains self-developed functions independently adapted from the Matlab air-sea toolbox 2.0
  (https://github.com/sea-mat/air-sea/blob/master/hfbulktc.m); some of them were not found in the previous modules, such as LKB, psiutc and psittc.
  
- Fourier_series_1.py: this is also a custom module, useful for analyzing a time series by creating its Fourier spectrum. AI-assisted generation.
