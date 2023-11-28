# Thermal_trait_parameterization
Phytoplankton thermal trait parameterization alters community structure and biogeochemical processes in a modeled ocean 

These scripts are provided in the interests of open science. If you have questions or find errors, please let us know.

Contact:<br/>
Stephanie I. Anderson<br/>
Department of Earth, Atmospheric and Planetary Sciences<br/>
Massachusetts Institute of Technology<br/>
siander@mit.edu<br/>


## Directory structure
- [data](data/): has the raw data
- [output](output/): files output by scripts 
- [scripts](scripts/): run using data and output, produce figures and tables
- [figures](figures/):  figures


## Before beginning, model output from the Darwin simulations must be obtained separately here:
## INCLUDE!!!!!!!!!!!!


## Analysis  workflow brief overview
- [Thermal_dependency_curves.R](scripts/Thermal_dependency_curves.R) to fit exponential curves to Phytoplankton Functional Type (PFT) growth data
- All files beginning with ‘q10’ to assess model output from the Darwin simulations
- [PHYSAT_MAREDAT.py](scripts/PHYSAT_MAREDAT.py) to compare with field and satellite observations
- [Figure8_simulation.ipynb](scripts/Figure8_simulation.ipynb) for a simplified 


## Script Overview

1. [Thermal_dependency_curves.R](scripts/Thermal_dependency_curves.R)
    1. Contains:
	1. Table 1. Model Parameters
	2. Table S1. a_PCmax
	3. Figure 1. Thermal Dependency curves for each model and PFT
	4. Figure S2 & S3. exponential curve fits for Kremer (S2) and Anderson (S3)
	5. Figure S4. allometrically scaled maximum growth curves
    2. Requires
       	 1. Phytoplankton strain information from [derived_traits_diazogreen.csv](data/derived_traits_diazogreen.csv).
	2. Phytoplankton growth rates from [growth_rates_diazogreen.csv](data/growth_rates_diazogreen.csv).
	3. Parameters from Anderson et al. 2021 ([data/Anderson_etal_2021_table1.csv](data/data/Anderson_etal_2021_table1.csv)).

2. [q10_PP.py](scripts/q10_PP.py)
    1. Contains:
	1. Figure 4. Mean biomass for each PFT from each model
	2. Figure S5. Difference in primary production (PP) and biomass between control simulations
	3. Figure 2. Zonal mean biomass and PP for each model
    2. Requires
       	1. Model grid ([grid_igsm.nc](data/grid_igsm.nc)).
	2. HARVARD DATAVERSE

3. [q10_diversity.py](scripts/q10_diversity.py)
    1. Contains:
	1. Figure 3 A&B. Phytoplankton community evenness and richness
    2. Requires
       	1. Model grid ([grid_igsm.nc](data/grid_igsm.nc)).
	2. HARVARD DATAVERSE

4. [q10_dominantPFT.py](scripts/q10_dominantPFT.py)
    1. Contains:
	1. Figure 3 C. Global maps of the dominant PFT for each model.
    2. Requires
       	1. Model grid ([grid_igsm.nc](data/grid_igsm.nc)).
	2. HARVARD DATAVERSE

5. [q10_export.py](scripts/q10_export.py)
    1. Contains:
	1. Figure S10. Mean biovolume
	2. Figure 5. Export processes
    2. Requires
       	1. Model grid ([grid_igsm.nc](data/grid_igsm.nc)).
	2. HARVARD DATAVERSE
	3. PFT biomass with climate change ([biomass_SQDUall.txt](output/biomass_SQDUall.txt) and [biomass_DQDUall.txt](output/biomass_DQDUall.txt))

6. [q10_bray_curtis.py](scripts/q10_bray_curtis.py)
    1. Contains:
	1. Figure 6. Bray-Curtis dissimilarity
	2. Figure 7. Mean change in depth-averaged PFT biomass
	3. Figure S6. Mean change in each PFT biomass
	4. PFT biomass with climate change ([biomass_SQDUall.txt](output/biomass_SQDUall.txt) and [biomass_DQDUall.txt](output/biomass_DQDUall.txt))
    2. Requires
       	1. Model grid ([grid_igsm.nc](data/grid_igsm.nc)).
	2. HARVARD DATAVERSE

7. [q10_warming_only.py](scripts/q10_warming_only.py)
    1. Contains:
	1. Figure S9. Primary production temp only vs. all physics changing
	2. Figure S11. Bray-Curtis dissimilarity
    2. Requires
       	1. Model grid ([grid_igsm.nc](data/grid_igsm.nc)).
	2. HARVARD DATAVERSE
	3. Requires PFT biomass from control simulations in [q10_bray_curtis.py](scripts/q10_bray_curtis.py)

8. [PHYSAT_MAREDAT.py](scripts/PHYSAT_MAREDAT.py)
    1. Contains:
	1. Figure S7. PHYSAT
	2. Figure S8. MAREDAT
    2. Requires
       	1. PHYSAT data ([physat_monthly_mapped_degree.nc](data/physat_monthly_mapped_degree.nc)).
	2. MAREDAT data downloaded from [PANGAEA](http://www.pangaea.de/search?&q=maredat )

9. [Figure8_simulation.ipynb](scripts/Figure8_simulation.ipynb)
    1. Contains:
	1. Figure 8
    2. Requires
       	1. Table S3:  ([ModelParameters.csv](data/ModelParameters.csv)).
