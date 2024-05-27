# Bias correction of ICAR CMIP data

This document describes the Bias correction process of ICAR data. This is a more generic version, that can be used to bias correct to any observational dataset (e.g. Livneh, GMET, ....).
A livneh-specific version can be found at https://github.com/bertjebertjek/BC2Liv


### Prerequisites
 - Postprocess the raw (or 3hr raw) ICAR output first: dir `CMIP_ICAR_postprocess` (on github: `https://github.com/bertjebertjek/CMIP_ICAR_postprocess`)
    This does:
        - correction of negative precipitation due to restart errors
        - aggregation into monthly (yearly) files w 3hr (24hr) timestep
        - removes the GCM's convective precipitation (cp) from the ICAR precipitation
        - removes unwanted output variables (optional)


## 1. Regrid the post-processed ICAR data to the output grid:
in subdir `1_Regrid`
This is required before we can bias correct to the output grid.
Set in- and output paths, models, scenarios, and timestep (daily or 3hr) in  `submit_regrid2outgrid.sh`
Set path to observational files (bc_grid_files) in `regrid2outgrid.py` ln38 and a name for the grid in ln39 (this will be used in the file naming.)



## 2. The Bias Correction

Bias correction is split into separate scripts /procedures for precipitation and temperature, because the runtime may otherwise exceed 12 hrs for large domains/small timesteps: `BC_Icar_5y_pcp.py` and `BC_Icar_5y_ta2m.py`

These script use Quantile mapping (in separate functions/files) with historical CMIPX-ICAR simulation as reference data, and Livneh as observations. QM is done in 5 year blocks of input data, where this 5 year period is excluded from the reference data, in order to preserve future extremes.

These functions can be called separately for each model & scenario from submit_BC_CMIP{5/6}_{pcp/ta2m}.sh, or for all models and scenarios from `loop_submit.sh`. This loop_submit script sets the correct model/scenario parameters by modifying a pbs template script (located in dir `auto_submit` ), and then launches 2 jobs per model/scenario combination; one for the precipitation bias correction, and a second job, dependent on the first, for the temperature bias correction.

Runtimes and memory requirements are heavily dependent on domain size and temporal resolution. For example, a small domain in the PNW at daily timestep requires 1h and 10GB per model/scenario, whereas the western US at 3hr and 6km resolution requires 12h and ~350 GB for each model/scenario combination.
Make sure to test and set appropriate memory use and time allocation in `2_Bias_Correct/auto_submit/template_submit_BC.sh`


### Further info: Precipitation BC
- There is a flag to add noise to the input and reference data during the bias correction procedure in `BC_Icar_5y_pcp.py`. Default=True


### Further info: Temperature BC
- relative humidity can be added as a variable to the output dataset calc_relhum. Default is True.
- drop_vars is a flag to drop variables listed in vars_to_drop. Default is true, but ideally the unwanted output variables have been removed in the post processing procedures `CMIP_ICAR_postprocess` (on github: `https://github.com/bertjebertjek/CMIP_ICAR_postprocess`) already.



### Post processing
- If there are still variables in the output we don;t need, use `remove_vars_and_float.py` to remove them and make all other variables float32 iso float64. If the preprocessing was done with the procedures in dir `../CMIP_ICAR_postprocess` this should not be necessary.

- the Bias correction produces subfolders 3hr_pcp (daily_pcp) and 3hr (daily). The *_pcp folders only have the precipitation corrected, and serve as input for the temperature bias correction. These should be deleted after checking that all bias correction has terminated succesfully, to avoid confusion and save space.




_Bert Kruyt , NCAR RAL, february 2024_

