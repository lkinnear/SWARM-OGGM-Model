# SWARM OGGM Model
 A repository containing scripts for running OGGM for the SWARM project. These scripts require a modified version of the oggm source code that can be found at:  https://github.com/lkinnear/oggm and installed via the method detailed on the oggm website at: **[install the dev version + get access to the OGGM code:](https://docs.oggm.org/en/v1.1/installing-oggm.html)**.

To run oggm for the SWARM project please refer to the **[oggm website](https://docs.oggm.org/en/v1.1/run_examples/run_mb_calibration.html)** about running from pre-pro level 1 i.e. with your own climate data.
Running oggm will require mass balance calibration for the climate dataset used as detailed in the link above and file structures will need to be modified based on the user's directories.

## Folder Contents
### CORDEX_runs: 
        Workflows for the runs using CORDEX data
        1. 'full_workflow_CORDEX_gcm.py' A full workflow as a template for running using gcm data (along with the reference period data)
        1. 'full_workflow_CORDEX_HADGEM.py' Workflow for the HadGEM2 data
        1. 'full_workflow_CORDEX_MPI.py' Workflow for the MPI data
        1. 'full_workflow_CORDEX_NCC.py' Workflow for the NCC data
        1. 'full_workflow_CORDEX.py' Old workflow to be deleted in future

### MSWEP_ERA_runs 
        Workflows for the runs using MSWEP+ERA data*
        1. 'full_workflow_mswep_era_cru_data.py' A full workflow for running the reference period using the CRU dataset (for comparison)
        1. 'full_workflow_mswep_era_run.py' A full workflow for the reference perdion using the MSWEP/ERA data

### Data_extraction_scripts 
        Scripts for extracting data from the OGGM run output folders
        1. 'data_extraction_for_fuse_gcm_data.py' Data extraction for gcm data runs
        1. 'data_extraction for fuse.py' Data extraction for generic runs extracting from *climate_historical.nc* files
        1. 'data_extraction_from_per_glacier.py' Old testing way for extracting data, can be deleted
###  mb_calibration_scripts 
        Scripts for running and testing the mass balance calibration
        1. 'cross_val.py' Script for carrying out the cross validation of massbalance calibration [from the oggm website](https://docs.oggm.org/en/v1.1/run_examples/run_mb_calibration.html#cross-validation)
        1. 'mass_balance_regional_ERA5_MSWEP_reference_period.py' Mass balance calibration for the reference period using the MSWEP/ERA data
### shape_files 
    The various shape files used
    1. 'gadm36_CHN_0.shp' Outline of China_test
    1. 'RiverBasins.shp' The river basins (as separate basins)
    1. 'RiverBasinsMerged.shp' The river basins (merged together to allow searching within using geopandas)
    1. 'swarm_grid.shp' The outline of the swarm grid
### figure_scripts 
    scripts for producing shape files
    1. 'basic_raster_plot.py' Initial script for plotting the output data as a raster


### prepro_directory_tester.py 
    A simple script to allow comparison between different prepro directory levels


