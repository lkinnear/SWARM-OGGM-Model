# SWARM OGGM Model
 A repository containing scripts for running OGGM for the SWARM project. These scripts require a modified version of the oggm source code that can be found at:  https://github.com/lkinnear/oggm and installed via the method detailed on the oggm website at: **[install the dev version + get access to the OGGM code:](https://docs.oggm.org/en/v1.1/installing-oggm.html)**.

To run oggm for the SWARM project please refer to the **[oggm website](https://docs.oggm.org/en/v1.1/run_examples/run_mb_calibration.html)** about running from pre-pro level 1 i.e. with your own climate data.
Running oggm will require mass balance calibration for the climate dataset used as detailed in the link above.

## Folder Contents
### 1. CORDEX_runs: workflows for the runs using CORDEX data
        1. test
 1. MSWEP_ERA_runs: workflows for the runs using MSWEP+ERA data
 1. Custom_Climate_runs: Needs renamed since this currently uses custom cliamte data but runs the OGGM forward model using run_random_climate
 1. data_extraction_scripts: scripts for extracting data from the OGGM run folders
 1. mb_calibration_scripts: scripts for running and testing the mass balance calibration
 1. shape_files: the various shape files used
 1. figure_scripts: scripts for producing shape files
 1. full_workflow: a folder containing full workflows (largely superseded and can probably be deleted)

 prepro_directory_tester.py a simple script to allow compaison between different prepro directory levels


