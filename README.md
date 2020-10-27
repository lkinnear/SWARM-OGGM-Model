# SWARM OGGM Model
 A repository containing scripts for running OGGM for the SWARM project. These scripts require a modified version of the oggm source code that can be found at:  https://github.com/lkinnear/oggm and installed via the method detailed at **[install the dev version + get access to the OGGM code:](https://docs.oggm.org/en/v1.1/installing-oggm.html)**. 
 Currently contains:
 1. CORDEX_runs: workflows for the runs using CORDEX data
 2. MSWEP_ERA_runs: workflows for the runs using MSWEP+ERA data
 3. Custom_Climate_runs: Needs renamed since this currently uses custom cliamte data but runs the OGGM forward model using run_random_climate
 4. data_extraction_scripts: scripts for extracting data from the OGGM run folders 
 5. mb_calibration_scripts: scripts for running and testing the mass balance calibration
 6. shape_files: the various shape files used
 7. figure_scripts: scripts for producing shape files
 8. full_workflow: a folder containing full workflows (largely superseded and can probably be deleted) 
 
 prepro_directory_tester.py a simple script to allow compaison between different prepro directory levels
 

