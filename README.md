# SWARM OGGM Model
 A repository containing scripts for running OGGM for the SWARM project
 Currently contains:
 1. CORDEX_runs: workflows for the runs using CORDEX data
 2. MSWEP_ERA_runs: workflows for the runs using MSWEP+ERA data
 3. Custom_Climate_runs: Needs renamed since this currently uses custom cliamte data but runs the OGGM forward model using run_random_climate
 4. data_extraction_scripts: scripts for extracting data from the OGGM run folders 
 5. mb_calibration_scripts: scripts for running and testing the mass balance calibration
 6. shape_files: the various shape files used
 
 prepro_directory_tester.py a simple script to allow compaison between different prepro directory levels
 
 The following scripts have been superseded and can be removed with time
 A script for selecting and outputting the RGI glacier information selected by a shapefile 
 A script for selecting ice inversion by shapefile
