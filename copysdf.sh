#!/bin/bash
dir_name=10Mth228
mkdir ${dir_name}

scratch_path=/scratch/t/tbhatta
destination=${scratch_path}/${dir_name}/merged/
#destination=${scratch_path}/${dir_name}/recon/

scp tbhatta@sdf-login.slac.stanford.edu:${destination}*.root ./${dir_name}/


