#!/bin/bash
dir_name=PXE01M
mkdir ${dir_name}

scratch_path=/scratch/t/tbhatta
destination=${scratch_path}/${dir_name}/merged/
#destination=${scratch_path}/${dir_name}/recon/

scp tbhatta@sdf-login.slac.stanford.edu:${destination}*.root ./${dir_name}/


