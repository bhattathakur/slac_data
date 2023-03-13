#!/bin/bash
dir_name=PZ1M
mkdir ${dir_name}

scratch_path=/scratch/t/tbhatta
destination=${scratch_path}/${dir_name}/merged/

scp tbhatta@sdf-login.slac.stanford.edu:${destination}*.root ./${dir_name}/


