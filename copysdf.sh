#!/bin/bash
dir_name=Th228PX
mkdir ${dir_name}

scratch_path=/scratch/t/tbhatta
destination=${scratch_path}/${dir_name}/g4/

scp tbhatta@sdf-login.slac.stanford.edu:${destination}*.root ./${dir_name}/


