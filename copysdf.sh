#!/bin/bash
dir_name=s2
mkdir ${dir_name}

scratch_path=/scratch/t/tbhatta
destination=${scratch_path}/${dir_name}/merged/
#destination=${scratch_path}/${dir_name}/recon/

sdf=tbhatta@sdf-login.slac.stanford.edu
scp ${sdf}:${scratch_path}/${dir_name}/${dir_name}.yaml ./${dir_name}/
scp ${sdf}:${destination}*.root ./${dir_name}/


