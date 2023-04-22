#!/bin/bash
source /etc/profile.d/z10_spack_environment.sh
export NEXOTOP=/sdf/home/t/tbhatta
source ${HOME}/nexo-offline-build/setup.sh 
cd ${HOME}/nexo-offline/Cards
python ./MergeIntoPlainTrees.py -c yamls/s10.yaml -r -k all -p Th228  > /scratch/t/tbhatta/s10/merge-jobs/s10_Th228_all.out 2> /scratch/t/tbhatta/s10/merge-jobs/s10_Th228_all.err
