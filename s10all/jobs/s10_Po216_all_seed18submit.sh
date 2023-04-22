#!/bin/bash
source /etc/profile.d/z10_spack_environment.sh
export NEXOTOP=/sdf/home/t/tbhatta
source ${HOME}/nexo-offline-build/setup.sh 
cd ${HOME}/nexo-offline/Cards
mkdir -p /scratch/t/tbhatta/s10/recon

python ./RunDetSim_new.py --run /scratch/t/tbhatta/s10/g4/s10_Po216_all_seed18.mac --evtmax 25000 --seed 18 --digioutput /scratch/t/tbhatta/s10/g4/s10_Po216_all_seed18.nEXOevents.root --padsize 6 --efield 380 --type bb0n --tilemap ../data/tilesMap_6mm.txt --localmap ../data/localChannelsMap_6mm.txt --wpFile ../data/singleStripWP6mm.root --noiselib ../data/noise_lib_1_2_us_100e.root --eleclife 10000 --sampling 0.5 --induc True --swf True --oversampleratio 50 --coType Pads --skipEThreshold 0.1 > /scratch/t/tbhatta/s10/g4/s10_Po216_all_seed18.nEXOevents.out 2> /scratch/t/tbhatta/s10/g4/s10_Po216_all_seed18.nEXOevents.err 

python ../Analysis/SensitivityRecon/share/RunSensRecon.py --input /scratch/t/tbhatta/s10/g4/s10_Po216_all_seed18.nEXOevents.root --evtmax 25000 --seed 18 --output /scratch/t/tbhatta/s10/recon/s10_Po216_all_seed18.reconTree.root > /scratch/t/tbhatta/s10/recon/s10_Po216_all_seed18.reconTree.out 2> /scratch/t/tbhatta/s10/recon/s10_Po216_all_seed18.reconTree.err 

