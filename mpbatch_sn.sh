#!/bin/bash

#SBATCH -J udm_sn

#SBATCH -o out_udm_sn

#SBATCH -e errorout_udm_sn

#SBATCH -n 8

###SBATCH -N 1         -----> reserve one node (24 cores)

###SBATCH --partition=test

#SBATCH --time=3-00:00:00

echo "Running MontePython for the UDM case with Gaussian prior on M and $N cores"

###module load OpenMPI/4.0.1-GCC-8.3.0-2.32
module load puck_openmpi

export OMP_NUM_THREADS=8

source /n/work02/tmiranda/montepython_public/code/plc_3.0/plc-3.01/bin/clik_profile.sh

###mpirun -np 8 -mca btl ^openib python3 /n/work02/tmiranda/montepython_public/montepython/MontePython.py run -o chains/udm_sn -p input/udm_sn.param -c chains/udm_sn_M_run1/udm_sn_M_run1.covmat -b chains/udm_sn_M_run1/udm_sn_M_run1.bestfit --superupdate 20 -N 300000 --silent

mpirun -n 8 -mca btl ^openib  python3 /n/work02/tmiranda/montepython_public/montepython/MontePython.py run -p input/udm_sn.param -o chains/udm_sn -r chains/udm_sn/2023-01-31_600000__1.txt -N 400000 --superupdate 20 --silent

###You can check if the code is running by typing: "squeue -u username".
