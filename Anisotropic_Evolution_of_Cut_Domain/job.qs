#!/bin/bash
#SBATCH --job-name=Anisotropic_Evolution_Cut_Domain
#SBATCH --partition=mathsci
#SBATCH --output=my_job_op%j.txt
#SBATCH --mem=100G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=nmirzaei@udel.edu
#SBATCH --mail-type=END


vpkg_require fenics


Sexec python3 Main_script.py
