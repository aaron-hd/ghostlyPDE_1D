#!/bin/bash
#SBATCH --job-name=freq-dep		# create a short name for your job
#SBATCH --nodes=1								# node count
#SBATCH --ntasks=1							# total number of tasks (one per Julia script)
#SBATCH --cpus-per-task=4				# cpu-cores per task
#SBATCH --mem-per-cpu=1G			# memory per task (4 tasks x 32GB = 128GB total)
#SBATCH --time=1:00:00					# total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin				# send email when job begins
#SBATCH --mail-type=end					# send email when job ends
#SBATCH --mail-user=aaron.held@phys.ens.fr

# Run four Julia scripts in parallel
cd 01_wave/mchi2\=0/
julia --threads 4 main.jl
cd ../../01_wave/mchi2\=0.01lam/
julia --threads 4 main.jl
cd ../../01_wave/mchi2\=0.1lam/
julia --threads 4 main.jl
cd ../../01_wave/mchi2\=1lam/
julia --threads 4 main.jl
cd ../../01_wave/mchi2\=10lam/
julia --threads 4 main.jl
cd ../../01_wave/mchi2\=100lam/
julia --threads 4 main.jl
cd ../../01_wave/effective/
julia --threads 4 main.jl
