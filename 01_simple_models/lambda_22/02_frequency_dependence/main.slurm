#!/bin/bash
#SBATCH --job-name=freq-dep		# create a short name for your job
#SBATCH --nodes=1								# node count
#SBATCH --ntasks=1							# total number of tasks (one per Julia script)
#SBATCH --cpus-per-task=8				# cpu-cores per task
#SBATCH --mem-per-cpu=32G			# memory per task (4 tasks x 32GB = 128GB total)
#SBATCH --time=48:00:00					# total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin				# send email when job begins
#SBATCH --mail-type=end					# send email when job ends
#SBATCH --mail-user=aaron.held@phys.ens.fr

# Run four Julia scripts in parallel
cd lam\=10/
cd k\=1/
julia --threads 8 main.jl
cd ../k\=8/
julia --threads 8 main.jl
cd ../k\=32/
julia --threads 8 main.jl
cd ../..

cd lam\=1/
cd k\=1/
julia --threads 8 main.jl
cd ../k\=8/
julia --threads 8 main.jl
cd ../k\=32/
julia --threads 8 main.jl
cd ../..

cd lam\=10_m\=10/
cd k\=1/
julia --threads 8 main.jl
cd ../k\=8/
julia --threads 8 main.jl
cd ../k\=32/
julia --threads 8 main.jl
cd ../..
