#!/bin/bash
#SBATCH --job-name=turb_512_1gpu
#SBATCH --partition=gpu -C a100-80gb
#SBATCH --nodes=1
#SBATCH --ntasks=1                         # Single task
#SBATCH --gres=gpu:1                       # 1 GPU
#SBATCH --cpus-per-task=32                 # Use more CPUs
#SBATCH --time=96:00:00                    # More time for 512³
#SBATCH --mem=0                            # Use all available memory
#SBATCH --output=turb_512_%j.out
#SBATCH --error=turb_512_%j.err

# Load modules (NO openmpi needed!)
module purge
module load modules/2.3-20240529 cuda/12.3

# CephFS optimizations
export LD_PRELOAD=/mnt/sw/fi/cephtweaks/lib/libcephtweaks.so
export CEPHTWEAKS_LAZYIO=1

# Set paths
athenak=/mnt/home/ptarancon/athenak
build=$athenak/build

# Create run directory
run_dir=/mnt/ceph/users/ptarancon/runs/turb_512_run_$(date +%Y%m%d_%H%M%S)
mkdir -p $run_dir
cd $run_dir

# Copy input file
cp $athenak/inputs/custom_tests/3d_turb.athinput ./input.athinput

echo "Starting 512³ turbulence simulation at $(date)"
echo "Running in directory: $run_dir"
echo "Running on 1 A100-80GB GPU"

# Run directly (NO mpirun!)
$build/src/athena -i input.athinput

echo "Simulation finished at $(date)"