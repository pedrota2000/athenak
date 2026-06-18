#!/bin/bash
#SBATCH --job-name=compress_vtk
#SBATCH --partition=gen
#SBATCH -C rome
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --output=compress_%j.out
#SBATCH --error=compress_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pedro.tarancon@fqa.ub.edu

# Go to your data directory
cd /mnt/ceph/users/ptarancon/runs/turb_128_run_20260223_130640

echo "Starting compression at $(date)"

# Compression command with progress feedback
tar -cf - vtk/ \
  --checkpoint=1000 \
  --checkpoint-action=exec='echo "Processed $TAR_CHECKPOINT files"' \
| zstd -3 -T0 -o data.tar.zst

echo "Finished at $(date)"