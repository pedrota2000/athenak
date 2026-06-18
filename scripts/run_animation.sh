#!/bin/bash
#SBATCH --job-name=animation
#SBATCH --partition=gen
#SBATCH -C rome
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem=128G
#SBATCH --output=animation_%j.out
#SBATCH --error=animation_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pedro.tarancon@fqa.ub.edu

module load ffmpeg

cd /mnt/home/ptarancon/athenak/scripts
source /mnt/home/ptarancon/python_env/neurodiffeq/bin/activate

python do_animation_vorticity.py