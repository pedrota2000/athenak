#!/bin/bash
#SBATCH --job-name=animation
#SBATCH --partition=gen
#SBATCH -C rome
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem=256G
#SBATCH --output=animation_%j.out
#SBATCH --error=animation_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pedro.tarancon@fqa.ub.edu

module purge

cd /mnt/home/ptarancon/athenak/vis/python
source /mnt/home/ptarancon/python_env/neurodiffeq/bin/activate

ls /mnt/ceph/users/ptarancon/runs/turb_1024_run_viscous_20260316_120419/bin/Turb.hydro_w.02010.bin \
    | sed 's/\.bin$//' \
    | xargs -P 1 -I {} python3 make_athdf.py {}