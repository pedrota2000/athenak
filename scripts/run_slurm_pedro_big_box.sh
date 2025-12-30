#!/bin/bash
#SBATCH --job-name=turb_1024_1gpu
#SBATCH --partition=gpuxl       
#SBATCH --nodes=1
#SBATCH --ntasks=4                         
#SBATCH --gpus-per-task=1                       
#SBATCH --cpus-per-task=24                 
#SBATCH --time=3-00:00:00                    
#SBATCH --mem=0                            
#SBATCH --output=turb_1024_%j.out
#SBATCH --error=turb_1024_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pedro.tarancon@fqa.ub.edu

# Load modules
module purge
module load modules/2.3-20240529 cuda/12.3

# CephFS optimizations
export LD_PRELOAD=/mnt/sw/fi/cephtweaks/lib/libcephtweaks.so
export CEPHTWEAKS_LAZYIO=1

# Set paths
athenak=/mnt/home/ptarancon/athenak
build=$athenak/build

# Create run directory
run_dir=/mnt/ceph/users/ptarancon/runs/turb_1024_run_$(date +%Y%m%d_%H%M%S)
mkdir -p $run_dir
cd $run_dir

# Copy input file
cp $athenak/inputs/custom_tests/3d_turb.athinput ./input.athinput

# Configuration for cleanup
KEEP_LAST_N_BIN=10     # Keep last 10 bin files
KEEP_LAST_N_RST=1      # Keep only last restart file
CHECK_INTERVAL=60      # Check every minute

# Cleanup function
cleanup_old_snapshots() {
    while true; do
        sleep ${CHECK_INTERVAL}
        
        # Cleanup bin files
        file_count=$(ls -1 ${run_dir}/bin/Turb.*.bin 2>/dev/null | wc -l)
        if [ $file_count -gt $KEEP_LAST_N_BIN ]; then
            deleted_count=$((file_count - KEEP_LAST_N_BIN))
            echo "[$(date)] Found ${file_count} bin files, deleting ${deleted_count} oldest..."
            ls -1v ${run_dir}/bin/Turb.*.bin 2>/dev/null | head -n ${deleted_count} | xargs rm -f
            echo "[$(date)] Kept ${KEEP_LAST_N_BIN} most recent bin files."
        fi
        
        # Cleanup restart files - KEEP ONLY LAST ONE
        rst_count=$(ls -1 ${run_dir}/rst/Turb.*.rst 2>/dev/null | wc -l)
        if [ $rst_count -gt $KEEP_LAST_N_RST ]; then
            deleted_count=$((rst_count - KEEP_LAST_N_RST))
            echo "[$(date)] Found ${rst_count} restart files, deleting ${deleted_count} oldest..."
            ls -1v ${run_dir}/rst/Turb.*.rst 2>/dev/null | head -n ${deleted_count} | xargs rm -f
            echo "[$(date)] Kept only the most recent restart file."
        fi
    done
}

# Start cleanup process in background
cleanup_old_snapshots &
CLEANUP_PID=$!

# Trap to ensure cleanup process is killed on exit
trap "kill $CLEANUP_PID 2>/dev/null; exit" SIGINT SIGTERM EXIT

echo "Starting 1024³ turbulence simulation at $(date)"
echo "Running in directory: $run_dir"
echo "Running on 4 tasks with 1 GPU each"
echo "Automatic cleanup enabled:"
echo "  - Keeping last ${KEEP_LAST_N_BIN} bin snapshots"
echo "  - Keeping last ${KEEP_LAST_N_RST} restart file"

# Run simulation
srun -n 4 $build/src/athena -i input.athinput

# Kill cleanup process when simulation finishes
kill $CLEANUP_PID 2>/dev/null

echo "Simulation finished at $(date)"
echo "Final cleanup..."

# Final cleanup for bin files
file_count=$(ls -1 ${run_dir}/bin/Turb.*.bin 2>/dev/null | wc -l)
if [ $file_count -gt $KEEP_LAST_N_BIN ]; then
    deleted_count=$((file_count - KEEP_LAST_N_BIN))
    ls -1v ${run_dir}/bin/Turb.*.bin | head -n ${deleted_count} | xargs rm -f
fi

# Final cleanup for restart files
rst_count=$(ls -1 ${run_dir}/rst/Turb.*.rst 2>/dev/null | wc -l)
if [ $rst_count -gt $KEEP_LAST_N_RST ]; then
    deleted_count=$((rst_count - KEEP_LAST_N_RST))
    ls -1v ${run_dir}/rst/Turb.*.rst | head -n ${deleted_count} | xargs rm -f
fi

echo "Kept ${KEEP_LAST_N_BIN} bin files and ${KEEP_LAST_N_RST} restart file."