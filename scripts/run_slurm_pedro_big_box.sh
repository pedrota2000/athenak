#!/bin/bash
#SBATCH --job-name=turb_1024_1gpu
#SBATCH --partition=gpuxl       # partition name
#SBATCH --nodes=1
#SBATCH --ntasks=4                         # Single task
#SBATCH --gpus-per-task=1                       # 4 GPUs
#SBATCH --cpus-per-task=24                 # Use more CPUs
#SBATCH --time=3-00:00:00                    # More time for 512³
#SBATCH --mem=0                            # Use all available memory
#SBATCH --output=turb_1024_%j.out
#SBATCH --error=turb_1024_%j.err

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
run_dir=/mnt/ceph/users/ptarancon/runs/turb_1024_run_$(date +%Y%m%d_%H%M%S)
mkdir -p $run_dir
cd $run_dir

# Copy input file
cp $athenak/inputs/custom_tests/3d_turb.athinput ./input.athinput

# Configuration for cleanup
KEEP_LAST_N=10
CHECK_INTERVAL=60  # Check every 5 minutes (300 seconds)

# Cleanup function
cleanup_old_snapshots() {
    while true; do
        sleep ${CHECK_INTERVAL}
        
        # Count bin files
        file_count=$(ls -1 ${run_dir}/bin/Turb.*.bin 2>/dev/null | wc -l)
        
        if [ $file_count -gt $KEEP_LAST_N ]; then
            # Delete oldest files by FILENAME (assuming numeric sequence in filename)
            deleted_count=$((file_count - KEEP_LAST_N))
            echo "[$(date)] Found ${file_count} files, deleting ${deleted_count} oldest..."
            
            # Sort by filename naturally (handles numbers correctly)
            ls -1v ${run_dir}/bin/Turb.*.bin 2>/dev/null | head -n ${deleted_count} | xargs rm -f
            
            echo "[$(date)] Cleanup complete. Kept ${KEEP_LAST_N} most recent snapshots."
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
echo "Running on 1 H100 GPU"
echo "Automatic cleanup enabled: keeping last ${KEEP_LAST_N} snapshots"

# Run directly 
srun -n 4 $build/src/athena -i input.athinput

# Kill cleanup process when simulation finishes
kill $CLEANUP_PID 2>/dev/null

echo "Simulation finished at $(date)"
echo "Final cleanup..."
# One final cleanup to ensure we're at exactly KEEP_LAST_N files
file_count=$(ls -1 ${run_dir}/bin/Turb.*.bin 2>/dev/null | wc -l)
if [ $file_count -gt $KEEP_LAST_N ]; then
    deleted_count=$((file_count - KEEP_LAST_N))
    ls -1v ${run_dir}/bin/Turb.*.bin | head -n ${deleted_count} | xargs rm -f
fi

echo "Kept ${KEEP_LAST_N} most recent snapshots."