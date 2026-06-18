#!/bin/bash
#SBATCH --job-name=turb_128_4gpu
#SBATCH --partition=gpuxl -C h100 --reservation=rocky9
#SBATCH --nodes=1               # 1 nodes = 4 GPUs
#SBATCH --ntasks-per-node=4     # 4 tasks per node
#SBATCH --cpus-per-task=16      
#SBATCH --gres=gpu:4            # 4 GPUs per node
#SBATCH --gpus-per-task=1       # 1 GPU per task
#SBATCH --time=3-00:00:00       
#SBATCH --mem=0                 
#SBATCH --output=turb_128_%j.out
#SBATCH --error=turb_128_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pedro.tarancon@fqa.ub.edu

# Load modules
module purge
module load modules/2.3-20240529 cuda/12.3

# CephFS optimizations
export LD_PRELOAD=/mnt/sw/fi/cephtweaks/lib/libcephtweaks.so
export CEPHTWEAKS_LAZYIO=1

# Set pathsS
athenak=/mnt/home/ptarancon/athenak
build=$athenak/build

# IMPORTANT: Set the directory with existing restart file
# For new run: leave blank or set to empty
# For restart: set to the previous run directory
RESTART_DIR=/mnt/ceph/users/ptarancon/runs/turb_128_run_20260118_085557

# Create run directory
run_dir=/mnt/ceph/users/ptarancon/runs/turb_128_run_$(date +%Y%m%d_%H%M%S)
mkdir -p $run_dir
cd $run_dir

# Check if we're restarting from a previous run
if [ -n "$RESTART_DIR" ] && [ -d "$RESTART_DIR/rst" ]; then
    # Find the most recent restart file
    RESTART_FILE=$(ls -1v ${RESTART_DIR}/rst/Turb.*.rst 2>/dev/null | tail -n 1)
    
    if [ -n "$RESTART_FILE" ]; then
        echo "Found restart file: $RESTART_FILE"
        # Copy restart file to new run directory
        mkdir -p $run_dir/rst
        cp $RESTART_FILE $run_dir/rst/
        RESTART_FILENAME=$(basename $RESTART_FILE)
        RESTART_MODE=true
    else
        echo "No restart file found in $RESTART_DIR/rst, starting from scratch"
        RESTART_MODE=false
    fi
else
    echo "No restart directory specified, starting new simulation"
    RESTART_MODE=false
fi

# Copy input file
cp $athenak/inputs/custom_tests/3d_turb_small.athinput ./input.athinput
sed -i 's/tlim\s*=.*/tlim = 2000.0/' ./input.athinput  # or whatever new limit you want

# Configuration for cleanup
KEEP_LAST_N_BIN=1000      # Keep last 10 bin files
KEEP_LAST_N_RST=1      # Keep only last restart file
CHECK_INTERVAL=60      # Check every minute

# Cleanup function
cleanup_old_snapshots() {
    while true; do
        sleep ${CHECK_INTERVAL}
        
        # Cleanup bin files
        file_count=$(ls -1 ${run_dir}/vtk/Turb.*.vtk 2>/dev/null | wc -l)
        if [ $file_count -gt $KEEP_LAST_N_BIN ]; then
            deleted_count=$((file_count - KEEP_LAST_N_BIN))
            echo "[$(date)] Found ${file_count} vtk files, deleting ${deleted_count} oldest..."
            ls -1v ${run_dir}/vtk/Turb.*.vtk 2>/dev/null | head -n ${deleted_count} | xargs rm -f
            echo "[$(date)] Kept ${KEEP_LAST_N_BIN} most recent vtk files."
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

echo "Starting 128³ turbulence simulation at $(date)"
echo "Running in directory: $run_dir"
echo "Running on 4 GPUs (1 nodes × 4 GPUs)"
echo "Automatic cleanup enabled:"
echo "  - Keeping last ${KEEP_LAST_N_BIN} bin snapshots"
echo "  - Keeping last ${KEEP_LAST_N_RST} restart file"

# Run simulation (with or without restart)
if [ "$RESTART_MODE" = true ]; then
    echo "RESTARTING from: $RESTART_FILENAME"
    srun -n 4 $build/src/athena -r rst/$RESTART_FILENAME -i input.athinput
else
    echo "STARTING NEW simulation"
    srun -n 4 $build/src/athena -i input.athinput
fi

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