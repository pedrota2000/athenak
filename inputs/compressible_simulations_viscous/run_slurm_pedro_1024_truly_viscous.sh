#!/bin/bash
#SBATCH --job-name=turb_1024_truly_visc_4gpu
#SBATCH -p gpuxl
#SBATCH -C h100
## NOTE: the rocky9 reservation ENDED 2026-08-31T00:00:00 (scontrol shows State=INACTIVE),
## so requesting it now would make the job unschedulable. Re-enable only if it is renewed.
##SBATCH --reservation=rocky9
#SBATCH --nodes=1               # 1 node = 4 GPUs
#SBATCH --ntasks-per-node=4     # 4 tasks per node
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:4            # 4 GPUs per node
#SBATCH --gpus-per-task=1       # 1 GPU per task
#SBATCH --time=3-00:00:00
#SBATCH --mem=0
#SBATCH --output=turb_1024_truly_viscous_%j.out
#SBATCH --error=turb_1024_truly_viscous_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pedro.tarancon@fqa.ub.edu

# Load modules
module purge
module load modules/2.4-20250724 cuda/12.5.1   # 2.3 tree + cuda/12.3 were retired; 12.3-built binary runs on the 12.5 runtime

# CephFS optimizations
export LD_PRELOAD=/mnt/sw/fi/cephtweaks/lib/libcephtweaks.so
export CEPHTWEAKS_LAZYIO=1

# --- single-node comms: sidestep InfiniBand memory registration ---
# ulimit -l is 8192 kB on these nodes (not unlimited), so UCX cannot register IB
# memory regions: it spins in a retry loop emitting GB of ib_md/mpool errors and
# never reaches the task list.  Every one of these jobs is --nodes=1, so shared
# memory + CUDA IPC is sufficient and avoids IB entirely.
ulimit -l unlimited 2>/dev/null || true
export UCX_TLS=self,sm,cuda_copy,cuda_ipc
export OMPI_MCA_btl=^openib

# Set paths
athenak=/mnt/home/ptarancon/athenak
build=$athenak/build
input=$athenak/inputs/compressible_simulations_viscous/3d_turb_1024_truly_viscous.athinput

# ---------------------------------------------------------------------------
# Fresh start, or resume.  Export RESUME_DIR=<existing run dir> before sbatch
# to continue that run from its newest restart file instead of starting over.
# ---------------------------------------------------------------------------
if [ -n "${RESUME_DIR:-}" ]; then
    run_dir="$RESUME_DIR"
    rst_file=$(ls -1v ${run_dir}/rst/Turb.*.rst 2>/dev/null | tail -n 1)
    if [ -z "$rst_file" ]; then
        echo "ERROR: no restart file found in ${run_dir}/rst" >&2
        exit 1
    fi
    cd "$run_dir" || exit 1
    athena_args="-r $rst_file"
    echo "RESUMING from $rst_file"
else
    run_dir=/mnt/ceph/users/ptarancon/runs/turb_1024_viscous_compressible_run_truly_viscous_$(date +%Y%m%d_%H%M%S)
    mkdir -p "$run_dir" || exit 1
    cd "$run_dir" || exit 1
    cp "$input" ./input.athinput
    athena_args="-i input.athinput"
fi

# Configuration for cleanup
KEEP_LAST_N_BIN=10   # same snapshot count as the matching inviscid run
KEEP_LAST_N_RST=1            # keep only last restart file
CHECK_INTERVAL=60            # check every minute

# Cleanup function
cleanup_old_snapshots() {
    while true; do
        sleep ${CHECK_INTERVAL}

        # Cleanup bin files
        file_count=$(ls -1 ${run_dir}/bin/Turb.*.bin 2>/dev/null | wc -l)
        if [ $file_count -gt $KEEP_LAST_N_BIN ]; then
            deleted_count=$((file_count - KEEP_LAST_N_BIN))
            echo "[$(date)] Found ${file_count} bin files, deleting ${deleted_count} oldest..."
            ls -1v ${run_dir}/bin/Turb.*.bin 2>/dev/null | head -n ${deleted_count} | xargs -r rm -f
            echo "[$(date)] Kept ${KEEP_LAST_N_BIN} most recent bin files."
        fi

        # Cleanup restart files - KEEP ONLY LAST ONE
        rst_count=$(ls -1 ${run_dir}/rst/Turb.*.rst 2>/dev/null | wc -l)
        if [ $rst_count -gt $KEEP_LAST_N_RST ]; then
            deleted_count=$((rst_count - KEEP_LAST_N_RST))
            echo "[$(date)] Found ${rst_count} restart files, deleting ${deleted_count} oldest..."
            ls -1v ${run_dir}/rst/Turb.*.rst 2>/dev/null | head -n ${deleted_count} | xargs -r rm -f
            echo "[$(date)] Kept only the most recent restart file."
        fi
    done
}

# Start cleanup process in background
cleanup_old_snapshots &
CLEANUP_PID=$!

# Trap to ensure cleanup process is killed on exit
trap "kill $CLEANUP_PID 2>/dev/null; exit" SIGINT SIGTERM EXIT

echo "Starting 1024^3 turbulence simulation (RESOLVED PHYSICAL VISCOSITY) at $(date)"
echo "Running in directory: $run_dir"
echo "Running on 4 GPUs (1 node x 4 GPUs)"
echo "nu from <hydro>/viscosity in: $input"
grep -E "^viscosity" "$input"
echo "Automatic cleanup enabled:"
echo "  - Keeping last ${KEEP_LAST_N_BIN} bin snapshots"
echo "  - Keeping last ${KEEP_LAST_N_RST} restart file"

# Run with 4 tasks (1 node x 4 tasks/node)
srun -n 4 $build/src/athena $athena_args

# Kill cleanup process when simulation finishes
kill $CLEANUP_PID 2>/dev/null

echo "Simulation finished at $(date)"
echo "Final cleanup..."

# Final cleanup for bin files
file_count=$(ls -1 ${run_dir}/bin/Turb.*.bin 2>/dev/null | wc -l)
if [ $file_count -gt $KEEP_LAST_N_BIN ]; then
    deleted_count=$((file_count - KEEP_LAST_N_BIN))
    ls -1v ${run_dir}/bin/Turb.*.bin | head -n ${deleted_count} | xargs -r rm -f
fi

# Final cleanup for restart files
rst_count=$(ls -1 ${run_dir}/rst/Turb.*.rst 2>/dev/null | wc -l)
if [ $rst_count -gt $KEEP_LAST_N_RST ]; then
    deleted_count=$((rst_count - KEEP_LAST_N_RST))
    ls -1v ${run_dir}/rst/Turb.*.rst | head -n ${deleted_count} | xargs -r rm -f
fi

echo "Kept ${KEEP_LAST_N_BIN} bin files and ${KEEP_LAST_N_RST} restart file."
