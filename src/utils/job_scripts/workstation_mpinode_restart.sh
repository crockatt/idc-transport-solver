#!/bin/bash -login

#SBATCH --job-name <__JOBNAME__>
#SBATCH --time=1:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --partition=batch

alias date='date +"%F  %H:%M:%S"'

# Prints basic job description (used for emails).
function job_stats {
    echo "[$(date)]  "$(hostname)
    echo "[$(date)]  "$(pwd)
    echo "[$(date)]  "${JOBPREFIX}
}

export JOBPREFIX="<__JOBNAME__>"

if [ ! -z "${SLURM_SUBMIT_DIR}" ]
then
    cd ${SLURM_SUBMIT_DIR}
fi

echo "[$(date)]  "$(hostname)
echo "[$(date)]  "$(pwd)

# Set library and executable path environment variables.
source setenv.sh

# Set OpenMP environment variables.
export OMP_NUM_THREADS=1
export OMP_PROC_BIND=false
export OMP_DYNAMIC=false
export OMP_SCHEDULE=static

# Run the code.
mpiexec \
    --bind-to core \
    --mca btl ^tcp \
    --npernode 39 \
    solver_<__SPACE_DIMS__>d.x ${JOBPREFIX}.deck

# Check return codes.
RET=$?
echo "[$(date)]  RET="${RET}

if [ ! -e ./finished_jobs/ ]; then mkdir -vp ./finished_jobs/ ; fi

if [[ ( "${RET}" -eq "0" ) ]]
then

    # Job completed without error.
    echo "[$(date)]  Run complete."

    rm ${JOBPREFIX}*.chk
    mv ${JOBPREFIX}.* ./finished_jobs/

elif [[ ( "${RET}" -eq "150" ) ]]
then

    # Job checkpointed.
    echo "[$(date)]  Job checkpointed. Linking checkpoint and resubmitting job:"

    rm -v ${JOBPREFIX}.chk
    ln -sv $( ls ${JOBPREFIX}.*.chk | sort | tail -n 1 ) ${JOBPREFIX}.chk

    if [ ! -z "${SLURM_SUBMIT_DIR}" ]
    then
        sbatch ${JOBPREFIX}.sub
    fi
    
else

    # Some error occurred.
    echo "[$(date)]  Run terminated with error ${RET}."
fi

if [ ! -s "${SLURM_JOB_ID}" ]
then
    mv -v slurm-${SLURM_JOB_ID}.out ./finished_jobs/
fi

exit 0
