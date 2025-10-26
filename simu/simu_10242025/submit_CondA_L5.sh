#!/bin/bash

# Get the absolute path of the current script directory
job_directory="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Create necessary directories
mkdir -p job
mkdir -p out
mkdir -p err

for i in {1..500}; do
    random_seed=$RANDOM

    job_file="${job_directory}/job/CondA_L5_${i}.job"

    echo "#!/bin/bash
#SBATCH -c 1
#SBATCH -t 0-12:00
#SBATCH -p short
#SBATCH --output=out/CondA_L5_${i}.out
#SBATCH --error=err/CondA_L5_${i}.err
#SBATCH -x compute-f-17-[09-25]
#SBATCH --mem=4G
#SBATCH --job-name=CondA_L5_${i}

module load rstudio_launcher/1.0 gcc/9.2.0 R/4.4.0
cd ${job_directory}
Rscript ${job_directory}/NewSimu_CondA_L5.R ${random_seed}" > $job_file

    sbatch $job_file
    
    # Optional: add a small delay to avoid overwhelming the scheduler
    if [ $((i % 50)) -eq 0 ]; then
        echo "Submitted $i jobs..."
        sleep 1
    fi
done

echo "All 500 jobs submitted successfully!"