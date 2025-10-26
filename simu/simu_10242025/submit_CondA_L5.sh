#!/bin/bash

# Get the absolute path of the current script directory
job_directory="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Create necessary directories
mkdir -p job
mkdir -p out
mkdir -p err

for i in {1..3}; do
    random_seed=${i}

    job_file="${job_directory}/job/CondA_L5_${i}.job"

    echo "#!/bin/bash
#SBATCH -c 1
#SBATCH -t 0-10:00
#SBATCH -p short
#SBATCH --output=out/CondA_L5_${i}.out
#SBATCH --error=err/CondA_L5_${i}.err
#SBATCH --mem=2G
#SBATCH --job-name=CondA_L5_${i}

module load gcc/14.2.0 R/4.4.2
cd ${job_directory}
Rscript ${job_directory}/../NewSimu_CondA_L5.R ${random_seed}" > $job_file

    sbatch $job_file
    
    # Optional: add a small delay to avoid overwhelming the scheduler
    if [ $((i % 1)) -eq 0 ]; then
        echo "Submitted $i jobs..."
        sleep 1
    fi
done

echo "All ${i} jobs submitted successfully!"