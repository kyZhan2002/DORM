#!/bin/bash

# Get the absolute path of the current script directory
job_directory="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
project_root="$( cd "${job_directory}/../.." && pwd )"

# Create necessary directories
mkdir -p "${job_directory}/job"
mkdir -p "${job_directory}/out"
mkdir -p "${job_directory}/err"

# Check if renv is set up (optional check)
if [ ! -d "${project_root}/renv/library" ]; then
    echo "WARNING: renv library not found. You may need to run setup first."
    echo "Run: sbatch ${job_directory}/setup_job.sh"
    read -p "Continue anyway? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

echo "Submitting simulation jobs..."

for i in {1..3}; do
    random_seed=${i}

    job_file="${job_directory}/job/CondA_L5_${i}.job"

    echo "#!/bin/bash
#SBATCH -c 1
#SBATCH -t 0-10:00
#SBATCH -p short
#SBATCH --output=${job_directory}/out/CondA_L5_${i}.out
#SBATCH --error=${job_directory}/err/CondA_L5_${i}.err
#SBATCH --mem=2G
#SBATCH --job-name=CondA_L5_${i}

# Load required modules
module load gcc/14.2.0 R/4.4.2

# Change to project root directory
cd ${project_root}

# Run simulation script
Rscript ${project_root}/simu/NewSimu_CondA_L5.R ${random_seed}" > $job_file

    sbatch $job_file
    
    # Optional: add a small delay to avoid overwhelming the scheduler
    if [ $((i % 10)) -eq 0 ]; then
        echo "Submitted $i jobs..."
        sleep 1
    fi
done

echo "All $i jobs submitted successfully!"