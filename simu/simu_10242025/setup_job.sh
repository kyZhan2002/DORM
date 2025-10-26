
#!/bin/bash
#SBATCH -c 1
#SBATCH -t 0-01:00
#SBATCH -p short
#SBATCH --mem=4G
#SBATCH --job-name=renv_setup
#SBATCH --output=setup.out
#SBATCH --error=setup.err

# Get the project root directory (two levels up from this script)
job_directory="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
project_root="$( cd "${job_directory}/../.." && pwd )"

echo "=== renv Setup Job ==="
echo "Project root: ${project_root}"
echo "Starting at: $(date)"

# Load required modules
module load gcc/14.2.0 R/4.4.2

# Change to project root
cd ${project_root}

# Run setup script
Rscript simu/simu_10242025/setup_cluster.R

echo "Completed at: $(date)"