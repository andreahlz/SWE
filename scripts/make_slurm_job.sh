#!/bin/bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="${script_dir%/scripts}"
pipeline_script="${project_root}/scripts/parallized_pipeline.nf"

Job_name=SWE_nextflow_pipeline
cpus=4
memory=8GB
time=1-00:00:00

mkdir -p "${project_root}/slurm_scripts" "${project_root}/logs"
mkdir -p "${project_root}/.nextflow" "${project_root}/nextflow_work"


job_file="${project_root}/slurm_scripts/slurm_${Job_name}.job"

cat > "${job_file}" <<EOF
#!/bin/bash
#
#SBATCH --job-name=${Job_name}
#SBATCH --cpus-per-task=${cpus}
#SBATCH --mem=${memory}
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT_50,TIME_LIMIT_80,TIME_LIMIT
#SBATCH --output=${project_root}/logs/%x-%j.out
#SBATCH --error=${project_root}/logs/%x-%j.err
#SBATCH --time=${time}

module load anaconda3

#conda enviroment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate readsdnaberts


export NXF_HOME="${project_root}/.nextflow"
export NXF_WORK="${project_root}/nextflow_work"

#into projct root
cd ${project_root}

#start simulation
nextflow run ${pipeline_script} -profile lisc -resume
EOF

echo "Slurm-Jobscript created: ${job_file}"
