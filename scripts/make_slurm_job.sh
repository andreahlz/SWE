#!/bin/bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="${script_dir%/scripts}"
pipeline_script="${project_root}/scripts/final_pipeline.nf"

Job_name=SWE_nextflow_pipeline
cpus=4
memory=32GB
time=2-00:00:00

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

#conda enviroment
module load Conda/Miniforge3
conda activate env


export NXF_HOME="${project_root}/.nextflow"
export NXF_WORK="${project_root}/nextflow_work"

#into projct root
cd ${project_root}

#start simulation
nextflow run ${pipeline_script} -profile lisc -resume
EOF

echo "Slurm-Jobscript created: ${job_file}"
