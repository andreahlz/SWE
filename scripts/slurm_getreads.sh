#!/bin/bash
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="${script_dir%/scripts}"
simulation_script="${project_root}/scripts/getfastp_all.sh"

Job_name=iss_sim
cpus=8
memory=20GB
time=02:00:00


mkdir -p "${project_root}/slurm_scripts"
job_file="${project_root}/slurm_scripts/slurm_${Job_name}.job"

cat > "${job_file}" <<EOF
#!/bin/bash
#
#SBATCH --job-name=${Job_name}
#SBATCH --cpus-per-task=${cpus}
#SBATCH --mem=${memory}
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT_50,TIME_LIMIT_80,TIME_LIMIT
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --time=${time}

#conda enviroment
module load Conda/Miniforge3
conda activate readsdnaberts

#into projct root
cd ${project_root}

#start simulation
bash ${simulation_script}
EOF

echo "Slurm-Jobscript created: ${job_file}"
