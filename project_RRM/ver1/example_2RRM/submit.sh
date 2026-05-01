#!/bin/bash

#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --partition=graphic          # ← 确认你们集群的 GPU 分区名
#SBATCH --constraint="gpu"
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=10
#SBATCH --mem=225000

#SBATCH --mail-type=END
#SBATCH --job-name="2RRM_openmm"
#SBATCH --output=./%j.out
#SBATCH --error=./%j.err

set -e

# 加载 CUDA
module purge
module load cuda/12.5

# 激活 conda 环境
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate allatom

# 运行模拟
cd "$(dirname "$0")"
python run_openmm.py

exit 0
