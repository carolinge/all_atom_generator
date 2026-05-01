#!/bin/bash
# 在服务器登录节点上执行此脚本，一键建立 allatom 环境
# 用法：bash setup_server.sh

set -e

echo "=== Step 1: 加载 CUDA 12.5 ==="
module load cuda/12.5

echo "=== Step 2: 建立 conda 环境 ==="
conda env create -f environment_server.yml

echo "=== Step 3: 验证安装 ==="
conda run -n allatom python -c "
import openmm
print('OpenMM version:', openmm.__version__)
from openmm import Platform
for i in range(Platform.getNumPlatforms()):
    p = Platform.getPlatform(i)
    print('  Platform:', p.getName())
"

echo ""
echo "=== 完成！==="
echo "运行模拟："
echo "  module load cuda/12.5"
echo "  conda activate allatom"
echo "  python run_openmm.py"
