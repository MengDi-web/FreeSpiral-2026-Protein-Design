#!/bin/bash
# ============================================================
# 部署到服务器 - 安全版
# 只在你的家目录 (~/FoldSynth/) 下操作，不影响其他人
# 从你的 MacBook 终端运行
# ============================================================

set -e

SERVER="mengd@172.20.73.14"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
LOCAL_DIR="$(dirname "$SCRIPT_DIR")"

echo "============================================"
echo " FoldSynth - 部署到服务器"
echo " 目标: $SERVER:~/FoldSynth/"
echo " 注意: 只使用你的家目录，不影响其他人"
echo "============================================"
echo ""

# 1) 先检查连通性
echo ">>> 检查服务器连通性..."
ssh -o ConnectTimeout=10 -o StrictHostKeyChecking=no "$SERVER" "echo 'OK: $(hostname)'" 2>&1 || {
    echo "❌ 无法连接服务器"
    exit 1
}

# 2) 在服务器上创建目录（必须提前创建，否则 scp 中文路径会失败）
echo ""
echo ">>> 创建目录..."
ssh "$SERVER" "mkdir -p ~/FoldSynth/scripts ~/FoldSynth/results ~/FoldSynth/data"
echo "  ✅ 目录已创建"

# 3) 复制核心文件（逐个复制，避免中文路径问题）
echo ""
echo ">>> 复制核心文件..."

echo "  submission.csv..."
scp "$LOCAL_DIR/submission.csv" "$SERVER:~/FoldSynth/"

echo "  designed_sequences.txt..."
scp "$LOCAL_DIR/designed_sequences.txt" "$SERVER:~/FoldSynth/"

echo "  README.md..."
scp "$LOCAL_DIR/README.md" "$SERVER:~/FoldSynth/"

echo "  design report..."
scp "$LOCAL_DIR/docs/design_report.pdf" "$SERVER:~/FoldSynth/docs/" 2>/dev/null || true

echo "  scripts..."
scp "$LOCAL_DIR/scripts/run_esm3_generation.sh" "$SERVER:~/FoldSynth/scripts/"
scp "$LOCAL_DIR/scripts/run_alphafold2.sh" "$SERVER:~/FoldSynth/scripts/"
scp "$LOCAL_DIR/scripts/run_md_simulations.sh" "$SERVER:~/FoldSynth/scripts/"
scp "$LOCAL_DIR/scripts/run_openmm_md.py" "$SERVER:~/FoldSynth/scripts/"
scp "$LOCAL_DIR/scripts/validate.py" "$SERVER:~/FoldSynth/scripts/"
scp "$LOCAL_DIR/scripts/design_sequences.py" "$SERVER:~/FoldSynth/scripts/"

# 4) 检查服务器环境
echo ""
echo ">>> 服务器环境检查..."
ssh "$SERVER" "
echo '--- GPU ---'
nvidia-smi 2>/dev/null | head -5 || echo '无GPU或nvidia-smi不可用'
echo '--- Python ---'
python3 --version 2>/dev/null || python --version 2>/dev/null
echo '--- 已安装工具 ---'
for cmd in alphafold colabfold_batch gmx; do
    which \$cmd 2>/dev/null && echo \"  \$cmd: 已安装\" || echo \"  \$cmd: 未安装\"
done
echo '--- Python包 ---'
python3 -c 'import openmm; print(\"  openmm: 已安装\")' 2>/dev/null || echo '  openmm: 未安装'
python3 -c 'import torch; print(f\"  torch: 已安装, CUDA={torch.cuda.is_available()}\")' 2>/dev/null || echo '  torch: 未安装'
"

# 5) 显示可用命令
echo ""
echo "============================================"
echo " ✅ 部署完成！"
echo "============================================"
echo ""
echo "现在执行:"
echo ""
echo "  ssh $SERVER"
echo "  cd ~/FoldSynth"
echo ""
echo "查看可用的:"
echo "  ls -la"
echo ""
echo "运行管线:"
echo "  1) ESM3 生成:     bash scripts/run_esm3_generation.sh"
echo "  2) AlphaFold2:     bash scripts/run_alphafold2.sh"
echo "  3) MD 模拟:        python3 scripts/run_openmm_md.py results results/md"
echo ""
echo "传回结果:"
echo "  scp -r $SERVER:~/FoldSynth/results/ ./results_from_server/"
