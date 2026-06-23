# 服务器执行指南

## ⚠️ 安全须知
- 所有操作只在你自己的家目录 (`/home/mengd/`) 下进行
- 不修改 `/etc/`, `/opt/`, `/usr/` 等系统目录
- 不读取或修改其他用户的文件
- 用完的资源（临时文件、大文件）及时清理

---

## 第一步：从 MacBook 连上服务器

打开 MacBook 终端，执行：

```bash
ssh mengd@172.20.73.14
密码: 1218@mandy
```

如果连接不上，检查 MacBook 是否连接了校园网/实验室网络。

---

## 第二步：把项目文件传到服务器

在 MacBook 的**另一个终端窗口**执行（不要关掉 SSH 连接）：

```bash
cd /Users/mandy/Documents/Codex_project/蛋白质设计大赛

# 创建目标目录（你自己的家目录）
ssh mengd@172.20.73.14 "mkdir -p ~/FoldSynth"

# 复制文件（只复制需要的内容）
scp submission.csv designed_sequences.txt mengd@172.20.73.14:~/FoldSynth/
scp -r scripts/ mengd@172.20.73.14:~/FoldSynth/
```

> **为什么要手动 scp 而不是整个目录？**
> - 避免把 `.git/`、`.DS_Store`、`__pycache__/` 等无关文件传到服务器
> - 只传需要计算的文件，节省带宽和服务器空间

或者用我们准备好的部署脚本：

```bash
bash scripts/deploy_to_server.sh
```

---

## 第三步：在服务器上安装依赖

```bash
# SSH 到服务器
ssh mengd@172.20.73.14
cd ~/FoldSynth

# 先看看服务器上有啥工具
echo "=== 已有工具 ==="
which python3
which nvidia-smi
ls /opt/ 2>/dev/null
ls /usr/local/ 2>/dev/null

# 安装需要的包（只安装在你的用户下）
pip3 install --user openmm pdbfixer numpy matplotlib
pip3 install --user esm-sdk  # 如果需要 ESM3 API
```

如果服务器有 conda，用 conda 环境更好：

```bash
# 查看已有环境
conda env list

# 创建自己的环境（不会动别人的）
conda create -n foldsynth python=3.10 -y
conda activate foldsynth
pip install openmm pdbfixer numpy matplotlib
# 如果需要 AlphaFold2
conda install -c bioconda colabfold -y
```

---

## 第四步：检查 GPU

```bash
nvidia-smi
```

看看有没有 GPU、显存多大、有没有人在用。如果有人正在跑大任务，可以先等一等，或者只跑轻量计算（比如只跑 ESM3 API 版，不跑 AlphaFold2 本地预测）。

---

## 第五步：运行计算管线

### 选项 A：ESM3 生成更多候选序列

```bash
cd ~/FoldSynth
mkdir -p results/esm3_generated
export OUTPUT_DIR=~/FoldSynth/results/esm3_generated
bash scripts/run_esm3_generation.sh
```

脚本会尝试用 ESM3 生成序列。如果 API key 不可用，会用模拟模式（基于已有亮度数据产生变异）。

> ⚠️ ESM3 需要 API key：https://forge.evolutionaryscale.ai/

### 选项 B：AlphaFold2 结构预测（需要 GPU）

```bash
# 方式1：用 ColabFold
conda activate foldsynth
colabfold_batch --num-recycle 3 designed_sequences.fasta results/alphafold2/

# 方式2：如果有完整 AlphaFold2
bash scripts/run_alphafold2.sh
```

### 选项 C：MD 模拟 - 72°C 热稳定性测试

轻量版（用 OpenMM，不需要 GROMACS）：

```bash
cd ~/FoldSynth
python3 scripts/run_openmm_md.py \
  results/alphafold2/foldsynth_designs \
  results/md_simulations
```

---

## 第六步：结果传回 MacBook

```bash
# 在 MacBook 终端执行
scp -r mengd@172.20.73.14:~/FoldSynth/results/ \
  ~/Documents/Codex_project/蛋白质设计大赛/results_from_server/
```

---

## 第七步：清理服务器（可选）

如果不再需要，可以删除服务器上的项目文件：

```bash
# 在服务器上执行
rm -rf ~/FoldSynth
# 或者如果用了 conda 环境
conda env remove -n foldsynth
```

---

## 结果解读指南

### AlphaFold2 pLDDT 分数
| 分数 | 含义 | 判断 |
|------|------|------|
| > 90 | 极高置信度 | ✅ 结构可靠 |
| 80-90 | 良好置信度 | ✅ 可以信赖 |
| 70-80 | 中等置信度 | ⚠️ 部分区域可能不准 |
| < 70 | 低置信度 | ❌ 结构不可靠 |

**重点关注发色团区域（TYG，残基65-67附近）的 pLDDT**，这个区域对荧光至关重要。

### MD 模拟结果
| RMSD 漂移 | 含义 | 判断 |
|-----------|------|------|
| < 0.3 nm | 72°C 下结构稳定 | ✅ 热稳定性好 |
| 0.3-0.5 nm | 有波动但可接受 | ⚠️ 中等稳定 |
| > 0.5 nm | 明显去折叠 | ❌ 热稳定性差 |

**发色团水化**（水分子进入发色团周围）是荧光衰减的前兆，需要特别关注。

## 最终提交

根据服务器计算结果：
1. 选择综合表现最好的序列（pLDDT 高 + MD 稳定）
2. 更新 submission.csv
3. 提交到竞赛系统
