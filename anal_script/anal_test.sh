#!/bin/bash
# =============================================================================
# Anal Test MD Workflow V6.0 (Ultimate Stabilized Version)
# 修复：TPR原子匹配、hbond/mindist标签、多重LIG组歧义、全路径调用
# =============================================================================

set -euo pipefail

# --- 环境与参数配置 ---
PDB_ROOT=$(realpath "${1:-.}")
LIGNAME=${2:-LIG}
GPU_ID=${3:-0}
# 锁定脚本所在的物理目录，用于精准调用同目录下的 python 绘图脚本
SCRIPT_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)

MMPBSA_TEMPLATE="${SCRIPT_ROOT}/mmpbsa.in"

GMX=${GMX:-gmx}
MPI_BIN=${MPI_BIN:-mpirun}
MMPBSA_BIN=${MMPBSA_BIN:-gmx_MMPBSA}

# 颜色输出函数
log_info() { echo -e "\033[0;34m[INFO]\033[0m $1"; }
log_success() { echo -e "\033[0;32m[SUCCESS]\033[0m $1"; }
log_warning() { echo -e "\033[0;33m[WARNING]\033[0m $1"; }
log_error() { echo -e "\033[0;31m[ERROR]\033[0m $1"; }

# =============================================================================
# 函数：分析单个 Replica
# =============================================================================
analyze_replica() {
    local r_dir="$1"
    local r_abs=$(realpath "$r_dir")
    local prod_tpr=$(find "$r_abs/results/prod" -name "md_*_prod.tpr" -print -quit)
    local clean_full=$(find "$r_abs/results/prod" -name "md_*_clean_full.xtc" -print -quit)

    if [[ -z "$prod_tpr" || -z "$clean_full" ]]; then return; fi

    local analysis_dir="$r_abs/analysis"
    mkdir -p "$analysis_dir"
    pushd "$analysis_dir" >/dev/null

    # 1. 精准合并：使用编号 1 (Protein) 和 13 (LIG) 确保含氢且无歧义
    log_info "正在提取纯净复合物 (Protein + LIG, 排除离子)..."
    "$GMX" make_ndx -f "$prod_tpr" -o raw.ndx <<EOF
1 | 13
name 22 Complex_Pure
q
EOF

    # 2. 生成 match.tpr 和 match.xtc (原子数必须为 4576)
    log_info "生成原子数为 4576 的分析文件..."
    echo "Complex_Pure" | "$GMX" convert-tpr -s "$prod_tpr" -n raw.ndx -o match.tpr > /dev/null 2>&1
    echo "Complex_Pure" | "$GMX" trjconv -s "$prod_tpr" -f "$clean_full" -n raw.ndx -o match.xtc > /dev/null 2>&1

    # 3. 创建最终索引
    "$GMX" make_ndx -f match.tpr -o Index4_anal.ndx <<EOF
r $LIGNAME
q
EOF

    # 4. 执行常规 GMX 分析任务 (统一使用 match.xtc)
    # 编号逻辑：0-System, 1-Protein, 3-C-alpha, 4-Backbone, 13-LIG (基于 4576 原子的 match.tpr)
    
    log_info "计算 RMSD (蛋白 Cα)..."
    echo "3 3" | "$GMX" rms -s match.tpr -f match.xtc -o rmsd_protein_ca.xvg -n Index4_anal.ndx -tu ns <<EOF
3
3
EOF

    log_info "计算 RMSD (配体相对蛋白)..."
    echo "3 13" | "$GMX" rms -s match.tpr -f match.xtc -o rmsd_ligand.xvg -n Index4_anal.ndx -tu ns <<EOF
3
13
EOF

    log_info "计算 RMSF (Backbone)..."
    echo "4" | "$GMX" rmsf -s match.tpr -f match.xtc -o rmsf.xvg -res <<EOF
4
EOF

    log_info "计算回转半径 (Rg)..."
    echo "1" | "$GMX" gyrate -s match.tpr -f match.xtc -o rg_total.xvg <<EOF
1
EOF

    log_info "计算 SASA (System)..."
    echo "0" | "$GMX" sasa -s match.tpr -f match.xtc -o sasa.xvg <<EOF
0
EOF

    log_info "计算氢键..."
    echo "1 13" | "$GMX" hbond -s match.tpr -f match.xtc -n Index4_anal.ndx -num hbnum.xvg -tu ns <<EOF
1
13
EOF

    log_info "计算最小距离..."
    echo "1 13" | "$GMX" mindist -s match.tpr -f match.xtc -n Index4_anal.ndx -od mindist.xvg -on numcont.xvg -d 0.5 -tu ns <<EOF
1
13
EOF

    # 6. 提取 2D 自由能景观数据 (RMSD vs Rg)
    log_info "提取 2D FES 数据..."
    python3 <<EOF
import pandas as pd
import os
def get_v(p):
    d = []
    if not os.path.exists(p): return None
    with open(p) as f:
        for l in f:
            if not l.startswith(('#','@')):
                pts = l.split()
                if len(pts)>=2: d.append([float(pts[0]), float(pts[1])])
    return d
r = get_v('rmsd_protein_ca.xvg')
g = get_v('rg_total.xvg')
if r and g:
    ml = min(len(r), len(g))
    pd.DataFrame({'Time':[x[0] for x in r[:ml]], 'RMSD_nm':[x[1] for x in r[:ml]], 'Rg_nm':[x[1] for x in g[:ml]]}).to_csv('fes_2d_data.csv', index=False)
EOF

    # 6. 执行 MMPBSA
    if [[ -f "$MMPBSA_TEMPLATE" ]]; then
        cp "$MMPBSA_TEMPLATE" "mmpbsa.in"
        log_info "运行 gmx_MMPBSA (受体=1, 配体=13, np=14)..."
        local top_file=$(find "$r_abs/results/prod" -name "topol.top" -print -quit)
        
        "$MPI_BIN" -np 14 "$MMPBSA_BIN" -O \
            -i mmpbsa.in \
            -cs match.tpr \
            -ct match.xtc \
            -ci Index4_anal.ndx \
            -cg 1 13 \
            -cp "$top_file" \
            -o mmpbsa_results.dat \
            -eo DECOMP_RESULTS_MMPBSA.csv

        # 执行子模块：残基 Occupancy 分析
        popd >/dev/null
        process_residue_occupancy "$r_abs" "$analysis_dir/mmgbsa_results.dat"
        pushd "$analysis_dir" >/dev/null
    fi

    popd >/dev/null
}

# =============================================================================
# 函数：合并数据 (修正核心逻辑)
# =============================================================================
# =============================================================================
# 函数：汇总所有数据并执行全指标绘图 (V25.0 增强版)
# =============================================================================
merge_and_plot_all() {
    local pdb_path="$1"
    local merged_dir="$pdb_path/merged_analysis"
    mkdir -p "$merged_dir"
    
    log_info "正在汇总合并数据、计算残基 Occupancy 并启动自动化绘图..."
    pushd "$pdb_path" >/dev/null

    # --- 第一部分：Python 数据聚合处理 ---
    python3 <<EOF
import pandas as pd
import glob
import os
import numpy as np

def read_xvg(path, cols=[0,1], is_rg=False):
    if not os.path.exists(path): return None
    data = []
    with open(path) as f:
        for line in f:
            if line.startswith(('#', '@')): continue
            pts = line.strip().split()
            if pts:
                row = [float(x) for x in pts]
                if is_rg: row[0] = row[0] / 1000.0 # ps -> ns
                data.append([row[i] for i in cols if i < len(row)])
    return pd.DataFrame(data) if data else None

def merge_rep(pattern, out_name, x_label='Time', is_rg=False):
    files = sorted(glob.glob(f'replica_*/analysis/{pattern}'))
    if not files: return
    df0 = read_xvg(files[0], is_rg=is_rg)
    if df0 is None: return
    merged = pd.DataFrame({x_label: df0[0]})
    for f in files:
        rep = f.split('/')[0]
        df = read_xvg(f, is_rg=is_rg)
        if df is not None:
            merged[rep] = df[1].values[:len(merged)]
    os.makedirs('merged_analysis', exist_ok=True)
    merged.to_csv(f'merged_analysis/{out_name}', index=False)

# 1. 执行常规全局指标合并
merge_rep('rmsd_protein_ca.xvg', 'rmsd_protein_ca_merged.csv')
merge_rep('rmsd_ligand.xvg', 'rmsd_ligand_merged.csv')
merge_rep('rmsf.xvg', 'rmsf_merged.csv', x_label='Residue')
merge_rep('hbnum.xvg', 'hbond_merged.csv')
merge_rep('mindist.xvg', 'mindist_merged.csv')
merge_rep('numcont.xvg', 'numcont_merged.csv')
merge_rep('sasa.xvg', 'sasa_merged.csv')
merge_rep('rg_total.xvg', 'rg_total_merged.csv', is_rg=True)

# 2. Rg 四轴分量 (基于副本1)
rg_f = glob.glob('replica_1/analysis/rg_total.xvg')
if rg_f:
    df_rg = read_xvg(rg_f[0], cols=[0,1,2,3,4], is_rg=True)
    if df_rg is not None:
        df_rg.columns = ['Time', 'Rg_Total', 'Rg_X', 'Rg_Y', 'Rg_Z']
        df_rg.to_csv('merged_analysis/rg_merged.csv', index=False)

# 3. 核心：聚合各氨基酸残基的 Occupancy (跨副本计算)
def aggregate_residue_occupancy():
    all_reps_data = []
    # 查找所有副本下的残基接触文件 (由之前的模块生成)
    rep_dirs = sorted(glob.glob('replica_*'))
    for rd in rep_dirs:
        cont_files = glob.glob(f'{rd}/analysis/res_occupancy/cont_*.xvg')
        for f in cont_files:
            res_id = f.split('_')[-1].replace('.xvg','')
            raw_data = read_xvg(f)
            if raw_data is not None:
                # 计算该副本该残基的 Occupancy: (接触帧数 / 总帧数) * 100
                contacts = raw_data[1].values
                occ = (np.count_nonzero(contacts > 0) / len(contacts)) * 100
                all_reps_data.append({
                    'Replica': rd,
                    'Residue': f'Res {res_id}',
                    'Occupancy': occ,
                    'ResNum': int(res_id)
                })
    
    if all_reps_data:
        df_all = pd.DataFrame(all_reps_data)
        # 跨副本计算平均值和标准差
        summary = df_all.groupby(['ResNum', 'Residue'])['Occupancy'].agg(['mean', 'std']).reset_index()
        summary.columns = ['ResNum', 'Residue', 'Occupancy_Mean', 'Occupancy_Std']
        summary = summary.sort_values('ResNum')
        os.makedirs('merged_analysis/res_analysis', exist_ok=True)
        summary.to_csv('merged_analysis/res_analysis/residue_occupancy_summary.csv', index=False)
        print("已生成残基 Occupancy 汇总表。")

aggregate_residue_occupancy()
EOF

    # --- 第二部分：调用绘图脚本 ---
    log_info "启动全量自动化绘图..."
    pushd "merged_analysis" >/dev/null

    # 基础分析图
    [[ -f "rmsd_protein_ca_merged.csv" ]] && python3 "$SCRIPT_ROOT/plot_rmsd_analysis.py" rmsd_protein_ca_merged.csv "Protein Cα RMSD" "RMSD (nm)"
    [[ -f "rmsd_ligand_merged.csv" ]] && python3 "$SCRIPT_ROOT/plot_rmsd_analysis.py" rmsd_ligand_merged.csv "Ligand RMSD" "RMSD (nm)"
    [[ -f "rmsf_merged.csv" ]] && python3 "$SCRIPT_ROOT/plot_rmsf_analysis.py" rmsf_merged.csv
    [[ -f "rg_merged.csv" ]] && python3 "$SCRIPT_ROOT/plot_rg_analysis.py" rg_merged.csv
    [[ -f "hbond_merged.csv" ]] && python3 "$SCRIPT_ROOT/plot_hbond_analysis.py" hbond_merged.csv
    [[ -f "sasa_merged.csv" ]] && python3 "$SCRIPT_ROOT/plot_sasa_analysis.py" sasa_merged.csv
    
    # --- 景观图分析 ---
    [[ -f "fes_2d_merged.csv" ]] && python3 "$SCRIPT_ROOT/plot_fes_2d.py" fes_2d_merged.csv


    # 2. 关键：解析 MMPBSA 分解能（基于你指定的 .dat 文件）
    if [[ -f "$DECOMP_DAT" ]]; then
        log_info "正在绘制残基结合能贡献图..."
        python3 "$SCRIPT_ROOT/plot_mmpbsa_decomp.py" "$DECOMP_DAT"
    else
        log_warning "未找到能量分解文件 $DECOMP_DAT，跳过贡献图绘制。"
    fi

    # 最小距离与接触数全局图
    [[ -f "mindist_merged.csv" ]] && python3 "$SCRIPT_ROOT/plot_mindist_analysis.py" mindist_merged.csv
    [[ -f "numcont_merged.csv" ]] && python3 "$SCRIPT_ROOT/plot_numcont_analysis.py" numcont_merged.csv

    # 关键：执行残基 Occupancy 柱状图绘制
    if [[ -f "res_analysis/residue_occupancy_summary.csv" ]]; then
        log_info "绘制残基级接触频率 (Occupancy) 柱状图..."
        python3 <<EOF
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('res_analysis/residue_occupancy_summary.csv')
plt.figure(figsize=(14, 7), dpi=300)
# 绘制带误差棒的柱状图
plt.bar(df['Residue'], df['Occupancy_Mean'], yerr=df['Occupancy_Std'], 
        color='#4C72B0', edgecolor='black', capsize=5, alpha=0.8)

plt.ylabel('Contact Occupancy (%)', fontsize=14)
plt.xlabel('Residues (Identified by MMPBSA)', fontsize=14)
plt.title('Binding Site Residue Interaction Frequency (Occupancy)', fontsize=16)
plt.xticks(rotation=45, ha='right')
plt.ylim(0, 110)
plt.grid(axis='y', linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig('res_analysis/occupancy_bar_chart.png')
EOF
    fi

    popd >/dev/null
    popd >/dev/null
    log_success "体系 $(basename "$pdb_path") 的所有数据合并与绘图任务已完成。"
}



# =============================================================================
# 函数：自动化绘图集成 (调用同目录下的脚本)
# =============================================================================
plot_all_results() {
    local pdb_path="$1"
    local merged_dir="$pdb_path/merged_analysis"
    log_info "正在调用 Python 脚本进行全指标自动化绘图..."

    pushd "$merged_dir" >/dev/null
    
    # --- 基础结构分析 ---
    [[ -f "rmsd_protein_ca_merged.csv" ]] && python3 "$SCRIPT_ROOT/plot_rmsd_analysis.py" rmsd_protein_ca_merged.csv "Protein Cα RMSD" "RMSD (nm)"
    [[ -f "rmsd_ligand_merged.csv" ]] && python3 "$SCRIPT_ROOT/plot_rmsd_analysis.py" rmsd_ligand_merged.csv "Ligand RMSD" "RMSD (nm)"
    [[ -f "rmsf_merged.csv" ]] && python3 "$SCRIPT_ROOT/plot_rmsf_analysis.py" rmsf_merged.csv
    [[ -f "rg_merged.csv" ]] && python3 "$SCRIPT_ROOT/plot_rg_analysis.py" rg_merged.csv
    [[ -f "sasa_merged.csv" ]] && python3 "$SCRIPT_ROOT/plot_sasa_analysis.py" sasa_merged.csv
    
    # --- 相互作用分析 ---
    [[ -f "hbond_merged.csv" ]] && python3 "$SCRIPT_ROOT/plot_hbond_analysis.py" hbond_merged.csv
    
    # --- 景观图分析 ---
    [[ -f "fes_2d_merged.csv" ]] && python3 "$SCRIPT_ROOT/plot_fes_2d.py" fes_2d_merged.csv


    # 2. 关键：解析 MMPBSA 分解能（基于你指定的 .dat 文件）
    if [[ -f "$DECOMP_DAT" ]]; then
        log_info "正在绘制残基结合能贡献图..."
        python3 "$SCRIPT_ROOT/plot_mmpbsa_decomp.py" "$DECOMP_DAT"
    else
        log_warning "未找到能量分解文件 $DECOMP_DAT，跳过贡献图绘制。"
    fi

    popd >/dev/null

}

# =============================================================================
# 主程序逻辑
# =============================================================================
main() {
    cd "$PDB_ROOT"
    for pdb_file in [!_]*.pdb; do
        [[ -f "$pdb_file" ]] || continue
        pdb_name=$(basename "$pdb_file" .pdb)
        [[ "$pdb_name" == *"_GMXMMPBSA_"* ]] && continue
        [[ -d "$pdb_name" ]] || continue

        log_info ">>> Processing: $pdb_name"
        for replica in "$pdb_name"/replica_*; do
            [[ -d "$replica" ]] && analyze_replica "$replica" "$pdb_name"
        done
        # 3. 数据汇总 (跨副本合并)
        merge_and_plot_system "$(realpath "$pdb_name")"
        
        # 4. 【Bug 修复】显式执行绘图调用
        plot_all_results "$(realpath "$pdb_name")"

        res_list=$(create_residue_index "match.tpr")
        analyze_residue_specifics "$analysis_dir" "$res_list"
        plot_residue_results "$merged_dir"
    done
    
    log_success "恭喜！所有体系的 500ns 模拟分析及绘图已圆满完成。"
}

main