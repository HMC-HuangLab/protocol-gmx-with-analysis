#!/bin/bash
# =============================================================================
# Anal Test MD Workflow V40.0 (Individual Replica Plotting Edition)
# Fix: Generate Hbond and 2D FES images independently for each replica
# Fix: gmx_MMPBSA usage with the -nogui option
# Feature: Automated parsing of GLY:193 format; replica-specific color scheme system
# =============================================================================

set -euo pipefail

# --- [1] 环境配置 ---
PDB_ROOT=$(realpath "${1:-.}")
LIGNAME=${2:-LIG}
SCRIPT_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)
MMPBSA_TEMPLATE="${SCRIPT_ROOT}/mmpbsa.in"

GMX=${GMX:-gmx}
MPI_BIN=${MPI_BIN:-mpirun}
MMPBSA_BIN=${MMPBSA_BIN:-gmx_MMPBSA}

log_info() { echo -e "\033[0;34m[INFO]\033[0m $1"; }
log_success() { echo -e "\033[0;32m[SUCCESS]\033[0m $1"; }
log_error() { echo -e "\033[0;31m[ERROR]\033[0m $1"; }

# =============================================================================
# [2] 函数：analyze_residue_specifics (残基深挖)
# =============================================================================
analyze_residue_specifics() {
    local r_abs="$1"
    local rep_idx="$2"
    local analysis_dir="$r_abs/analysis"
    local target_dat="${analysis_dir}/FINAL_DECOMP_MMPBSA.dat"

    # 目标目录：merged_analysis/replica_X/
    local merged_base="$(dirname "$r_abs")/merged_analysis"
    local rep_merged_dir="${merged_base}/replica_${rep_idx}"
    local dist_merged_dir="${rep_merged_dir}/distance_plots"
    mkdir -p "$dist_merged_dir"

    [[ -f "$target_dat" ]] || { log_error "找不到分解文件: $target_dat"; return; }

    [[ -f "$target_dat" ]] || return
    pushd "$analysis_dir" >/dev/null
    mkdir -p res_plots

    # 颜色硬编码
    local cur_col="#a1c9f4"
    [[ "$rep_idx" -eq 2 ]] && cur_col="#ffd1dc"
    [[ "$rep_idx" -eq 3 ]] && cur_col="#baffc9"

    # 1. 提取残基列表 (ARG:137 格式)
    local res_info_list
    res_info_list=$(awk -F',' '/DELTAS:/,/Sidechain Energy Decomposition/ { if($1 ~ /^R::/) print $1 }' "$target_dat" | sed 's/R:://' | xargs || echo "")
    
    local res_count=$(echo "$res_info_list" | wc -w)
    log_info "副本 $rep_idx: 识别到 $res_count 个贡献残基，正在生成 1000 DPI 高清图..."

    # 2. 绘制 MMPBSA 能量分解图 (1000 DPI)
    # 2. 绘制当前副本的 MMPBSA 贡献图并移动到 merged_analysis
    if [[ -f "$SCRIPT_ROOT/plot_mmpbsa_decomp.py" ]]; then
        python3 "$SCRIPT_ROOT/plot_mmpbsa_decomp.py" "$target_dat"
        # 假设绘图脚本生成的文件名是 FINAL_DECOMP_MMPBSA_contribution.png
        if [[ -f "FINAL_DECOMP_MMPBSA_contribution.png" ]]; then
            mv "FINAL_DECOMP_MMPBSA_contribution.png" "${rep_merged_dir}/mmpbsa_energy_decomp_rep${rep_idx}.png"
            log_success "副本 $rep_idx 能量分解图已重命名保存至 merged_analysis"
        fi
    fi
    # mv mmpbsa_decomp_temp.png "${rep_merged_dir}/mmpbsa_energy_decomp_rep${rep_idx}.png"

    # 3. 生成索引 (使用 r 指令)
    { for info in $res_info_list; do echo "r $(echo $info | cut -d: -f2)"; done; echo "r $LIGNAME"; echo "q"; } | "$GMX" make_ndx -f match.tpr -o INdex4res.ndx > /dev/null 2>&1
    local lig_val=$(($(grep -c "^\[ " INdex4res.ndx) - 1))

    # 4. 全量 Mindist 循环
    for info in $res_info_list; do
        local r_name=$(echo "$info" | cut -d: -f1)
        local r_num=$(echo "$info" | cut -d: -f2)
        local res_lab="${r_name}${r_num}"
        local res_grp_line=$(grep "^\[ " INdex4res.ndx | grep -n -w "r_${r_num}" | head -n 1 | cut -d: -f1 || echo "")
        
        if [[ -n "$res_grp_line" ]]; then
            local res_val=$((res_grp_line - 1))
            echo "$res_val $lig_val" | "$GMX" mindist -s match.tpr -f match.xtc -n INdex4res.ndx \
                -od "res_plots/dist_${res_lab}.xvg" -on "res_plots/cont_${res_lab}.xvg" -d 0.25 -tu ns > /dev/null 2>&1
            
            python3 <<EOF
import pandas as pd, matplotlib.pyplot as plt, os
xvg = 'res_plots/dist_${res_lab}.xvg'
if os.path.exists(xvg):
    d = []
    with open(xvg) as f:
        for l in f:
            if not l.startswith(('#','@')):
                p = l.split()
                if len(p)>=2: d.append([float(p[0]), float(p[1])])
    if d:
        df = pd.DataFrame(d, columns=['T','D'])
        plt.figure(figsize=(8,4))
        plt.plot(df['T'], df['D'], color='$cur_col', lw=1.2)
        plt.title('$res_lab Interaction - Replica $rep_idx')
        plt.xlabel('Time (ns)'); plt.ylabel('Distance (nm)')
        plt.ylim(0, 1.2); plt.grid(alpha=0.3); plt.tight_layout()
        plt.savefig('dist_temp.png', dpi=1000)
        plt.close()
EOF
            mv dist_temp.png "${dist_merged_dir}/plot_dist_${res_lab}_rep${rep_idx}.png"
        fi
    done

    # 5. 副本 Occupancy 汇总柱状图
    python3 <<EOF
import pandas as pd, glob, matplotlib.pyplot as plt, numpy as np
occ_list = []
for f in sorted(glob.glob('res_plots/cont_*.xvg')):
    lab = f.split('cont_')[-1].replace('.xvg','')
    v = []
    with open(f) as o:
        for l in o:
            if not l.startswith(('#','@')):
                p = l.split(); v.append(float(p[1]))
    if v:
        occ = (np.count_nonzero(np.array(v) > 0) / len(v)) * 100
        occ_list.append({'Res': lab, 'Occ': occ})
if occ_list:
    df = pd.DataFrame(occ_list).sort_values('Occ', ascending=False)
    plt.figure(figsize=(14,6))
    plt.bar(df['Res'], df['Occ'], color='$cur_col', edgecolor='black', alpha=0.8)
    plt.ylabel('Occupancy (%)'); plt.title('Residue Contact Occupancy - Replica $rep_idx')
    plt.xticks(rotation=45, ha='right', fontsize=8)
    plt.ylim(0, 105); plt.tight_layout()
    plt.savefig('occ_temp.png', dpi=1000)
    plt.close()
EOF
    mv occ_temp.png "${rep_merged_dir}/summary_occupancy_rep${rep_idx}.png"
    popd >/dev/null
}


# =============================================================================
# [3] 函数：analyze_replica (核心分析与独立绘图)
# =============================================================================
analyze_replica() {
    local r_dir="$1"
    local rep_idx="$2"
    local r_abs=$(realpath "$r_dir")
    local analysis_dir="$r_abs/analysis"

    # 路径修复
    local merged_base="$(dirname "$r_abs")/merged_analysis"
    local rep_merged_dir="${merged_base}/replica_${rep_idx}"
    mkdir -p "$rep_merged_dir"

    log_info ">>> 处理副本 $rep_idx: $(basename "$r_abs")"

    # 同步拓扑至 prod
    local prod_dir="$r_abs/results/prod"
    for f in "topol.top" "posre_ligand.itp" "posre.itp" "ligand.itp" "index.ndx"; do
        [[ -f "$r_abs/$f" ]] && cp "$r_abs/$f" "$prod_dir/"
    done

    local prod_tpr=$(find "$prod_dir" -name "*.tpr" | head -n 1)
    local clean_xtc=$(find "$prod_dir" -name "*_clean_full.xtc" | head -n 1)
    [[ -z "$prod_tpr" || -z "$clean_xtc" ]] && return

    mkdir -p "$analysis_dir"
    pushd "$analysis_dir" >/dev/null

    # A. 轨迹精简
    "$GMX" make_ndx -f "$prod_tpr" -o raw.ndx <<EOF
1 | r $LIGNAME
name 22 Complex_Pure
q
EOF
    echo "Complex_Pure" | "$GMX" convert-tpr -s "$prod_tpr" -n raw.ndx -o match.tpr > /dev/null 2>&1
    echo "Complex_Pure" | "$GMX" trjconv -s "$prod_tpr" -f "$clean_xtc" -n raw.ndx -o match.xtc > /dev/null 2>&1
    "$GMX" make_ndx -f match.tpr -o Index4_anal.ndx <<EOF
r $LIGNAME
q
EOF

    # 2. 常规计算
    echo "3 3" | "$GMX" rms -s match.tpr -f match.xtc -o rmsd_protein_ca.xvg -n Index4_anal.ndx -tu ns > /dev/null 2>&1
    echo "3 13" | "$GMX" rms -s match.tpr -f match.xtc -o rmsd_ligand.xvg -n Index4_anal.ndx -tu ns > /dev/null 2>&1
    echo "4" | "$GMX" rmsf -s match.tpr -f match.xtc -o rmsf.xvg -res > /dev/null 2>&1
    echo "1" | "$GMX" gyrate -s match.tpr -f match.xtc -o rg_total.xvg > /dev/null 2>&1
    echo "1 13" | "$GMX" hbond -s match.tpr -f match.xtc -n Index4_anal.ndx -num hbnum.xvg -tu ns > /dev/null 2>&1
    echo "0" | "$GMX" sasa -s match.tpr -f match.xtc -o sasa.xvg -tu ns > /dev/null 2>&1

    # --- C. Hbond & FES 绘图逻辑 ---
    local cur_col="#a1c9f4"
    [[ "$rep_idx" -eq 2 ]] && cur_col="#ffd1dc"
    [[ "$rep_idx" -eq 3 ]] && cur_col="#baffc9"

    python3 <<EOF
import pandas as pd, numpy as np, matplotlib.pyplot as plt, os
def load_xvg(p, is_rg=False):
    d = []
    if os.path.exists(p):
        with open(p) as f:
            for l in f:
                if not l.startswith(('#','@')):
                    p = l.split(); d.append([float(p[0]), float(p[1])])
    if not d: return None
    df = pd.DataFrame(d, columns=['T','V'])
    if is_rg: df['T'] = df['T'] / 1000.0
    return df

# Hbond
hb = load_xvg('hbnum.xvg')
if hb is not None:
    plt.figure(figsize=(10,4))
    plt.plot(hb['T'], hb['V'], color='$cur_col', lw=1.2)
    plt.title('Hydrogen Bonds - Replica $rep_idx'); plt.xlabel('Time (ns)'); plt.ylabel('Count')
    plt.savefig('hb_temp.png', dpi=1000, bbox_inches='tight'); plt.close()

# FES 绘图增强版 (包含 Colorbar 和 1000 DPI 优化)
rmsd = load_xvg('rmsd_protein_ca.xvg')
rg = load_xvg('rg_total.xvg', is_rg=True)
if rmsd is not None and rg is not None:
    ml = min(len(rmsd), len(rg))
    x, y = rmsd['V'][:ml].values, rg['V'][:ml].values
    plt.figure(figsize=(8, 6.5))
    cnt, x_e, y_e = np.histogram2d(x, y, bins=40)
    fe = -np.log((cnt / np.sum(cnt)) + 1e-10); fe -= np.min(fe)
    
    # 绘图并捕获对象用于 Colorbar
    cf = plt.contourf(x_e[:-1], y_e[:-1], fe.T, levels=20, cmap='jet')
    
    # 关键修正：这里的 \$ 符号在 Bash 中必须被转义
    cbar = plt.colorbar(cf)
    cbar.set_label(r'\$\Delta G\$ (\$k_BT\$)', fontsize=12, labelpad=10)
    
    plt.title('2D Free Energy Surface - Replica $rep_idx', fontsize=14, pad=15)
    plt.xlabel('RMSD (nm)', fontsize=12); plt.ylabel('Radius of Gyration (nm)', fontsize=12)
    plt.tight_layout()
    plt.savefig('fes_temp.png', dpi=1000, bbox_inches='tight'); plt.close()
EOF

    mv hb_temp.png "${rep_merged_dir}/hbond_rep${rep_idx}.png" 2>/dev/null || true
    mv fes_temp.png "${rep_merged_dir}/fes2d_rep${rep_idx}.png" 2>/dev/null || true

    # C. MMPBSA (-nogui)
    if [[ -f "$MMPBSA_TEMPLATE" ]]; then
        cp "$MMPBSA_TEMPLATE" "mmpbsa.in"
        "$MPI_BIN" -np 14 "$MMPBSA_BIN" -O -i mmpbsa.in -cs match.tpr -ct match.xtc -ci Index4_anal.ndx \
            -cg 1 13 -cp "$prod_dir/topol.top" -o mmgbsa_results.dat -eo FINAL_DECOMP_MMPBSA.dat 
    fi
    popd >/dev/null

    # D. 深度挖掘
    analyze_residue_specifics "$r_abs" "$rep_idx"
}

# =============================================================================
# [4] 全局汇总与绘图 (Python)
# =============================================================================
merge_and_plot_all() {
    local pdb_path="$1"
    local merged_dir="$pdb_path/merged_analysis"
    mkdir -p "$merged_dir"
    log_info "聚合全局数据..."
    pushd "$pdb_path" >/dev/null

    find . -name "replica_*_hbond.png" -o -name "replica_*_fes2d.png" -exec cp {} "$merged_dir/" \;

    python3 <<EOF
import pandas as pd, glob, os
def read_xvg(path, is_rg=False):
    if not os.path.exists(path): return None
    data = []
    with open(path) as f:
        for l in f:
            if l.startswith(('#','@')): continue
            p = l.split()
            if len(p)>=2:
                t = float(p[0])/1000.0 if is_rg else float(p[0])
                data.append([t, float(p[1])])
    return pd.DataFrame(data) if data else None

tasks = [('rmsd_protein_ca.xvg', 'rmsd_merged.csv'), ('rmsd_ligand.xvg', 'lig_rmsd_merged.csv'),
         ('rmsf.xvg', 'rmsf_merged.csv'), ('rg_total.xvg', 'rg_merged.csv'), ('sasa.xvg', 'sasa_merged.csv')]

for pat, out in tasks:
    rg_f = True if 'rg_total' in pat else False
    files = sorted(glob.glob(f'replica_*/analysis/{pat}'))
    if files:
        df0 = read_xvg(files[0], is_rg=rg_f)
        if df0 is not None:
            x_col = 'Residue' if 'rmsf' in pat else 'Time'
            merged = pd.DataFrame({x_col: df0[0]})
            for f in files:
                d = read_xvg(f, is_rg=rg_f)
                if d is not None: merged[f.split('/')[0]] = d[1].values[:len(merged)]
            os.makedirs('merged_analysis', exist_ok=True); merged.to_csv(f'merged_analysis/{out}', index=False)
EOF

    # 搬运结果与绘图
    pushd "merged_analysis" >/dev/null
    [[ -f "rmsd_merged.csv" ]] && python3 "$SCRIPT_ROOT/plot_rmsd_analysis.py" rmsd_merged.csv "Protein RMSD" "RMSD (nm)"
    [[ -f "lig_rmsd_merged.csv" ]] && python3 "$SCRIPT_ROOT/plot_rmsd_analysis.py" lig_rmsd_merged.csv "Ligand RMSD" "RMSD (nm)"
    [[ -f "rmsf_merged.csv" ]] && python3 "$SCRIPT_ROOT/plot_rmsf_analysis.py" rmsf_merged.csv
    [[ -f "rg_merged.csv" ]] && python3 "$SCRIPT_ROOT/plot_rg_analysis.py" rg_merged.csv
    [[ -f "sasa_merged.csv" ]] && python3 "$SCRIPT_ROOT/plot_sasa_analysis.py" sasa_merged.csv

    # for f in $(find "$pdb_path" -name "plot_dist_*.png" -o -name "occupancy_replica_*.png"); do cp "$f" .; done
    popd >/dev/null; popd >/dev/null
}

# =============================================================================
# [5] Main
# =============================================================================
main() {
    cd "$PDB_ROOT"
    for pdb_file in [!_]*.pdb; do
        [[ -f "$pdb_file" ]] || continue
        pdb_name=$(basename "$pdb_file" .pdb)
        [[ "$pdb_name" == *"_GMXMMPBSA_"* ]] && continue
        [[ -d "$pdb_name" ]] || continue

        log_info "################################################################"
        log_info ">>> Processing: $pdb_name"
        
        rep_count=1
        for replica in "$pdb_name"/replica_*; do
            if [[ -d "$replica" ]]; then
                analyze_replica "$replica" "$rep_count"
                rep_count=$((rep_count + 1))
            fi
        done
        merge_and_plot_all "$(realpath "$pdb_name")"
    done
    log_success "所有分析已完成，请在各体系下的 merged_analysis 文件夹查看残基图表。"
}

main