#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os

def plot_mmpbsa_decomp(dat_file):
    """解析 MMPBSA .dat 文件并绘制残基能量贡献图"""
    if not os.path.exists(dat_file):
        print(f"错误：找不到文件 {dat_file}")
        return

    results = []
    start_parsing = False
    
    with open(dat_file, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            # 定位到 DELTAS 下的 Total Energy Decomposition 部分
            if 'DELTAS:' in line:
                start_parsing = True
                continue
            
            if start_parsing and 'R::' in line and 'TOTAL' not in line:
                # 提取 Residue 和 TOTAL 的 Avg 值
                # 格式: R::ARG:137,Internal(avg,std,err),...,TOTAL(avg,std,err)
                parts = line.strip().split(',')
                res_id = parts[0].replace('R::', '') # 变为 ARG:137
                try:
                    # TOTAL Avg 位于倒数第三个元素
                    total_avg = float(parts[-3])
                    results.append({'Residue': res_id, 'Energy': total_avg})
                except (ValueError, IndexError):
                    continue
            
            # 如果遇到下一个大段落则停止（可选）
            if start_parsing and 'Sidechain Energy Decomposition' in line:
                break

    if not results:
        print("未在 DELTAS 截面中提取到有效残基数据。")
        return

    df = pd.DataFrame(results)
    # 按残基编号排序
    df['ResNum'] = df['Residue'].str.split(':').str[-1].astype(int)
    df = df.sort_values('ResNum')

    # 绘图
    plt.figure(figsize=(14, 6), dpi=300)
    # 结合能贡献通常为负值（越负贡献越大），使用不同颜色区分
    colors = ['#ff7f0e' if x > 0 else '#1f77b4' for x in df['Energy']]
    
    plt.bar(df['Residue'], df['Energy'], color=colors, edgecolor='black', alpha=0.8)
    
    plt.axhline(0, color='black', linewidth=1)
    plt.ylabel('Binding Energy Contribution (kcal/mol)', fontsize=14)
    plt.xlabel('Residue', fontsize=14)
    plt.title('Per-residue Energy Decomposition (DELTAS TOTAL)', fontsize=16)
    plt.xticks(rotation=45, ha='right', fontsize=10)
    plt.grid(axis='y', linestyle='--', alpha=0.4)
    
    plt.tight_layout()
    out_img = dat_file.replace('.dat', '_contribution.png')
    plt.savefig(out_img, bbox_inches='tight')
    print(f"能量分解贡献图已保存：{out_img}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("用法：python plot_mmpbsa_decomp.py <FINAL_DECOMP_MMPBSA.dat>")
    else:
        plot_mmpbsa_decomp(sys.argv[1])