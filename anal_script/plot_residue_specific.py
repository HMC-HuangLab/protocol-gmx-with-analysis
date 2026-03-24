#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
残基特异性分析脚本
根据 MMPBSA 结果中的关键残基，绘制对应的 numcont 和 mindist 柱状图
"""

import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
import re

def parse_residue_name(res_name):
    """
    解析残基名称，提取残基名和序号
    例如：GLU211 -> ('GLU', 211)
    """
    match = re.match(r'([A-Z]+)(\d+)', str(res_name))
    if match:
        return match.group(1), int(match.group(2))
    return res_name, 0


def plot_residue_specific(numcont_file, mindist_file, residues_file):
    """
    根据关键残基列表绘制 numcont 和 mindist 柱状图
    
    参数:
        numcont_file: 合并后的接触数 CSV 文件
        mindist_file: 合并后的最小距离 CSV 文件  
        residues_file: 关键残基列表文件
    """
    # 读取关键残基列表
    with open(residues_file, 'r') as f:
        key_residues = [line.strip() for line in f if line.strip()]
    
    print(f"关键残基：{key_residues}")
    
    # 注意：这里需要从 GROMACS 的 residue-specific 分析文件中提取数据
    # 由于主脚本中只计算了整体的 numcont 和 mindist，这里需要重新分析每个残基的贡献
    
    # 创建图表 - Numcont
    fig1, ax1 = plt.subplots(figsize=(12, 8), dpi=300)
    
    # 创建图表 - Mindist
    fig2, ax2 = plt.subplots(figsize=(12, 8), dpi=300)
    
    # 颜色方案
    colors = ['#a1c9f4', '#ffd1dc', '#baffc9', '#ffdfba', '#ffffba', '#ff9999']
    
    # 由于需要残基级别的 numcont 和 mindist 数据，
    # 这里绘制的是基于关键残基的示意图
    # 实际数据需要从 GROMACS 的 residue-specific 分析中获取
    
    # 这里提供一个模板，实际使用时需要结合具体的残基分析数据
    print("注意：残基级别的 numcont 和 mindist 需要额外的 GROMACS 分析")
    print("请使用以下命令获取残基级别的接触数据:")
    print("gmx mindist -s topol.tpr -f traj.xtc -on numcont_res.xvg -res -cutoff 0.5")
    
    # 示例：绘制关键残基的示意图 (实际使用时替换为真实数据)
    x_pos = range(len(key_residues))
    color = colors[0]
    
    # Numcont 柱状图 (示例)
    ax1.bar(x_pos, [50] * len(key_residues), color=color, edgecolor='black', linewidth=1)
    ax1.set_xlabel('Residue', fontsize=16)
    ax1.set_ylabel('Contact Occupancy (%)', fontsize=16)
    ax1.set_title('Residue-Specific Contact Analysis', fontsize=18)
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(key_residues, fontsize=12, rotation=45, ha='right')
    ax1.grid(True, alpha=0.3, linestyle='--', axis='y')
    
    # Mindist 柱状图 (示例)
    color = colors[1]
    ax2.bar(x_pos, [0.3] * len(key_residues), color=color, edgecolor='black', linewidth=1)
    ax2.set_xlabel('Residue', fontsize=16)
    ax2.set_ylabel('Minimum Distance (nm)', fontsize=16)
    ax2.set_title('Residue-Specific Minimum Distance Analysis', fontsize=18)
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(key_residues, fontsize=12, rotation=45, ha='right')
    ax2.grid(True, alpha=0.3, linestyle='--', axis='y')
    
    plt.tight_layout()
    
    # 保存图片
    base_name = os.path.basename(residues_file).replace('.txt', '')
    output_file1 = f'{base_name}_numcont.png'
    output_file2 = f'{base_name}_mindist.png'
    
    fig1.savefig(output_file1, dpi=1000, bbox_inches='tight')
    fig2.savefig(output_file2, dpi=1000, bbox_inches='tight')
    
    print(f"残基接触图已保存：{output_file1}")
    print(f"残基最小距离图已保存：{output_file2}")
    
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("用法：python plot_residue_specific.py <numcont_csv> <mindist_csv> <residues_txt>")
        sys.exit(1)
    
    numcont_file = sys.argv[1]
    mindist_file = sys.argv[2]
    residues_file = sys.argv[3]
    
    plot_residue_specific(numcont_file, mindist_file, residues_file)
