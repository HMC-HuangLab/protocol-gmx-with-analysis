#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
残基级别接触数和最小距离分析脚本
从 GROMACS 轨迹中计算每个残基与配体的接触数和最小距离
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import subprocess
import re

def analyze_residue_contacts(traj_dir, pdb_name, ligand_name='LIG'):
    """
    分析每个残基与配体的接触情况
    
    参数:
        traj_dir: 轨迹目录 (包含 analysis 子目录)
        pdb_name: PDB 文件名 (不含扩展名)
        ligand_name: 配体名称
    """
    
    analysis_dir = os.path.join(traj_dir, 'analysis')
    os.makedirs(analysis_dir, exist_ok=True)
    
    tpr_file = os.path.join(traj_dir, f'results/prod/md_{pdb_name}_prod.tpr')
    xtc_file = os.path.join(traj_dir, f'results/prod/md_{pdb_name}_clean_nowat.xtc')
    
    if not os.path.exists(tpr_file) or not os.path.exists(xtc_file):
        print(f"找不到轨迹文件，跳过 {traj_dir}")
        return
    
    print(f"分析 {traj_dir} 的残基级别接触...")
    
    # 使用 gmx mindist 计算残基级别的接触
    cmd = f"""
echo "Protein {ligand_name}" | gmx mindist -s {tpr_file} -f {xtc_file} -on {analysis_dir}/numcont_res.xvg -res -cutoff 0.5
"""
    
    try:
        subprocess.run(cmd, shell=True, check=True, capture_output=True)
        print(f"残基接触数计算完成：{analysis_dir}/numcont_res.xvg")
    except subprocess.CalledProcessError as e:
        print(f"gmx mindist 失败：{e}")
        return
    
    # 使用 gmx mindist 计算残基级别的最小距离
    cmd = f"""
echo "Protein {ligand_name}" | gmx mindist -s {tpr_file} -f {xtc_file} -od {analysis_dir}/mindist_res.xvg -res
"""
    
    try:
        subprocess.run(cmd, shell=True, check=True, capture_output=True)
        print(f"残基最小距离计算完成：{analysis_dir}/mindist_res.xvg")
    except subprocess.CalledProcessError as e:
        print(f"gmx mindist 失败：{e}")
        return
    
    # 处理结果并绘制图表
    plot_residue_analysis(analysis_dir)


def parse_xvg_residue_data(filename):
    """
    解析残基级别的 XVG 文件
    
    返回：DataFrame，包含 Residue, Value 列
    """
    if not os.path.exists(filename):
        return None
    
    data = []
    res_names = []
    
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith(('#', '@')):
                # 从注释行提取残基名称
                if line.startswith('#'):
                    match = re.search(r'Res\s+(\d+)\s+(\w+)', line)
                    if match:
                        res_names.append(f"{match.group(2)}{match.group(1)}")
                continue
            
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    res_num = int(float(parts[0]))
                    value = float(parts[1])
                    data.append([res_num, value])
                except ValueError:
                    continue
    
    if not data:
        return None
    
    df = pd.DataFrame(data, columns=['Residue_num', 'Value'])
    
    # 尝试添加残基名称
    if res_names and len(res_names) == len(data):
        df['Residue'] = res_names
    else:
        df['Residue'] = df['Residue_num'].apply(lambda x: f"Res{x}")
    
    return df


def plot_residue_analysis(analysis_dir):
    """
    绘制残基级别的分析图
    """
    
    # 读取 numcont 数据
    numcont_file = os.path.join(analysis_dir, 'numcont_res.xvg')
    mindist_file = os.path.join(analysis_dir, 'mindist_res.xvg')
    
    # 计算平均接触数/距离
    numcont_data = []
    mindist_data = []
    
    if os.path.exists(numcont_file):
        with open(numcont_file, 'r') as f:
            for line in f:
                if line.startswith(('#', '@')):
                    continue
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        res_num = int(float(parts[0]))
                        occupancy = float(parts[1])
                        numcont_data.append([res_num, occupancy])
                    except ValueError:
                        continue
    
    if os.path.exists(mindist_file):
        with open(mindist_file, 'r') as f:
            for line in f:
                if line.startswith(('#', '@')):
                    continue
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        res_num = int(float(parts[0]))
                        dist = float(parts[1])
                        mindist_data.append([res_num, dist])
                    except ValueError:
                        continue
    
    # 绘制 numcont 柱状图
    if numcont_data:
        df_numcont = pd.DataFrame(numcont_data, columns=['Residue_num', 'Occupancy'])
        
        # 按占有率排序，选择前 20 个残基
        df_sorted = df_numcont.sort_values(by='Occupancy', ascending=False).head(20)
        
        fig, ax = plt.subplots(figsize=(12, 10), dpi=300)
        
        y_pos = range(len(df_sorted))
        colors = 'steelblue'
        
        ax.barh(y_pos, df_sorted['Occupancy'], color=colors, edgecolor='black', linewidth=0.5)
        
        ax.set_xlabel('Contact Occupancy (%)', fontsize=16)
        ax.set_title(f'Residue Contact Analysis - {os.path.basename(analysis_dir)}', fontsize=18)
        ax.set_yticks(y_pos)
        ax.set_yticklabels([f"Res{r}" for r in df_sorted['Residue_num']], fontsize=12)
        ax.invert_yaxis()
        ax.grid(True, alpha=0.3, linestyle='--', axis='x')
        
        plt.tight_layout()
        output_file = os.path.join(analysis_dir, 'numcont_residue_bar.png')
        plt.savefig(output_file, dpi=1000, bbox_inches='tight')
        print(f"残基接触图已保存：{output_file}")
        plt.show()
        
        # 保存为 CSV
        df_sorted.to_csv(os.path.join(analysis_dir, 'numcont_residue.csv'), index=False)
    
    # 绘制 mindist 柱状图
    if mindist_data:
        df_mindist = pd.DataFrame(mindist_data, columns=['Residue_num', 'Distance'])
        
        # 按距离排序 (从小到大)，选择前 20 个残基
        df_sorted = df_mindist.sort_values(by='Distance', ascending=True).head(20)
        
        fig, ax = plt.subplots(figsize=(12, 10), dpi=300)
        
        y_pos = range(len(df_sorted))
        colors = 'coral'
        
        ax.barh(y_pos, df_sorted['Distance'], color=colors, edgecolor='black', linewidth=0.5)
        
        ax.set_xlabel('Minimum Distance (nm)', fontsize=16)
        ax.set_title(f'Residue Minimum Distance Analysis - {os.path.basename(analysis_dir)}', fontsize=18)
        ax.set_yticks(y_pos)
        ax.set_yticklabels([f"Res{r}" for r in df_sorted['Residue_num']], fontsize=12)
        ax.invert_yaxis()
        ax.grid(True, alpha=0.3, linestyle='--', axis='x')
        
        plt.tight_layout()
        output_file = os.path.join(analysis_dir, 'mindist_residue_bar.png')
        plt.savefig(output_file, dpi=1000, bbox_inches='tight')
        print(f"残基最小距离图已保存：{output_file}")
        plt.show()
        
        # 保存为 CSV
        df_sorted.to_csv(os.path.join(analysis_dir, 'mindist_residue.csv'), index=False)


def merge_replica_residue_analysis(pdb_name, output_dir='merged_analysis'):
    """
    合并所有 replica 的残基级别分析结果
    """
    
    os.makedirs(output_dir, exist_ok=True)
    
    # 收集所有 replica 的数据
    replica_dirs = [d for d in os.listdir('.') if d.startswith('replica_') and pdb_name in d]
    
    numcont_data = {}
    mindist_data = {}
    
    for replica_dir in replica_dirs:
        analysis_dir = os.path.join(replica_dir, 'analysis')
        replica_num = replica_dir.split('_')[1]
        
        # 读取 numcont
        numcont_csv = os.path.join(analysis_dir, 'numcont_residue.csv')
        if os.path.exists(numcont_csv):
            df = pd.read_csv(numcont_csv)
            numcont_data[f'Replica_{replica_num}'] = df.set_index('Residue_num')['Occupancy']
        
        # 读取 mindist
        mindist_csv = os.path.join(analysis_dir, 'mindist_residue.csv')
        if os.path.exists(mindist_csv):
            df = pd.read_csv(mindist_csv)
            mindist_data[f'Replica_{replica_num}'] = df.set_index('Residue_num')['Distance']
    
    # 合并 numcont 数据
    if numcont_data:
        df_merged = pd.DataFrame(numcont_data)
        df_merged.index.name = 'Residue_num'
        df_merged.to_csv(os.path.join(output_dir, 'numcont_residue_merged.csv'))
        
        # 绘制合并图
        plot_merged_residue_numcont(df_merged, output_dir)
    
    # 合并 mindist 数据
    if mindist_data:
        df_merged = pd.DataFrame(mindist_data)
        df_merged.index.name = 'Residue_num'
        df_merged.to_csv(os.path.join(output_dir, 'mindist_residue_merged.csv'))
        
        # 绘制合并图
        plot_merged_residue_mindist(df_merged, output_dir)


def plot_merged_residue_numcont(df, output_dir):
    """
    绘制合并的残基 numcont 柱状图
    """
    fig, ax = plt.subplots(figsize=(14, 10), dpi=300)
    
    # 计算平均占有率
    mean_occupancy = df.mean(axis=1).sort_values(ascending=False).head(20)
    
    colors = ['#a1c9f4', '#ffd1dc', '#baffc9', '#ffdfba', '#ffffba', '#ff9999']
    
    y_pos = range(len(mean_occupancy))
    ax.barh(y_pos, mean_occupancy.values, color=colors[0], edgecolor='black', linewidth=0.5)
    
    ax.set_xlabel('Contact Occupancy (%)', fontsize=16)
    ax.set_title('Merged Residue Contact Analysis (All Replicas)', fontsize=18)
    ax.set_yticks(y_pos)
    ax.set_yticklabels([f"Res{r}" for r in mean_occupancy.index], fontsize=12)
    ax.invert_yaxis()
    ax.grid(True, alpha=0.3, linestyle='--', axis='x')
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'numcont_residue_merged.png')
    plt.savefig(output_file, dpi=1000, bbox_inches='tight')
    print(f"合并残基接触图已保存：{output_file}")
    plt.show()


def plot_merged_residue_mindist(df, output_dir):
    """
    绘制合并的残基 mindist 柱状图
    """
    fig, ax = plt.subplots(figsize=(14, 10), dpi=300)
    
    # 计算平均最小距离
    mean_dist = df.mean(axis=1).sort_values(ascending=True).head(20)
    
    colors = ['#a1c9f4', '#ffd1dc', '#baffc9', '#ffdfba', '#ffffba', '#ff9999']
    
    y_pos = range(len(mean_dist))
    ax.barh(y_pos, mean_dist.values, color=colors[1], edgecolor='black', linewidth=0.5)
    
    ax.set_xlabel('Minimum Distance (nm)', fontsize=16)
    ax.set_title('Merged Residue Minimum Distance Analysis (All Replicas)', fontsize=18)
    ax.set_yticks(y_pos)
    ax.set_yticklabels([f"Res{r}" for r in mean_dist.index], fontsize=12)
    ax.invert_yaxis()
    ax.grid(True, alpha=0.3, linestyle='--', axis='x')
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'mindist_residue_merged.png')
    plt.savefig(output_file, dpi=1000, bbox_inches='tight')
    print(f"合并残基最小距离图已保存：{output_file}")
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("用法：python plot_residue_analysis_full.py <traj_dir> <pdb_name> [ligand_name]")
        print("示例：python plot_residue_analysis_full.py . FAK_Karmadock_top_3k LIG")
        sys.exit(1)
    
    traj_dir = sys.argv[1]
    pdb_name = sys.argv[2]
    ligand_name = sys.argv[3] if len(sys.argv) > 3 else 'LIG'
    
    analyze_residue_contacts(traj_dir, pdb_name, ligand_name)
