#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MMPBSA 总结合能柱状图绘制脚本
"""

import pandas as pd
import matplotlib.pyplot as plt
import sys

def plot_mmpbsa_bar(csv_file):
    """
    绘制 MMPBSA 总结合能柱状图
    
    参数:
        csv_file: MMPBSA 总结果 CSV 文件 (包含 Replica, Component, Mean, Std 列)
    """
    df = pd.read_csv(csv_file)
    
    # 创建图表
    fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
    
    # 颜色 (与现有脚本风格一致 - 绿色)
    uniform_color = 'mediumseagreen'
    
    # 绘制柱状图
    replicas = df['Replica'].unique()
    x_pos = range(len(replicas))
    means = [df[df['Replica'] == r]['Mean'].values[0] for r in replicas]
    stds = [df[df['Replica'] == r]['Std'].values[0] for r in replicas]
    
    bars = ax.bar(x_pos, means, yerr=stds, color=uniform_color, 
                  capsize=5, edgecolor='black', linewidth=1.5)
    
    # 设置标签和标题
    ax.set_xlabel('Replica', fontsize=16)
    ax.set_ylabel('Binding Energy (kJ/mol)', fontsize=16)
    ax.set_title('MMPBSA Binding Free Energy', fontsize=18)
    
    # 设置刻度字体
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(replicas, fontsize=14)
    
    # 添加网格
    ax.grid(True, alpha=0.3, linestyle='--', axis='y')
    
    plt.tight_layout()
    
    # 保存图片
    output_file = csv_file.replace('.csv', '.png')
    plt.savefig(output_file, dpi=1000, bbox_inches='tight')
    print(f"MMPBSA 柱状图已保存：{output_file}")
    
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("用法：python plot_mmpbsa_bar.py <csv_file>")
        sys.exit(1)
    
    csv_file = sys.argv[1]
    plot_mmpbsa_bar(csv_file)
