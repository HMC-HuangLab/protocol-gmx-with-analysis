#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SASA (溶剂可及表面积) 分析绘图脚本
"""

import pandas as pd
import matplotlib.pyplot as plt
import sys

def plot_sasa_analysis(csv_file):
    """
    绘制 SASA 分析图
    
    参数:
        csv_file: 合并后的 SASA 数据 CSV 文件
    """
    df = pd.read_csv(csv_file)
    
    # 创建图表
    fig, ax = plt.subplots(figsize=(12, 7), dpi=300)
    
    # 颜色方案
    colors = ['#a1c9f4', '#ffd1dc', '#baffc9', '#ffdfba', '#ffffba', '#ff9999']
    
    # 绘制每个 replica 的数据
    columns = [col for col in df.columns if col != 'Time']
    
    for i, col in enumerate(columns):
        color = colors[i % len(colors)]
        ax.plot(df['Time'], df[col], label=col, color=color, linewidth=2, alpha=0.8)
    
    # 设置标签和标题
    ax.set_xlabel('Time (ns)', fontsize=16)
    ax.set_ylabel('SASA (nm²)', fontsize=16)
    ax.set_title('Solvent Accessible Surface Area Analysis', fontsize=18)
    
    # 设置刻度字体
    ax.tick_params(axis='both', which='major', labelsize=14)
    
    # 添加图例
    ax.legend(loc='best', fontsize=12, framealpha=0.9)
    
    # 添加网格
    ax.grid(True, alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    
    # 保存图片
    output_file = csv_file.replace('.csv', '.png')
    plt.savefig(output_file, dpi=1000, bbox_inches='tight')
    print(f"SASA 图已保存：{output_file}")
    
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("用法：python plot_sasa_analysis.py <csv_file>")
        sys.exit(1)
    
    csv_file = sys.argv[1]
    plot_sasa_analysis(csv_file)
