#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2D 自由能景观图绘制脚本
基于 RMSD 和 Rg 绘制 2D 自由能景观
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys

def plot_fes_2d(csv_file):
    """
    绘制 2D 自由能景观图
    
    参数:
        csv_file: 合并后的 2D FES 数据 CSV 文件 (包含 Time, RMSD_nm, Rg_nm 列)
    """
    df = pd.read_csv(csv_file)
    
    # 提取 RMSD 和 Rg 数据
    rmsd = df['RMSD_nm'].values
    rg = df['Rg_nm'].values
    
    # 创建 2D 直方图 (用于计算自由能)
    nbins = 100
    H, xedges, yedges = np.histogram2d(rmsd, rg, bins=nbins, 
                                        range=[[rmsd.min(), rmsd.max()], 
                                               [rg.min(), rg.max()]])
    
    # 计算自由能 G = -kT * ln(P)
    # 假设 T = 300K, k = 0.008314 kJ/(mol·K)
    kT = 0.008314 * 300  # kJ/mol
    
    # 将计数转换为概率
    H = H / H.sum()
    
    # 避免 log(0)，添加小值
    H = np.ma.masked_array(H, mask=(H == 0))
    H = np.ma.masked_where(H <= 0, H)
    
    # 计算自由能 (kJ/mol)
    G = -kT * np.log(H)
    
    # 设置最小自由能为 0
    G = G - G.min()
    
    # 创建图表
    fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
    
    # 绘制 2D 自由能景观
    X, Y = np.meshgrid(xedges[:-1], yedges[:-1])
    im = ax.pcolormesh(X, Y, G.T, cmap='jet', shading='auto', 
                       norm=LogNorm(vmin=0.1, vmax=G.max()))
    
    # 添加颜色条
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Free Energy (kJ/mol)', fontsize=14)
    
    # 设置标签和标题
    ax.set_xlabel('RMSD (nm)', fontsize=16)
    ax.set_ylabel('Radius of Gyration Rg (nm)', fontsize=16)
    ax.set_title('2D Free Energy Landscape', fontsize=18)
    
    # 设置刻度字体
    ax.tick_params(axis='both', which='major', labelsize=14)
    
    plt.tight_layout()
    
    # 保存图片
    output_file = csv_file.replace('.csv', '.png')
    plt.savefig(output_file, dpi=1000, bbox_inches='tight')
    print(f"2D 自由能景观图已保存：{output_file}")
    
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("用法：python plot_fes_2d.py <csv_file>")
        sys.exit(1)
    
    csv_file = sys.argv[1]
    plot_fes_2d(csv_file)
