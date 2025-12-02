# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 09:32:24 2025

@author: USER
"""

import numpy as np

def build_cell_index_map(Nx, Ny):
    """
    對每個 ix (cell 編號) 回傳它包含的原子全域 index。
    每個 cell 有 2*Ny 個原子，依照 lattice 生成順序：
    for ix in range(Nx):
        for iy in range(Ny):
            A(ix,iy)
            B(ix,iy)
    """
    cell_to_indices = []
    for ix in range(Nx):
        start = ix * (2 * Ny)
        end   = start + (2 * Ny)
        cell_to_indices.append(np.arange(start, end))
    return cell_to_indices  # list of length Nx, each is array of length (2*Ny)

def extract_blocks_from_ribbon(H_full, Nx, Ny,
                               left_cells=1,
                               center_cells=4,
                               right_cells=1):
    """
    切三塊：
    - 左 lead:    最左的 left_cells 個 cell
    - 中央 region: 中間 center_cells 個 cell
    - 右 lead:    最右的 right_cells 個 cell

    回傳：
      HL, HC, HR,
      VL_block, VR_block,
      left_idx, center_idx, right_idx,
      cell_map
    """

    assert Nx >= left_cells + center_cells + right_cells, "帶太短"

    cell_map = build_cell_index_map(Nx, Ny)

    left_range_start   = 0
    left_range_end     = left_cells
    center_range_start = left_range_end
    center_range_end   = left_range_end + center_cells
    right_range_start  = center_range_end
    right_range_end    = center_range_end + right_cells

    left_cells_list   = list(range(left_range_start,  left_range_end))
    center_cells_list = list(range(center_range_start, center_range_end))
    right_cells_list  = list(range(right_range_start, right_range_end))

    left_idx   = np.concatenate([cell_map[c] for c in left_cells_list])
    center_idx = np.concatenate([cell_map[c] for c in center_cells_list])
    right_idx  = np.concatenate([cell_map[c] for c in right_cells_list])

    HL = H_full[np.ix_(left_idx, left_idx)]
    HC = H_full[np.ix_(center_idx, center_idx)]
    HR = H_full[np.ix_(right_idx, right_idx)]

    # 嘗試抓 interface block（可能是 0，如果幾何沒有直接最近鄰）
    left_edge_idx    = cell_map[left_cells_list[-1]]
    center_left_idx  = cell_map[center_cells_list[0]]
    VL_block = H_full[np.ix_(left_edge_idx, center_left_idx)]

    center_right_idx = cell_map[center_cells_list[-1]]
    right_edge_idx   = cell_map[right_cells_list[0]]
    VR_block = H_full[np.ix_(center_right_idx, right_edge_idx)]

    return (HL, HC, HR,
            VL_block, VR_block,
            left_idx, center_idx, right_idx,
            cell_map)
