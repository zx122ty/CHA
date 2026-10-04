#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
01.py — 读取指定目录下所有 *.mzML 文件的文件名（含后缀），
      将结果保存为 CSV 文件（列名: file_name）。

用法:
    python 01.py <目录路径> [输出名称]

参数:
    <目录路径>    要扫描的目录（必填）
    [输出名称]    输出文件名的标识部分，默认为 "negative"

示例:
    python 01.py /home/data/shareData/16tdisk/Project/WuTouTang/ZX-WTD/成分
    python 01.py /path/to/dir positive
    python 01.py /path/to/dir my_custom_name
"""

import os
import sys
import csv
import glob

OUTPUT = 'negative'

def get_mzml_filenames(directory: str) -> list[str]:
    """
    递归扫描 directory 下所有 *.mzML 文件，返回文件名（含后缀）列表。
    按文件名自然排序（字母顺序）。
    """
    # 使用 glob 递归搜索所有 .mzML 文件
    pattern = os.path.join(directory, "**", "*.mzML")
    files = glob.glob(pattern, recursive=True)

    # 提取文件名（含后缀），不包含路径
    filenames = [os.path.basename(f) for f in files]

    # 去重并排序
    filenames = sorted(set(filenames))
    return filenames


def save_to_csv(filenames: list[str], output_path: str) -> None:
    """将文件名列表写入 CSV 文件，列名为 file_name。"""
    with open(output_path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow(["file_name"])
        for name in filenames:
            writer.writerow([name])
    print(f"✅ 已保存 {len(filenames)} 个文件名到: {output_path}")


def main():
    # 获取目标目录路径
    if len(sys.argv) < 2:
        print("❌ 请指定目录路径。")
        print(f"用法: python {os.path.basename(sys.argv[0])} <目录路径> [输出名称]")
        sys.exit(1)

    target_dir = sys.argv[1]

    # 获取输出名称（可选参数，默认为 'negative'）
    output_name = sys.argv[2] if len(sys.argv) >= 3 else OUTPUT

    # 检查目录是否存在
    if not os.path.isdir(target_dir):
        print(f"❌ 目录不存在或不是有效目录: {target_dir}")
        sys.exit(1)

    # 搜索 *.mzML 文件
    print(f"🔍 正在扫描目录: {target_dir}")
    filenames = get_mzml_filenames(target_dir)

    if not filenames:
        print("⚠️  未找到任何 *.mzML 文件。")
        sys.exit(0)

    print(f"📄 找到 {len(filenames)} 个 *.mzML 文件。")

    # 输出 CSV 文件路径（与脚本同目录）
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_csv = os.path.join(script_dir, f"mzml_file_list_{output_name}.csv")

    save_to_csv(filenames, output_csv)


if __name__ == "__main__":
    main()