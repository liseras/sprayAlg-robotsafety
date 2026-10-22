# -*- coding: utf-8 -*-
"""WPF Viewer 实际渲染截图对比：同一重建，行距 12 / 20 / 60。"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "SimSun"]
plt.rcParams["axes.unicode_minus"] = False

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    items = (
        ("ui_row12_iso.png", "行距 12mm（默认）：曲面 12 条行 / 全件 29 行"),
        ("ui_row20_iso.png", "行距 20mm：曲面 8 条行 / 全件 16 行"),
        ("ui_row60_iso.png", "行距 60mm：曲面 3 条行 / 全件 7 行"),
    )
    fig, axes = plt.subplots(1, 3, figsize=(19.5, 5.9))
    for ax, (fn, title) in zip(axes, items):
        ax.imshow(mpimg.imread(os.path.join(HERE, fn)))
        ax.set_title(title, fontsize=12)
        ax.axis("off")
    fig.suptitle("WPF Viewer 实测：行距 12/20/60 重建网格完全相同（同为 3 个面片），只有布行数随工艺行距变化",
                 fontsize=15)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out = os.path.join(HERE, "ui_compare.png")
    fig.savefig(out, dpi=130)
    print("saved", out)


if __name__ == "__main__":
    main()
