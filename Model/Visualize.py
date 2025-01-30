import matplotlib.pyplot as plt
from matplotlib import patches
import numpy as np

class Visualize:
    def __init__(self,fem_model):
        self.fem_model = fem_model
        self.mesh = self.fem_model.mesh

    def plot_mesh(self):
        
        node_color = '#FFB6C1'  # 淡いピンク
        element_color = '#FFE4B5'  # クリーム色
        node_number_color = '#2F4F4F'  # ダークグリーン
        element_number_color = '#8A2BE2'  # 紫色
        
        fig, ax = plt.subplots()
        for key,elem in self.mesh.elements.items():
            patch = patches.Polygon(xy=elem.xy[:4],facecolor = element_color , edgecolor=node_color)
            ax.add_patch(patch)
            text_xy = np.mean(elem.xy, axis=0)
            ax.text(text_xy[0], text_xy[1], elem.id, fontsize=12,color=element_number_color, va='center', ha='center')
        for node in self.mesh.nodes.values():
            ax.scatter(node.x, node.y, fc=node_color, s=100)
            ax.text(node.x, node.y, node.id, fontsize=8, color=node_number_color, va='center', ha='center')
        ax.autoscale()
        ax.set_aspect('equal', 'box')
        plt.show()

