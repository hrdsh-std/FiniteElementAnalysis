import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
import matplotlib.tri as tri

class Visualize:

    def __init__(self,solver):
         self.solver = solver

    def plot_deform(self,ratio=50):
            x = np.array([[node.x, node.y] for idx, node in self.solver.fem_model.ori4node.items()])
            x_new = x + self.solver.d.reshape(len(self.solver.fem_model.ori4node), 2) * ratio
            fig, ax = plt.subplots()
            for v in self.solver.mesh.elements:
                patch = patches.Polygon(xy=v.xy[:4], ec='black', alpha=0.3, fill=False,linestyle='--')
                ax.add_patch(patch)
            for v in self.solver.mesh.elements:
                xy_new = [(x_new[idx, 0], x_new[idx, 1]) for idx in v.node[:4]]
                patch = patches.Polygon(xy=xy_new, fc='red', ec='red', alpha=0.3, fill=False)
                ax.add_patch(patch)
            ax.autoscale()
            ax.set_aspect('equal', 'box')
            plt.show()  

    def plot_stress(self, stress_idx=3, ratio=20):
        """
        応力プロットを作成。

        Args:
            stress_idx (int(0~3)): 0-x応力, 1-y応力, 2-せん断応力, 3-von Mises
            ratio (float): 変形の倍率
        """
        x = np.array([[node.x, node.y] for idx, node in self.solver.fem_model.ori4node.items()])
        x_new = x + self.solver.d.reshape(len(self.solver.fem_model.ori4node), 2) * ratio
        nodal_values = []
        for i in range(len(self.solver.fem_model.ori4node)):
            nodal_values.append(self.solver.stress_dict[i][stress_idx])

        nodes_x = x_new[:, 0]
        nodes_y = x_new[:, 1]
        elements_tris = []
        for key,elm in self.solver.mesh.elements.items():
            elements_tris.append([elm.node[0], elm.node[1], elm.node[2]])
            elements_tris.append([elm.node[0], elm.node[2], elm.node[3]])
        triangulation = tri.Triangulation(nodes_x, nodes_y, elements_tris)
        fig, ax = plt.subplots()
        result = ax.tricontourf(triangulation, nodal_values,cmap='jet')
        for key,v in self.solver.mesh.elements.items():#要素境界の表示
            xy_new = [(x_new[idx, 0], x_new[idx, 1]) for idx in v.node[:4]]
            patch = patches.Polygon(xy=xy_new,ec='grey', fill=False)
            ax.add_patch(patch)
        ax.autoscale()
        ax.set_aspect('equal', 'box')
        ax.axis('off')
        fig.colorbar(result, ax=ax)
        max_dx,max_dy = self.solver.calc_max_d()
        fig.text(0.1, 0.1,f"max_dx:{max_dx:.4f} \nmax_dy:{max_dy:.4f}")

        plt.show()
