import numpy as np
import Model.Visualize as modelv
from Model.Mesh import Mesh
from Model.Node import Node
from Model.Element import Element
from Model.Material import Material
from Analysis.Solver import Solver
import Result.Visualize as Resultv

import matplotlib.pyplot as plt
#定数
MESH_SIZE = 5  #mm
LENGTH = 100. #mm
HEIGHT = 10. #mm
THICKNESS = 10. #mm
E = 205000 #N/mm2
NU = 0.3 

FORCE = 10 #N

#材料定義
material = Material(E,NU)

#メッシュの作成
nodes, elems = modelv.create_mesh(LENGTH, HEIGHT, MESH_SIZE, THICKNESS)
mesh = Mesh(nodes, elems,material)
# for node in mesh.nodes.values():
#     print(node)
# for elem in mesh.elements:
#     print(elem)
# modelv.plot_mesh(mesh)



num_mesh_len = int(LENGTH / MESH_SIZE)
num_mesh_hei = int(HEIGHT / MESH_SIZE)
if num_mesh_len == 0:
    num_mesh_len = 1
if num_mesh_hei == 0:
    num_mesh_hei = 1

for k,v in mesh.nodes.items():
    # x変位の境界条件を設定
    if v.x == 0.0:
        v.bc_dist[0] = 0.0
    # y変位の境界条件を設定
    if v.x == 0.0:
        v.bc_dist[1] = 0.0
    # 荷重の設定
    if v.x == LENGTH and v.y == HEIGHT/2:
        v.bc_force = [0.0,FORCE / num_mesh_hei]

#解析実行
solver = Solver(mesh)

#ポスト処理
visualize = Resultv.Visualize(solver)
visualize.plot_stress()

def create_mesh(length, height, mesh_size, thickness):
    """
    モデルを格子点に分割し、全ての点のNodeクラスおよびそれらを使用した4角形1次要素を作成する.
    ただし、モデル形状は四角形とする。
    
    Args:
        length, height (float): モデルの長さ、高さ
        mesh_size (float): 分割したいメッシュのサイズ
        thickness (float): 要素の厚み(今回は全て一定とする)
        
    Returns:
        nodes_dict ({Node.id: Node}): 作成したNodeの辞書
        elsms (list): 作成したElementクラスをまとめたリスト
    """
    # Nodeの作成
    num_mesh_len = int(length / mesh_size)
    num_mesh_hei = int(height / mesh_size)
    if num_mesh_len == 0:
        num_mesh_len = 1
    if num_mesh_hei == 0:
        num_mesh_hei = 1
    x = np.linspace(0, length, num_mesh_len + 1)
    y = np.linspace(0, height, num_mesh_hei + 1)
    X, Y = np.meshgrid(x, y)
    X = X.ravel()
    Y = Y.ravel()

    nodes_dict = {}
    for i, coord in enumerate(zip(X, Y)):
        nodes_dict[i] = Node(i, coord[0], coord[1])
    
    # Elementの作成
    node_idx = 0
    elem_idx = 0
    elems = []
    for i in range(num_mesh_hei):
        for j in range(num_mesh_len + 1):
            if (node_idx + 1) % (num_mesh_len + 1) == 0:
                node_idx += 1
                continue
            else:
                node_idxs = [node_idx, node_idx + 1,
                         node_idx + num_mesh_len + 2, node_idx + num_mesh_len + 1]
                elems.append(Element(elem_idx, node_idxs, thickness))
                node_idx += 1
                elem_idx += 1
    return nodes_dict, elems