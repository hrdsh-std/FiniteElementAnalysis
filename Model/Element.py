import numpy as np

class Element:
    """
    4角形1次要素Elementのインデックス、要素内のNode情報および要素の厚みをまとめておくクラス
    
    Attributes:
        id (int): Elementのインデックス番号(重複不可)
        node ([node_id, node_id, node_id, node_id]): 要素内のNodeのインデックスのリスト(左下から反時計回りで指定)
        sections ({id;thickness}}): 断面の辞書
        xy (nparray(size: 2, 4)): Element内のNodeのxy座標まとめ
    """

    def __init__(self,idx,node_idxs,section):
        self.id = idx
        self.node = node_idxs
        self.cornernode = node_idxs
        self.section = section
        self.xy = None
    
    def __str__(self):
        return "Element {}: Nodes {}, Thickness {}".format(self.id,self.node,self.section.thickness)
    
    def get_cordination(self,nodes_dict):
        """
        要素の各節点のx,y座標を求めてself.xyに代入。

        Args:
            nodes_dict ({Node.id: Node}): Nodeの辞書
        """
        res = []
        for node_idx in self.node:
            res.append([nodes_dict[node_idx].x,nodes_dict[node_idx].y])
        self.xy = np.array(res)

    def add_node(self,node_dict):
        """
        Nodeを追加

        Args:
            nodes_dict ({Node.id:np.array(x,y)}): Nodeの辞書
        """
        for key , xy in node_dict.items():
            self.node.append(key)
            self.xy = np.append(self.xy,xy,axis=0)