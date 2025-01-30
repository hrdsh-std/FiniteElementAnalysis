
class Node:
    """
    Nodeのインデックス、座標および境界条件をまとめておくクラス(2次元用)

    Attributes:
        id (int): Nodeのインデックス番号(重複不可)
        x, y (float): Nodeのx, y座標
    """
    def __init__(self,idx,x,y):
        self.id = idx
        self.x = x
        self.y = y
        self.bc_dist = [None, None] # 節点拘束[x方向、y方向](None:拘束なし)
        self.bc_force = [0, 0] # 節点荷重[x方向、y方向]

    def __str__(self):
        return "Node {}: x {}, y {}".format(self.id,self.x,self.y)
    
    def set_bc_dist(self,xy):
        self.bc_dist = xy
    
    def set_bc_force(self,xy):
        self.bc_force = xy