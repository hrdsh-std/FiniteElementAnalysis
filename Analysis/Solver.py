import numpy as np
from Analysis.CalcStiffness import CalcStifness

class Solver:
    """
    有限要素法の計算を行うクラス
    
    Attributes:
        mesh (Mesh): Meshクラス
        material (Material): Materialクラス
        D (array 3x3): Dマトリクス
    """
    def __init__(self,fem_model):
        self.fem_model = fem_model
        self.mesh = self.fem_model.mesh
        self.f , self.bc_dist_dict = self.read_boundary_cond()
        #self.D = self.calc_D()
        #self.K = self.calc_K()
        #self.d = self.calc_d()

        self.gps = CalcStifness.set_gps(self.fem_model)
        self.D = CalcStifness.calc_D(self.fem_model)
        self.K = CalcStifness.calc_K(self.fem_model,self.D,self.gps)
        self.d = CalcStifness.calc_d(self.fem_model,self.f , self.bc_dist_dict,self.K)
        self.dx = self.d[::2]
        self.dy = self.d[1::2]
        self.stress_dict = self.calc_stress()
        
    def read_boundary_cond(self):
        """
        境界条件を読み込み。

        Returns:
            f (array 2node_num): 全体荷重ベクトル
            bc_dist_dict (dict): 全体変位ベクトルのうち、変位拘束のあるindをまとめた辞書
        """
        bc_dist_dict = {}
        f = []
        for k,v in self.mesh.nodes.items():
            f += v.bc_force
            if v.bc_dist[0] is not None:
                bc_dist_dict[k * 2] = v.bc_dist[0]
            if v.bc_dist[1] is not None:
                bc_dist_dict[k * 2 + 1] = v.bc_dist[1]
        return np.array(f), bc_dist_dict

    def calc_d(self):
        for i in range(len(self.mesh.nodes) * 2):
            if i in self.bc_dist_dict.keys():
                self.K[i, :] = np.zeros(len(self.mesh.nodes) * 2)
                self.K[:, i] = np.zeros(len(self.mesh.nodes) * 2)
                self.K[i, i] = 1.0
        d = np.dot(np.linalg.inv(self.K), self.f.T)
        return d

    def calc_stress(self):
        """
        節点応力を算出。

        Args:
            elem (Element): 対象のElementクラス

        Returns:
            stress_e (array 3x1): 節点応力
        """
        stress_dict = {i : np.zeros(3) for i in range(len(self.fem_model.ori4node))}
        counts_dict = {i : 0 for i in range(len(self.fem_model.ori4node))}

        for key,elem in self.mesh.elements.items():
            nodes_ind = []
            for ind in elem.node:
                nodes_ind.append(ind * 2)
                nodes_ind.append(ind * 2 + 1)
                d_elem = self.d[nodes_ind]

            stress_e = []
            for i, gp in enumerate(self.gps):
                xi, eta = gp[0], gp[1]
                B = CalcStifness.calc_B(self.fem_model,elem, xi, eta,self.D,self.gps) #プロパティに設定してもよい
                stress_e.append(np.dot(np.dot(self.D, B), d_elem.T))#積分点での応力がもとまる(9)x3materix
            #最小二乗法で節点力を算出
            
            node_stress_e = self.calc_stress_by_lstsq(stress_e,self.gps)
            
            #節点位置での平均応力を算出
            
            for i in range(4):
                stress_dict[elem.cornernode[i]] += node_stress_e[i]
                counts_dict[elem.cornernode[i]] += 1

        for i in range(len(self.fem_model.ori4node)):
            stress_dict[i] = stress_dict[i]  / counts_dict[i]
        
        
        ##mises応力
        for key,stress in stress_dict.items():
            mises_stress = (stress[0]**2-stress[0]*stress[1]+stress[1]**2+3*stress[2]**2)**(1/2)
            stress_dict[key] = np.append(stress,mises_stress)

        return stress_dict
    
    def calc_max_d(self):
        """
        最大変位を算出。
        Returns:
            max_dx,max_dy: 最大変位
        """
        max_dx = max(self.dx)
        max_dy = max(self.dy)

        return max_dx,max_dy

    def calc_stress_by_lstsq(self,stress_e,gps):
        integration_points = np.array(gps) #(xi,eta,wi,wj)
        stress_e = np.array(stress_e)
        sigX = np.array(stress_e[:,0])
        sigY = np.array(stress_e[:,1])
        tauXY = np.array(stress_e[:,2])

        A = np.column_stack([
            np.ones(len(integration_points)),  # 定数項
            integration_points[:, 0],          # x
            integration_points[:, 1],          # y
            integration_points[:, 0]**2,       # x^2
            integration_points[:, 1]**2,       # y^2
            integration_points[:, 0] * integration_points[:, 1]  # xy
        ])

        # 最小二乗法で係数 a を求める
        coefficients_sigX, _, _, _ = np.linalg.lstsq(A, sigX, rcond=None)
        coefficients_sigY, _, _, _ = np.linalg.lstsq(A, sigY, rcond=None)
        coefficients_tauXY, _, _, _ = np.linalg.lstsq(A, tauXY, rcond=None)

        corner_nodes = np.array([
            [-1, -1], [1, -1], [1, 1], [-1, 1]
        ])
        
        A_corners = np.column_stack([
            np.ones(len(corner_nodes)),  
            corner_nodes[:, 0],  
            corner_nodes[:, 1],  
            corner_nodes[:, 0]**2,  
            corner_nodes[:, 1]**2,  
            corner_nodes[:, 0] * corner_nodes[:, 1]
        ])

        sigX_nodes = A_corners @ coefficients_sigX
        sigY_nodes = A_corners @ coefficients_sigY
        tauXY_nodes = A_corners @ coefficients_tauXY

        node_stress_e = np.column_stack([sigX_nodes,sigY_nodes,tauXY_nodes])
        node_stress_e = node_stress_e.tolist()
        
        return node_stress_e #4x3
        