import numpy as np
from Analysis.Calc_Q4 import Calc_Q4
from Analysis.Calc_Q8 import Calc_Q8
from Analysis.Calc_Q4Incomp import Calc_Q4Incomp

class CalcStifness:
    @staticmethod
    def calc_D(fem_model):
        """
        Dマトリクスを算出。

        Returns:
            D (array 3x3): Dマトリクス
        """
        mesh = fem_model.mesh
        D = np.array([
                    [1,mesh.material.nu,0],
                    [mesh.material.nu,1,0],
                    [0,0,(1-mesh.material.nu)/2]  
                        ]) * mesh.material.E / (1-mesh.material.nu**2)
        return D
    @staticmethod
    def calc_K(fem_model,D,gps):
        if fem_model.analysis_setting["element_type"] == "Quad_4node":
            return Calc_Q4.calc_K(fem_model,D,gps)
        elif fem_model.analysis_setting["element_type"] == "Quad_8node":
            return Calc_Q8.calc_K(fem_model,D,gps)
        elif fem_model.analysis_setting["element_type"] == "Quad_4node_Incomp":
            return  Calc_Q4Incomp.calc_K(fem_model,D,gps)
        else:
            return
    @staticmethod
    def calc_B(fem_model,elem, xi, eta):
        if fem_model.analysis_setting["element_type"] == "Quad_4node":
            return Calc_Q4.calc_B(elem, xi, eta)
        elif fem_model.analysis_setting["element_type"] == "Quad_8node":
            return Calc_Q8.calc_B(elem, xi, eta)
        else:
            return
    @staticmethod
    def set_gps(fem_model):
        elmtype = fem_model.analysis_setting["element_type"]
        gps = None
        if elmtype in ("Quad_4node","Quad_4node_Incomp"):
            gps = ((-1 / np.sqrt(3), -1 / np.sqrt(3), 1, 1),
            (1 / np.sqrt(3), -1 / np.sqrt(3), 1, 1),
            (1 / np.sqrt(3), 1 / np.sqrt(3), 1, 1),
            (-1 / np.sqrt(3), 1 / np.sqrt(3), 1, 1),
            ) #(xi,eta,wi,wj)
        elif elmtype == "Quad_8node":
            gps = (
            (-np.sqrt(3/5), -np.sqrt(3/5),5/9, 5/9),
            (0, -np.sqrt(3/5),8/9, 5/9),
            (np.sqrt(3/5), -np.sqrt(3/5),5/9, 5/9),
            (-np.sqrt(3/5), 0,5/9, 8/9),
            (0, 0,8/9, 8/9),
            (np.sqrt(3/5), 0,5/9, 8/9),
            (-np.sqrt(3/5), np.sqrt(3/5),5/9, 5/9),
            (0, np.sqrt(3/5),8/9, 5/9),
            (np.sqrt(3/5), np.sqrt(3/5),5/9, 5/9)
            ) #(xi,eta,wi,wj)
        return gps
    @staticmethod
    def calc_d(fem_model,f,bc_dist_dict,K):
        mesh = fem_model.mesh
        for i in range(len(mesh.nodes) * 2):
            if i in bc_dist_dict.keys():
                K[i, :] = np.zeros(len(mesh.nodes) * 2)
                K[:, i] = np.zeros(len(mesh.nodes) * 2)
                K[i, i] = 1.0
        d = np.dot(np.linalg.inv(K), f.T)
        return d
