import numpy as np

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

    def calc_K(fem_model,D,gps):
        if fem_model.analysis_setting["element_type"] == "Quad_4node":
            return CalcStifness.calc_4node_K(fem_model,D,gps)
        elif fem_model.analysis_setting["element_type"] == "Quad_8node":
            return CalcStifness.calc_8node_K(fem_model,D,gps)
        else:
            return
        
    def calc_B(fem_model,elem, xi, eta):
        if fem_model.analysis_setting["element_type"] == "Quad_4node":
            return CalcStifness.calc_4node_B(elem, xi, eta)
        elif fem_model.analysis_setting["element_type"] == "Quad_8node":
            return CalcStifness.calc_8node_B(elem, xi, eta)
        else:
            return
            
    @staticmethod
    def set_gps(fem_model):
        elmtype = fem_model.analysis_setting["element_type"]
        gps = None
        if elmtype == "Quad_4node":
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

    def calc_4node_K(fem_model,D,gps):
        mesh = fem_model.mesh
        K = np.zeros((len(mesh.nodes) * 2, len(mesh.nodes) * 2))
        for key,elem in mesh.elements.items():
            Ke = CalcStifness.calc_4node_Ke(elem,D,gps)
            for i in range(4):
                for j in range(4):
                    K[2 * elem.node[i]:2 * elem.node[i]+2,
                      2 * elem.node[j]:2 * elem.node[j]+2] += Ke[2 * i:2 * i+2, 2 * j:2 * j+2]
        return K
    
    @staticmethod
    def calc_4node_Ke(elem,D, gps):
        """
        4nodeのKeマトリクスを算出。

        Args:
            elem (Element): 対象のElementクラス
            
        Returns:
            Ke (array 8x8): 要素剛性マトリクス
        """  
        Ke = np.zeros((8, 8))
        for xi, eta, wi, wj in gps:
            B, J = CalcStifness.calc_4node_B(elem, xi, eta)
            Ke += wi * wj * np.dot(B.T, np.dot(D, B)) * np.linalg.det(J) * elem.section.thickness
        return Ke
    
    @staticmethod
    def calc_4node_B(elem, xi, eta):
        """
        Bマトリクス及びヤコビアンを算出。

        Args:
            elem (Element): 対象のElementクラス
            xi, eta (float): Bマトリクス算出用の座標

        Returns:
            B (array 3x8): Bマトリクス
            J (array 2x2): ヤコビアン
        """
        dndxi = np.array([-1+eta,1-eta,1+eta,-1-eta])/4
        dndeta = np.array([-1+xi,-1-xi,1+xi,1-xi])/4

        x = elem.xy[:,0]
        y = elem.xy[:,1]

        dxdxi =  np.dot(dndxi,x)
        dydxi =  np.dot(dndxi,y)
        dxdeta = np.dot(dndeta,x)
        dydeta = np.dot(dndeta,y)

        J = np.array([[dxdxi,dydxi],[dxdeta,dydeta]])
        B = np.zeros((3,8))
        for i in range(4):
            Bi = np.dot(np.linalg.inv(J),np.array([dndxi[i],dndeta[i]]).T)
            B[0,2*i] = Bi[0]
            B[1,2*i+1] = Bi[1]
            B[2,2*i] = Bi[1]
            B[2,2*i+1] = Bi[0]
        return B , J
    
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
    
    @staticmethod
    def calc_8node_K(fem_model,D,gps):
        mesh = fem_model.mesh
        K = np.zeros((len(mesh.nodes) * 2, len(mesh.nodes) * 2))
        for key,elem in mesh.elements.items():
            Ke = CalcStifness.calc_8node_Ke(elem,D,gps)
            for i in range(8):
                for j in range(8):
                    K[2 * elem.node[i]:2 * elem.node[i]+2,
                      2 * elem.node[j]:2 * elem.node[j]+2] += Ke[2 * i:2 * i+2, 2 * j:2 * j+2]
        return K
    
    @staticmethod
    def calc_8node_Ke(elem,D, gps):
        """
        4nodeのKeマトリクスを算出。

        Args:
            elem (Element): 対象のElementクラス
            
        Returns:
            Ke (array 16x16): 要素剛性マトリクス
        """
        Ke = np.zeros((16, 16))
        for xi, eta, wi, wj in gps:
            B, J = CalcStifness.calc_8node_B(elem, xi, eta)
            Ke += wi * wj * np.dot(B.T, np.dot(D, B)) * np.linalg.det(J) * elem.section.thickness
        return Ke
    
    @staticmethod
    def calc_8node_B(elem, xi, eta):
        """
        Bマトリクス及びヤコビアンを算出。

        Args:
            elem (Element): 対象のElementクラス
            xi, eta (float): Bマトリクス算出用の座標

        Returns:
            B (array 3x8): Bマトリクス
            J (array 2x2): ヤコビアン
        """
        dndxi = np.array([(1-eta)*(2*xi+eta),(1-eta)*(2*xi-eta),(1+eta)*(2*xi+eta),(1+eta)*(2*xi-eta),\
                            -4*xi*(1-eta),2*(1-eta**2),-4*xi*(1+eta),-2*(1-eta**2)])/4
        dndeta = np.array([(1-xi)*(xi+2*eta),(1+xi)*(-xi+2*eta),(1+xi)*(xi+2*eta),(1-xi)*(-xi+2*eta),\
                            -2*(1-xi**2),-4*eta*(1+xi),2*(1-xi**2),-4*eta*(1-xi)])/4

        x = elem.xy[:,0]
        y = elem.xy[:,1]

        dxdxi =  np.dot(dndxi,x)
        dydxi =  np.dot(dndxi,y)
        dxdeta = np.dot(dndeta,x)
        dydeta = np.dot(dndeta,y)

        J = np.array([[dxdxi,dydxi],[dxdeta,dydeta]])
        B = np.zeros((3,16))
        for i in range(8):
            Bi = np.dot(np.linalg.inv(J),np.array([dndxi[i],dndeta[i]]).T)
            B[0,2*i] = Bi[0]
            B[1,2*i+1] = Bi[1]
            B[2,2*i] = Bi[1]
            B[2,2*i+1] = Bi[0]
        return B , J
     