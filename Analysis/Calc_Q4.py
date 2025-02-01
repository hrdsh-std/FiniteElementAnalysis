import numpy as np

class Calc_Q4:
    @staticmethod
    def calc_K(fem_model,D,gps):
        mesh = fem_model.mesh
        K = np.zeros((len(mesh.nodes) * 2, len(mesh.nodes) * 2))
        for key,elem in mesh.elements.items():
            Ke = Calc_Q4.calc_Ke(elem,D,gps)
            for i in range(4):
                for j in range(4):
                    K[2 * elem.node[i]:2 * elem.node[i]+2,
                        2 * elem.node[j]:2 * elem.node[j]+2] += Ke[2 * i:2 * i+2, 2 * j:2 * j+2]
        return K
    
    @staticmethod
    def calc_Ke(elem,D, gps):
        """
        4nodeのKeマトリクスを算出。

        Args:
            elem (Element): 対象のElementクラス
            
        Returns:
            Ke (array 8x8): 要素剛性マトリクス
        """  
        Ke = np.zeros((8, 8))
        for xi, eta, wi, wj in gps:
            B, J = Calc_Q4.calc_B(elem, xi, eta)
            Ke += wi * wj * np.dot(B.T, np.dot(D, B)) * np.linalg.det(J) * elem.section.thickness
        return Ke

    @staticmethod
    def calc_B(elem, xi, eta):
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
        