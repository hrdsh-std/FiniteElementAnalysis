import numpy as np

class Calc_Q4Incomp:
    @staticmethod
    def calc_K(fem_model,D,gps):
        mesh = fem_model.mesh
        K = np.zeros((len(mesh.nodes) * 2, len(mesh.nodes) * 2))
        for key,elem in mesh.elements.items():
            Ke,_ = Calc_Q4Incomp.calc_KeBe(elem,0,0,D,gps)
            for i in range(4):
                for j in range(4):
                    K[2 * elem.node[i]:2 * elem.node[i]+2,
                        2 * elem.node[j]:2 * elem.node[j]+2] += Ke[2 * i:2 * i+2, 2 * j:2 * j+2]
        return K

    @staticmethod
    def calc_KeBe(elem, _xi, _eta,D,gps):
        """
        4nodeのKeマトリクスを算出。

        Args:
            elem (Element): 対象のElementクラス

        Returns:
            Ke (array 8x8): 要素剛性マトリクス
        """
        Kcc = np.zeros((8, 8))
        Kci = np.zeros((8,4))
        Kic = np.zeros((4,8))
        Kii = np.zeros((4,4))
        _Bc,_Bi ,_= Calc_Q4Incomp.calc_BJ(elem, _xi, _eta)
        for xi, eta, wi, wj in gps:
            Bc,Bi, J = Calc_Q4Incomp.calc_BJ(elem, xi, eta)
            Kcc += wi * wj * (Bc.T @ D @ Bc) * np.linalg.det(J) * elem.section.thickness
            Kci += wi * wj * (Bc.T @ D @ Bi) * np.linalg.det(J) * elem.section.thickness
            Kic += wi * wj * (Bi.T @ D @ Bc) * np.linalg.det(J) * elem.section.thickness
            Kii += wi * wj * (Bi.T @ D @ Bi) * np.linalg.det(J) * elem.section.thickness

        Ke = Kcc - Kci @ np.linalg.inv(Kii) @ Kic
        Be  =  _Bc - _Bi @ np.linalg.inv(Kii) @ Kic
        return Ke ,Be

    def calc_BJ(elem, xi, eta):
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

        Bc = np.zeros((3,8))
        invJ = np.linalg.inv(J)
        for i in range(4):
            Bj = np.dot(invJ,np.array([dndxi[i],dndeta[i]]).T)
            Bc[0,2*i] = Bj[0]
            Bc[1,2*i+1] = Bj[1]
            Bc[2,2*i] = Bj[1]
            Bc[2,2*i+1] = Bj[0]


        Bi = np.array([
            [-2*xi*invJ[0,0],-2*eta*invJ[0,1],0,0],
            [0,0,-2*xi*invJ[1,0],-2*eta*invJ[1,1]],
            [-2*xi*invJ[1,0],-2*eta*invJ[1,1],-2*xi*invJ[0,0],-2*eta*invJ[0,1]]
        ])

        return Bc,Bi,J
