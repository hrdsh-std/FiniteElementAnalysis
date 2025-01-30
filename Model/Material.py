class Material:
    """
    材料特性をまとめておくクラス
    
    Attributes:
        E (float): ヤング率
        nu (float): ポアソン比
    """

    def __init__(self, id,name, E, nu):
        self.id = id
        self.name = name
        self.E = E
        self.nu = nu

    def __str__(self):
        return "Material {}:name{}, E {}, Nu {}".format(self.id,self.name,self.E,self.nu)