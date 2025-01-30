class Section:
    """
    断面の情報をまとめておくクラス
    
    Attributes:
        id (int): id
        thickness (float): 厚さ
    """

    def __init__(self,id,thickness):
        self.id = id
        self.thickness = thickness

    def __str__(self):
        return "Section {}: thickness {}".format(self.id,self.thickness)
