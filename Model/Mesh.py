class Mesh:
    """
    NodeオブジェクトとElementオブジェクトをまとめておくクラス
    
    Attributes:
        nodes ({Node.id: Node}): Nodeの辞書
        elements ({Element.id: Element}): Elementの辞書
    """
    def __init__(self,nodes,elements,material):
        self.nodes = nodes
        self.elements = elements
        self.material = material
        self.get_element_coord()

    def get_element_coord(self):
        for key,elm in self.elements.items():
            elm.get_cordination(self.nodes)
