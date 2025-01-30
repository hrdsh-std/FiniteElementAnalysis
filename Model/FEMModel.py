import json
from Model.Material import Material
from Model.Section import Section
from Model.Node import Node
from Model.Element import Element
from Model.Mesh import Mesh
import numpy as np

class FEMModel:

    def __init__(self,json_file):
        with open(json_file,"r",encoding="utf-8") as f:
            self.data = json.load(f)

        self.materials = self.read_materials()
        self.sections = self.read_sections()
        self.nodes = self.read_nodes()
        self.ori4node = self.nodes
        self.elements = self.read_elements()
        self.read_boundary_conditions()
        self.read_loads()
        self.mesh = Mesh(self.nodes,self.elements,self.materials[0])#materialはElementが保持すべき？
        self.analysis_setting = self.read_analysis_settings()
        self.generate_mid_nodes()

    def read_materials(self):
        material_dict = {}
        for material_data in self.data["materials"]:
            material_dict[material_data["id"]] = Material(material_data["id"],material_data["name"],material_data["E"],material_data["nu"])
        return material_dict

    def read_sections(self):
        section_dict = {}
        for section_data in self.data["sections"]:
            section_dict[section_data["id"]] = Section(section_data["id"],section_data["thickness"])
        return section_dict
    
    def read_nodes(self):
        node_dict = {}
        for node_data in self.data["nodes"]:
            node_dict[node_data["id"]] = Node(node_data["id"],node_data["x"],node_data["y"])
        return node_dict

    def read_elements(self):
        element_dict = {}
        for element_data in self.data["elements"]:
            section = self.sections[element_data["section_id"]]
            element_dict[element_data["id"]] = Element(element_data["id"],element_data["nodes"],section)
        return element_dict

    def read_boundary_conditions(self):
        for bc_data in self.data["boundary_conditions"]:
            node = self.nodes[bc_data["node"]]
            node.set_bc_dist(bc_data["type"])

    def read_loads(self):
        for load_data in self.data["loads"]:
            node = self.nodes[load_data["node"]]
            node.set_bc_force(load_data["value"])

    def read_analysis_settings(self):
        return self.data["analysis_settings"][0]

    def generate_mid_nodes(self):
        if self.analysis_setting["element_type"] != "Quad_8node":
            return
        
        mid_nodes = {}  # 中間節点の座標を管理 {id:(x,y)}
        next_node_id = max(self.nodes.keys()) + 1  # 新規節点のID

        for key,elem in self.elements.items():
            corner_nodes = [self.nodes[id] for id in elem.node[:4]]  # コーナー節点（前4つ）
            mid_point_ids = []
            mid_point_xy = []
            node_dict = {}
            # 各辺の中間節点を生成
            for i in range(4):
                n1 = corner_nodes[i]
                n2 = corner_nodes[(i + 1) % 4]
                mid_x = (n1.x + n2.x) / 2
                mid_y = (n1.y + n2.y) / 2
                
                # 座標がすでに登録されている場合は共有
                if (mid_x, mid_y) in mid_nodes:
                    mid_node_id = mid_nodes[(mid_x, mid_y)]
                else:
                    mid_node_id = next_node_id
                    self.nodes[mid_node_id] = Node(mid_node_id, mid_x, mid_y)
                    mid_nodes[(mid_x, mid_y)] = mid_node_id
                    next_node_id += 1

                mid_point_ids.append(mid_node_id)
                mid_point_xy.append(np.array([[mid_x, mid_y]]))
                node_dict = {id_: coord for id_,coord in zip(mid_point_ids,mid_point_xy)}

            # 要素の節点リストを更新（コーナー節点 + 中間節点）
            elem.nodes = elem.node[:4] + mid_point_ids
            elem.add_node(node_dict)