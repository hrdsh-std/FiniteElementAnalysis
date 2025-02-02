from Model.FEMModel import FEMModel
from Model.Visualize import Visualize as ModelV
from Analysis.Solver import Solver
from Result.Visualize import Visualize as ResultV
import json

if __name__ == "__main__":
     json_file = r"C:\Users\syuhe\Desktop\02_開発\02_Finite Element Analysis\inputFromGH.json"
    #json = r"C:\Users\syuhe\Desktop\02_開発\02_Finite Element Analysis\8-node quadrilateral element\input\input.json"
     with open(json_file,"r",encoding="utf-8") as f:
        data = json.load(f)
        fem_model = FEMModel(data)
        #print(fem_model.analysis_setting)
        #modelv = ModelV(fem_model)
        #modelv.plot_mesh()
        solver = Solver(fem_model)
        visualize = ResultV(solver)
        visualize.plot_stress()