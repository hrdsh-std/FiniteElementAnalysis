import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

integration_points = np.array([
    [-np.sqrt(3/5), -np.sqrt(3/5),5/9, 5/9],
    [0,-np.sqrt(3/5) ,8/9, 5/9],
    [np.sqrt(3/5), -np.sqrt(3/5),5/9, 5/9],
    [-np.sqrt(3/5), 0,5/9, 8/9],
    [0, 0,8/9, 8/9],
    [np.sqrt(3/5), 0,5/9, 8/9],
    [-np.sqrt(3/5), np.sqrt(3/5),5/9, 5/9],
    [0, np.sqrt(3/5),8/9, 5/9],
    [np.sqrt(3/5), np.sqrt(3/5),5/9, 5/9]
]) #(xi,eta,wi,wj)

# 積分点での応力（例: σ_x の値を仮定）
sigma_values = np.array([100, 110, 120,
                         130, 140, 150,
                         160, 170, 180])  # 任意の応力データ

A = np.column_stack([
    np.ones(len(integration_points)),  # 定数項
    integration_points[:, 0],          # x
    integration_points[:, 1],          # y
    integration_points[:, 0]**2,       # x^2
    integration_points[:, 1]**2,       # y^2
    integration_points[:, 0] * integration_points[:, 1]  # xy
])

# 最小二乗法で係数 a を求める
coefficients, _, _, _ = np.linalg.lstsq(A, sigma_values, rcond=None)

# 🔹 3Dプロットのための座標グリッドを作成（-1 ≦ x, y ≦ 1 の範囲）
x_vals = np.linspace(-1, 1, 30)
y_vals = np.linspace(-1, 1, 30)
X, Y = np.meshgrid(x_vals, y_vals)

# 🔹 2次多項式の応答曲面を計算
Z = (coefficients[0] + coefficients[1] * X + coefficients[2] * Y + 
     coefficients[3] * X**2 + coefficients[4] * Y**2 + coefficients[5] * X * Y)

# 🔹 3Dプロットの作成
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')

# 🔹 応答曲面の描画
ax.plot_surface(X, Y, Z, cmap='viridis', alpha=0.6)

# 🔹 積分点の散布図をプロット
ax.scatter(integration_points[:, 0], integration_points[:, 1], sigma_values, 
           color='red', s=50, label="積分点の応力")

# 🔹 ラベルを設定
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("応力 σ")
ax.set_title("最小二乗法による応答曲面の可視化")
ax.legend()

plt.show()

# 4つの隅部節点の座標（ξ, η）
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

sigma_nodes = A_corners @ coefficients

for i, sigma in enumerate(sigma_nodes):
    print(f"節点 {i+1} (ξ={corner_nodes[i,0]}, η={corner_nodes[i,1]}): σ = {sigma:.2f}")
    
    