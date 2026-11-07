import numpy as np, itertools
phi=(1+5**0.5)/2
# 正二十面体（辺 2）の頂点
V=[]
for s1 in (1,-1):
  for s2 in (1,-1):
    V+= [(0,s1,s2*phi),(s1,s2*phi,0),(s2*phi,0,s1)]
V=np.array(V,float)
E=[(i,j) for i,j in itertools.combinations(range(12),2) if abs(np.linalg.norm(V[i]-V[j])-2)<1e-9]
F=[t for t in itertools.combinations(range(12),3) if all(((a,b) in E) or ((b,a) in E) for a,b in itertools.combinations(t,2))]
C=np.array([V[list(t)].mean(0) for t in F])          # 20面の中心 → 正十二面体の20頂点
d=sorted(set(np.round([np.linalg.norm(C[i]-C[j]) for i,j in itertools.combinations(range(20),2)],9)))
print("正二十面体: 頂点",len(V),"辺",len(E),"面",len(F))
print("面の中心 20 点の中心からの距離:",set(np.round(np.linalg.norm(C,axis=1),9)))
print("面の中心どうしの最短距離（十二面体の辺）:",d[0], " 各点の最近接数:",[int(sum(abs(np.linalg.norm(C-C[i],axis=1)-d[0])<1e-9)) for i in range(3)])
R_ico=np.linalg.norm(V[0]); r_dod=np.linalg.norm(C[0])
print(f"二十面体の外接半径 {R_ico:.4f}, 内側の十二面体の外接半径 {r_dod:.4f}, 比 {r_dod/R_ico:.4f}")
print(f"辺の比（十二面体/二十面体） {d[0]/2:.4f} = φ/3 ? {phi/3:.4f}")
# ペロブスカイト ABX3 の許容因子 t=(rA+rX)/(√2(rB+rX))（Shannon のイオン半径、Å）
rO=1.40
for name,rA,rB,note in [("CaTiO3",1.34,0.605,"灰チタン石（八面体が傾く）"),("SrTiO3",1.44,0.605,"立方の基準"),("BaTiO3",1.61,0.605,"強誘電体（Ti が中心からずれる）")]:
    t=(rA+rO)/(2**0.5*(rB+rO)); print(f"{name}: t = {t:.3f}  {note}")
