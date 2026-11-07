import numpy as np
k=8.9875517923e9; e=1.602176634e-19; me=9.1093837e-31; hbar=1.054571817e-34
ke2=k*e**2; a0=hbar**2/(me*ke2); Ry=ke2/(2*a0)/e
# He: 核(+2e)が中心、電子2個が半径rの円の反対側
r=a0/(7/4)
att=2*ke2/r**2; rep=ke2/(2*r)**2
print(f"He r={r/a0:.4f}a0={r*1e10:.4f}Å  引力{att:.4e} 反発{rep:.4e} 合計{att-rep:.4e} N")
E=-2*(7/4)**2*Ry; print(f"He 全エネルギー {E:.2f} eV（実測 -79.00）")
# H2 平衡点
rho=4/(3*np.sqrt(3)-1)*a0; R=2*rho/np.sqrt(3); s=R
att2=2*ke2/s**2*(rho/s); rep2=ke2/(2*rho)**2
print(f"H2 ρ={rho*1e10:.4f}Å R={R*1e10:.4f}Å 陽子-電子距離s={s*1e10:.4f}Å")
print(f"H2 引力(半径方向){att2:.4e} 反発{rep2:.4e} 合計{att2-rep2:.4e} N")
print("係数 He 7/4=1.75, H2", 3*np.sqrt(3)/4-1/4)
