import numpy as np
from scipy.optimize import brentq
k=8.9875517923e9; e=1.602176634e-19; me=9.1093837e-31; hbar=1.054571817e-34
a0=hbar**2/(me*k*e**2); ke2=k*e**2
# 陽子 z=±R/2、電子は中央面の半径ρの円上で向かい合う
def F_proton(R,rho):  # 陽子1個にかかる軸方向の力（+＝外向き＝反発）
    return ke2/R**2 - 2*ke2*(R/2)/(rho**2+R**2/4)**1.5
def F_electron(R,rho): # 電子1個にかかる半径方向の力（+＝外向き）
    return -2*ke2*rho/(rho**2+R**2/4)**1.5 + ke2/(2*rho)**2
# 1) 陽子の釣り合い ρ(R)
for R_A in (0.5,0.583,0.74):
    R=R_A*1e-10; rho=np.sqrt(3)/2*R
    print(f"R={R_A}Å ρ=√3/2R={rho*1e10:.4f}Å  F_proton={F_proton(R,rho):.2e} N")
# 2) 電子の円運動 + 角運動量 ħ
def cond(rho):
    R=2*rho/np.sqrt(3); return -F_electron(R,rho) - hbar**2/(me*rho**3)
rho=brentq(lambda x: cond(x*a0),0.1,5,xtol=1e-14)*a0; R=2*rho/np.sqrt(3)
print(f"ρ={rho/a0:.4f}a0={rho*1e10:.4f}Å  R={R/a0:.4f}a0={R*1e10:.4f}Å")
print(f"釣り合い点の力: 陽子に {ke2/R**2:.4e} N（反発）と {2*ke2*(R/2)/(rho**2+R**2/4)**1.5:.4e} N（電子の引力）")
print(f"電子の向心力 {-F_electron(R,rho):.4e} N")
# 結合エネルギー
v=hbar/(me*rho); T=2*0.5*me*v**2
V=ke2/R + ke2/(2*rho) - 4*ke2/np.sqrt(rho**2+R**2/4)
Etot=(T+V)/e; print(f"H2全エネルギー {Etot:.3f} eV, 2H=-27.21 eV → 結合エネルギー {-(Etot+27.211):.3f} eV（実測4.75 eV、ゼロ点込み4.52）")
