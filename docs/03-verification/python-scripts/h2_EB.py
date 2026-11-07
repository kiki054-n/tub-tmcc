import numpy as np
from scipy.optimize import brentq
k=8.9875517923e9; e=1.602176634e-19; me=9.1093837e-31; hbar=1.054571817e-34
c=2.99792458e8; mu0=4e-7*np.pi; ke2=k*e**2; a0=hbar**2/(me*ke2)
def solve(mag):
    def f(x):
        rho=x*a0; s=2*rho/np.sqrt(3); v=hbar/(me*rho)
        rep=ke2/(2*rho)**2*(1+(v/c)**2 if mag else 1)
        return 2*ke2*rho/s**3 - rep - hbar**2/(me*rho**3)
    rho=brentq(f,0.5,2,xtol=1e-15)*a0; return rho, 2*rho/np.sqrt(3)
rho,R=solve(False); v=hbar/(me*rho)
print(f"電場のみ: ρ={rho*1e10:.6f}Å R={R*1e10:.6f}Å  v/c={v/c:.5f}")
# 電子リングの磁場（軸上）
I=2*e*v/(2*np.pi*rho)
for z,name in [(0,"中心"),(R/2,"陽子の位置")]:
    print(f"B({name}) = {mu0*I*rho**2/(2*(rho**2+z**2)**1.5):.3f} T  向き: 陽子-陽子の軸線方向")
# 電子1にかかる力の分解
s=R; FE_p=2*ke2*rho/s**3; FE_e=ke2/(2*rho)**2
FB_e=mu0*e**2*v**2/(4*np.pi*(2*rho)**2)
print(f"電子: 電場の引力(陽子2個) {FE_p:.4e} N, 電場の反発(電子) {FE_e:.4e} N, 磁場の反発(電子) {FB_e:.4e} N  比 {FB_e/FE_e:.2e}")
print("陽子: 電場の軸方向の合力 0、磁場の力 0（陽子は止まっている：v=0）")
rho2,R2=solve(True)
print(f"磁場込み: R={R2*1e10:.6f}Å  変化 {(R2-R)/R:.2e}（実測 0.741 Å との差 {(0.741e-10-R)/R:.2%}）")
