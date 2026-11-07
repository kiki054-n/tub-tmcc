# 電子を「雲（1s軌道、大きさ1/ζ）」として扱う水素分子（Heitler-London / Wang）
import numpy as np
from scipy.special import exp1
from scipy.optimize import minimize
Ha=27.211386; a0A=0.529177; g=np.euler_gamma
def E(R,z):
    D=z*R
    S=np.exp(-D)*(1+D+D**2/3); Sp=np.exp(D)*(1-D+D**2/3)
    Jb=1/D-np.exp(-2*D)*(1+1/D)          # ∫a²/r_b
    Kx=np.exp(-D)*(1+D)                   # ∫ab/r_a
    J2=1/D-np.exp(-2*D)*(1/D+11/8+3*D/4+D**2/6)
    K2=0.2*(-np.exp(-2*D)*(-25/8+23*D/4+3*D**2+D**3/3)
            +6/D*(S**2*(g+np.log(D))+Sp**2*(-exp1(4*D))-2*S*Sp*(-exp1(2*D))))
    haa=z**2/2 - z - z*Jb
    hab=(-z**2*S/2+z**2*Kx) - 2*z*Kx
    return 1/R+(2*haa+z*J2+2*S*hab+z*K2)/(1+S**2)
# 1) Heitler-London（ζ=1：電子雲の大きさは水素原子と同じ）
r=minimize(lambda x:E(x[0],1.0),[1.6]); R=r.x[0]
print(f"HL  ζ=1     R={R:.3f}a0={R*a0A:.3f}Å  結合E={-(r.fun+1)*Ha:.2f} eV  (文献 1.64a0, 3.14 eV)")
# 2) Wang（雲の大きさも最適化）
r=minimize(lambda x:E(x[0],x[1]),[1.4,1.17]); R,z=r.x
print(f"Wang ζ={z:.3f} R={R:.3f}a0={R*a0A:.3f}Å  結合E={-(r.fun+1)*Ha:.2f} eV  (文献 ζ=1.166, 1.41a0, 3.78 eV)")
print(f"雲の半径 a0/ζ = {a0A/z:.3f} Å（水素原子 0.529 Å）")
print("実測: R=0.741 Å, 結合E=4.75 eV")
# He（雲、有効電荷）
Zeff=27/16; print(f"He 雲: Zeff={Zeff:.4f} 半径 a0/Zeff={a0A/Zeff:.3f}Å  E={-Zeff**2*Ha:.2f} eV（実測 -79.00）")
