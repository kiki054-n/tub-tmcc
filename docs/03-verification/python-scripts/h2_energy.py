import numpy as np
k=8.9875517923e9; e=1.602176634e-19; me=9.1093837e-31; hbar=1.054571817e-34; c=2.99792458e8
ke2=k*e**2; a0=hbar**2/(me*ke2); Ha=ke2/a0/e
rho=4/(3*np.sqrt(3)-1)*a0; R=2*rho/np.sqrt(3); s=R; v=hbar/(me*rho)
T=2*0.5*me*v**2/e
V_pp=ke2/R/e; V_ee=ke2/(2*rho)/e; V_pe=-4*ke2/s/e
VE=V_pp+V_ee+V_pe
VB=ke2*v**2/(2*c**2*2*rho)/e   # Darwin項：逆向きに動く電子2個
tot=T+VE+VB
for n,x in [("運動エネルギー(力が電子にした仕事)",T),("電場: 陽子-陽子",V_pp),("電場: 電子-電子",V_ee),("電場: 陽子-電子×4",V_pe),("電場 合計",VE),("磁場: 電子-電子",VB),("全エネルギー",tot)]:
    print(f"{n:22s} {x:+12.6f} eV")
print(f"ビリアル T/|V_E| = {T/abs(VE):.4f}")
print(f"H原子2個 = {-Ha:.4f} eV → 結合エネルギー {-(tot+Ha):.4f} eV")
print(f"実測: H2 全電子エネルギー {-1.17447*Ha:.3f} eV, 結合エネルギー 4.748 eV")
