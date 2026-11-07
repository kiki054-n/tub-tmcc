import numpy as np
hbar=1.054571817e-34; c=2.99792458e8; G=6.67430e-11; k=8.9875517923e9; e=1.602176634e-19
me=9.1093837e-31; mp=1.67262192e-27
lP=np.sqrt(hbar*G/c**3); mP=np.sqrt(hbar*c/G); FP=c**4/G; alpha=k*e**2/(hbar*c)
print(f"ℓP={lP:.5e} m  F_P=c⁴/G={FP:.5e} N  α={alpha:.9f} 1/α={1/alpha:.6f}")
a0=hbar**2/(me*k*e**2)
cases=[("H a0",a0,1),("H2 実測 R",0.741e-10,1),("H2 ボーア模型 R",0.5825e-10,1),("He 電子半径",0.3024e-10,2),("He 電子間 2r",0.6048e-10,1)]
for name,r,Z in cases:
    n=r/lP; F=Z*alpha*FP/n**2
    print(f"{name:16s} n={n:.6e}  F=Zα·F_P/n²={F:.4e} N  (ke²Z/r²={Z*k*e**2/r**2:.4e})")
print(f"a0/ℓP = (mP/me)/α = {mP/me/alpha:.6e}")
print(f"重力 p-e at a0: F_P·(mp/mP)(me/mP)/n² = {FP*(mp/mP)*(me/mP)/(a0/lP)**2:.4e} N")
print(f"整数化による相対ずれの上限 ~ 1/n = {lP/0.741e-10:.1e}")
