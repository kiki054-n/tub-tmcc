# 万有引力の質・量拡張方程式 v0.1 の検証
G=6.6743e-11; k=8.9875517923e9; e=1.602176634e-19; a0=5.29177210903e-11
mp=1.67262192e-27; me=9.1093837e-31; mn=1.67492750e-27; u=1.66053907e-27; eV=e
Mp,Me,Mn=1500,500,2000

# 0. 文書の定数を再現
G_M=G*mp*me/(Mp*Me); E_H=13.6057*eV
G_E=(k*e**2/a0**2)/(E_H**2/a0**2)
print("0. G_M=%.4e  G_E=%.4e (文書: 1.3559e-73, 4.8592e7)"%(G_M,G_E))

# 1. 重力：G_M を p-e で決めたとき p-p, n-n, e-e を予測
for name,M1,M2,m1,m2 in [("p-p",Mp,Mp,mp,mp),("n-n",Mn,Mn,mn,mn),("e-e",Me,Me,me,me)]:
    print("1. %s 予測/実際 = %.4g"%(name,(G_M*M1*M2)/(G*m1*m2)))

# 2. 等価原理：1 u あたりの「質」の組成依存
def M_per_u(Z,N,mass_u): return (Z*Mp+N*Mn+Z*Me)/mass_u
H=M_per_u(1,0,1.007825); Fe=M_per_u(26,30,55.934936); Be=M_per_u(4,5,9.012183); Ti=M_per_u(22,26,47.947942)
print("2. M/u: H=%.2f Fe=%.2f  差=%.3e"%(H,Fe,abs(H-Fe)/Fe))
print("   Be-Ti 差=%.3e（MICROSCOPE 上限 ~1e-15）"%(abs(Be-Ti)/Ti))

# 3. 水素様イオン Z=1..3：F_E=G_E*E^2/r^2 (E=電離エネルギー, r=a0/Z) と実際のクーロン力
for Z in (1,2,3):
    r=a0/Z; E=E_H*Z**2
    print("3. Z=%d 予測/実際 = %.3g"%(Z,(G_E*E**2/r**2)/(k*Z*e**2/r**2)))

# 4. E が距離で変わる場合 (E=ke^2/2r) の距離依存
for r in (0.5*a0,a0,2*a0):
    E=k*e**2/(2*r); print("4. r=%.1fa0 F_E/Coulomb = %.3g"%(r/a0,(G_E*E**2/r**2)/(k*e**2/r**2)))

# 5. 中性物質：E に符号がないと打ち消されない（H原子1対あたりの比）
print("5. F_E/F_重力 (H原子対, 同距離) = %.3e"%((G_E*E_H**2)/(G*(mp+me)**2)))
