# 2原子分子・H2O・CO2：結合距離・振動の力・3相・臨界点・電離（プラズマ）の比較
# 値は NIST Chemistry WebBook / Huber-Herzberg の標準値（コミット前に要照合）
import numpy as np
c=2.99792458e10; h=6.62607015e-34; kB=1.380649e-23; u=1.66053907e-27; eV=1.602176634e-19; NA=6.02214076e23; R=8.314462618
#        名前  m1   m2    re(Å)  ωe(cm-1) D0(eV) IE(eV) Ttp(K) Ptp(kPa) Tb(K) Tc(K)  Pc(MPa)
D=[("H2", 1.00783,1.00783,0.7414,4401.2,4.478,15.426,13.80, 7.04,  20.27, 33.145,1.2964),
   ("N2",14.00307,14.00307,1.0977,2358.6,9.759,15.581,63.15,12.52,  77.36,126.19,3.3958),
   ("O2",15.99491,15.99491,1.2075,1580.2,5.116,12.070,54.36, 0.146, 90.19,154.58,5.043),
   ("F2",18.99840,18.99840,1.4119, 916.6,1.602,15.697,53.48, 0.25,  85.03,144.41,5.172)]
print("分子 re(Å) k(N/m) ħω(eV) D0(eV) IE(eV) | Ttp  Tb   Tc(K) Pc(MPa) | kTc(eV) ħω/kTc σ(Å) σ/re")
for n,m1,m2,re,w,D0,IE,Ttp,Ptp,Tb,Tc,Pc in D:
    mu=m1*m2/(m1+m2)*u; k=mu*(2*np.pi*c*w)**2; hw=h*c*w/eV
    kTc=kB*Tc/eV; b=R*Tc/(8*Pc*1e6)            # vdW b (m3/mol)
    sig=(3*b/(2*np.pi*NA))**(1/3)*1e10         # b=(2/3)πNσ³
    print(f"{n:3s} {re:6.4f} {k:7.1f} {hw:6.3f} {D0:6.3f} {IE:6.2f} | {Ttp:5.1f} {Tb:5.1f} {Tc:6.2f} {Pc:6.3f} | {kTc:.4f} {hw/kTc:6.1f} {sig:5.2f} {sig/re:4.2f}")
# He（単原子）
Tc,Pc=5.1953,0.22746; b=R*Tc/(8*Pc*1e6); sig=(3*b/(2*np.pi*NA))**(1/3)*1e10
print(f"He  単原子  Tc={Tc} K Pc={Pc} MPa kTc={kB*Tc/eV:.5f} eV σ={sig:.2f}Å  （常圧で固体にならない、三重点なし）")
print()
# H2O, CO2
M=[("H2O","O-H 0.9575 Å, 角 104.48°",[3657,1595,3756],5.10,12.62,273.16,0.611657,373.12,647.096,22.064),
   ("CO2","C=O 1.1600 Å, 直線 180°",[1333,667,2349],5.45,13.777,216.58,517.95,194.7,304.13,7.3773)]
for n,geo,ws,D0,IE,Ttp,Ptp,Tb,Tc,Pc in M:
    b=R*Tc/(8*Pc*1e6); sig=(3*b/(2*np.pi*NA))**(1/3)*1e10
    hw=[h*c*w/eV for w in ws]
    print(f"{n}: {geo}; 振動 ħω={['%.3f'%x for x in hw]} eV; 最初の結合を切る D0≈{D0} eV; IE={IE} eV")
    print(f"     三重点 {Ttp} K / {Ptp} kPa, 常圧の沸点(昇華点) {Tb} K, 臨界 {Tc} K / {Pc} MPa, kTc={kB*Tc/eV:.4f} eV, σ={sig:.2f} Å")
print("\n正四面体の角 =",np.degrees(np.arccos(-1/3)))
