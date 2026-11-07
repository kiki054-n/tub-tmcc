# C（グラファイト・ダイヤモンド）と BN（hBN・cBN）：格子定数から結合距離・密度を計算
import numpy as np
u=1.66053907e-24  # g
mC=12.011; mBN=10.81+14.007
rows=[]
# 層状（六方）：結合距離 = a/√3、層間 = c/2、1層の単位胞に2原子、c方向に2層 → 単位胞4原子
for name,a,c,m_per_atom in [("グラファイト",2.4612,6.7079,mC),("hBN",2.5040,6.6612,mBN/2)]:
    d=a/np.sqrt(3); V=np.sqrt(3)/2*a**2*c*1e-24; rho=4*m_per_atom*u/V
    rows.append((name,"層状（線3本・120°）",d,c/2,rho,V/4*1e24))
# 立方（ダイヤモンド・閃亜鉛鉱）：結合距離 = a√3/4、単位胞8原子
for name,a,m_per_atom in [("ダイヤモンド",3.5668,mC),("cBN",3.6157,mBN/2)]:
    d=a*np.sqrt(3)/4; V=(a*1e-8)**3; rho=8*m_per_atom*u/V
    rows.append((name,"立方（線4本・109.47°）",d,None,rho,a**3/8))
print("物質        形                    結合距離(Å) 層間(Å)  密度(g/cm³) 1原子の体積(Å³)")
for n,f,d,l,rho,v in rows:
    print(f"{n:8s} {f:20s} {d:8.4f}  {('%.3f'%l) if l else '  -  ':>7s}  {rho:8.3f}   {v:7.3f}")
print("\n層状→立方 で結合距離が伸びる割合: C %.2f%%, BN %.2f%%"%((1.5445/1.4210-1)*100,(1.5656/1.4457-1)*100))
print("層状→立方 で密度が増える割合: C %.1f%%, BN %.1f%%"%((3.513/2.267-1)*100,(3.488/2.279-1)*100))
print("電気陰性度 B 2.04, C 2.55, N 3.04（差: C–C 0, B–N 1.00）")
