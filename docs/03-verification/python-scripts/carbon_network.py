# なぜ C だけが強固な固体（網目）になり、N・O・F は2原子分子の気体になるか
# 結合エネルギー（kJ/mol、標準的な平均値）→ eV
k=1/96.485
E={"C":(347,614,839),"N":(163,418,945),"O":(146,498,None),"F":(159,None,None)}
lone={"C":0,"N":1,"O":2,"F":3}; v={"C":4,"N":3,"O":2,"F":1}
print("元素 線の数(結合) 孤立電子対 単結合 二重 三重 | 網目(単結合のみ) 分子(多重結合)  勝つ方   二重/単  三重/単")
for X in "CNOF":
    s,d,t=E[X]; s*=k; d=d*k if d else None; t=t*k if t else None
    net=v[X]/2*s                          # 1原子あたり：v本の単結合を2原子で分け合う
    mol={4:t,3:t,2:d,1:s}[v[X]]/2         # 2原子分子：結合次数vの結合1本を2原子で分ける（Cは三重で比較）
    win="網目（固体）" if net>mol*1.01 else ("分子（気体）" if mol>net*1.01 else "同じ（線1本では網目は作れない）")
    r2=f"{d/s:.2f}" if d else "-"; r3=f"{t/s:.2f}" if t else "-"
    print(f" {X}   {v[X]}            {lone[X]}       {s:.2f}  {d if d is None else round(d,2)}  {t if t is None else round(t,2)} |  {net:.2f} eV         {mol:.2f} eV      {win}  {r2}  {r3}")
import numpy as np
a=3.5668
print(f"\nダイヤモンド：格子定数 a={a} Å、C–C 距離 a√3/4 = {a*np.sqrt(3)/4:.4f} Å、結合角 {np.degrees(np.arccos(-1/3)):.2f}°")
print("ダイヤモンドの凝集エネルギー 7.37 eV/原子（網目の見積り 7.19 eV と比較）")
print(f"充填率 π√3/16 = {np.pi*np.sqrt(3)/16:.3f}（最密充填 0.740）")
print("グラファイト：層内 C–C 1.42 Å（線3本、120°）、層間 3.35 Å（弱いファンデルワールス力）")
