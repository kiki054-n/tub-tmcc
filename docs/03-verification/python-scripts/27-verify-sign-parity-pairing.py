"""
TUB-TMCC 27：符号の交互と対形成（電子・核子）の検証（2026-11-09）

この記録の数値はすべてこのスクリプトで出す。
  pip install periodictable
  python3 27-verify-sign-parity-pairing.py

1. 符号の交互：n 個の ±1 の積を (-1)^n と置いたときの点の数（2^(n-1)）
2. 電子：分子の電子数の偶奇と、対になっていない電子（不対電子）の数
3. 核子：天然に存在する核種を、陽子数 Z・中性子数 N の偶奇で4つに分けて数える
"""
import itertools
import math

import periodictable as pt

print("=" * 60)
print("1. 符号の交互（積 = (-1)^n）")
for n in range(2, 8):
    pts = [v for v in itertools.product((1, -1), repeat=n) if math.prod(v) == (-1) ** n]
    print(f"  n={n}: 積={(-1) ** n:+d}  {'閉' if n % 2 == 0 else '開'}  点の数={len(pts)}")

print("=" * 60)
print("2. 電子：電子数の偶奇と不対電子の数（基底状態、文献値）")
# (分子, 電子数, 不対電子の数)
molecules = [
    ("H2", 2, 0), ("N2", 14, 0), ("H2O", 10, 0), ("CO2", 22, 0), ("CH4", 10, 0),
    ("O2", 16, 2),  # 偶数なのに不対電子2個（三重項）
    ("OH", 9, 1), ("CH3", 9, 1), ("NO", 15, 1), ("NO2", 23, 1), ("ClO2", 33, 1),
]
for name, ne, unpaired in molecules:
    parity = "奇" if ne % 2 else "偶"
    note = ""
    if ne % 2 == 0 and unpaired:
        note = "  ← 偶数なのに開いている（反例）"
    print(f"  {name:5s} 電子 {ne:2d}（{parity}）  不対電子 {unpaired}{note}")
print("  奇数個の電子は、数え方から必ず1個以上が対にならない（論理的に必然で、検証にはならない）。")
print("  検証になるのは偶数側：偶数でも O2 は閉じない。")

print("=" * 60)
print("3. 核子：天然に存在する核種（periodictable の存在比 > 0）")
groups = {"偶-偶": [], "偶-奇": [], "奇-偶": [], "奇-奇": []}
for el in pt.elements:
    if el.number == 0:
        continue
    for A in el.isotopes:
        if getattr(el[A], "abundance", 0) > 0:
            Z, N = el.number, A - el.number
            key = ("奇" if Z % 2 else "偶") + "-" + ("奇" if N % 2 else "偶")
            groups[key].append((f"{el.symbol}-{A}", Z, N))
total = sum(len(v) for v in groups.values())
for k, v in groups.items():
    print(f"  Z-N が {k}: {len(v):3d}  ({len(v) / total:6.1%})")
print(f"  合計 {total}")

oo = groups["奇-奇"]
print("  奇-奇の核種:", ", ".join(f"{s}(Z={z},N={n})" for s, z, n in oo))
radioactive_long_lived = {"K-40", "V-50", "La-138", "Lu-176"}
stable_oo = [s for s, _, _ in oo if s not in radioactive_long_lived]
print("  うち長寿命の放射性（K-40, V-50, La-138, Lu-176）を除くと:", ", ".join(stable_oo))
light = [(s, z, n) for s, z, n in oo if s in ("H-2", "Li-6", "B-10", "N-14")]
print("  軽い奇-奇はすべて Z = N:", all(z == n for _, z, n in light))
