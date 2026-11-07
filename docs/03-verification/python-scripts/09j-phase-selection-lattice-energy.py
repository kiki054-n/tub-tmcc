#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
09j：相の選択をエネルギー比較で決められるか（格子エネルギーの必要条件と格子線の数え方）

09i までの統一エネルギー関数（σ なし、A・B 共通）では、
「格子軸として正四面体の 4 軸を使うか、XYZ の 6 軸を使うか」がまだ仮定として残っている。
ここでは同じエネルギー関数・同じ定数（09i 案(b)：H2O・H2S で決定）で、
各分子を四面体格子と直交格子の両方で最小化し、エネルギーを比べる。

[1] ΔU = U_min(直交格子) − U_min(四面体格子)
[2] 直交格子が格子エネルギー λ L^p だけ有利と仮定し、
    周期2（NH3, H2O）が四面体、周期3以降（6 分子）が直交になる λ の範囲を p ごとに求める
[3] p = 1（回復力の項 −L Σcos δ と同じスケーリング）で、TTT の格子線の数え方ごとに
    格子線 1 本あたりのエネルギー ε L の範囲を求める（ε は臂 1 本の整列エネルギー L に対する比）
CH4 は直交格子に 4 本の H を同等に置けないため比較から除く。
09i のスクリプトを読み込み、その geometry・定数の決め方をそのまま使う。
"""

import importlib.util
import io
import contextlib
import math
import os
import numpy as np
from scipy.optimize import minimize, fsolve

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "m09i", os.path.join(HERE, "09i-lp-tetra-angle-rule-fbi-frame.py"))
m = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(m)

VARIANT = "b"
A, B = fsolve(lambda p: [m.solve(n, p[0], p[1], VARIANT) - m.MOL[n][3] for n in m.REF],
              [426.0, 722.0], xtol=1e-10)
NAMES = ["NH3", "H2O", "PH3", "H2S", "AsH3", "H2Se", "SbH3", "H2Te"]
TETRA_SET = {"NH3", "H2O"}


def u_min(name, ortho):
    _, el, nH, _ = m.MOL[name]
    L = m.arm(el)
    axes, lps, w = m.geometry(nH, ortho, VARIANT)

    def U(x):
        H = [m.unit(x[3*i:3*i+3]) for i in range(nH)]
        u = -L * sum(float(np.dot(H[i], axes[i])) for i in range(nH))
        u += (A / L) * sum(1 / (2 * math.sin(m.gam(H[i], H[j]) / 2))
                           for i in range(nH) for j in range(i + 1, nH))
        u += (B / L) * w * sum(1 / (2 * math.sin(m.gam(h, lp) / 2)) for h in H for lp in lps)
        return u

    r = minimize(U, np.concatenate(axes), method="BFGS", options={"gtol": 1e-11, "maxiter": 20000})
    return r.fun, L, nH


data = {}
for n in NAMES:
    ut, L, nH = u_min(n, False)
    uo, _, _ = u_min(n, True)
    data[n] = (L, uo - ut, nH)


def window(weight):
    """weight(L, nH) を掛けた格子エネルギー差が相を正しく分ける係数の範囲 (下限, 上限)"""
    lo = max(dU / weight(L, nH) for n, (L, dU, nH) in data.items() if n not in TETRA_SET)
    hi = min(dU / weight(L, nH) for n, (L, dU, nH) in data.items() if n in TETRA_SET)
    return lo, hi


print("=" * 80)
print(" 09j：相の選択をエネルギー比較で決める（09i 案(b) の定数 A = %.1f, B = %.1f）" % (A, B))
print("=" * 80)
print(" [1] ΔU = U_min(直交格子) − U_min(四面体格子)")
print(f"{'分子':<6}{'臂長L':>8}{'ΔU':>9}  実際の相")
for n in NAMES:
    L, dU, _ = data[n]
    print(f"{n:<6}{L:>8.2f}{dU:>+9.3f}  {'四面体' if n in TETRA_SET else '直交回帰'}")
print("   全分子で ΔU > 0：斥力と回復力だけでは常に四面体格子が低い → 格子エネルギーの項が必要")
print("   ΔU は L とともに減る（直交格子の斥力の不利は 1/L で弱まる）")
print("=" * 80)

print(" [2] 直交格子が λ L^p だけ有利とした場合に、8 分子を正しく分ける λ の範囲")
print(f"{'p':>6}{'λ 下限':>12}{'λ 上限':>12}   判定")
for p in (-1, -0.5, 0, 0.5, 1, 2):
    lo, hi = window(lambda L, nH, p=p: L ** p)
    ok = f"成立（幅 {100 * (hi - lo) / lo:.0f}%）" if hi > lo else "不成立"
    print(f"{p:>6}{lo:>12.5f}{hi:>12.5f}   {ok}")
print("   p ≥ 0 で成立、p < 0（格子エネルギー差が L とともに減る）では不成立")
print("=" * 80)

print(" [3] p = 1、格子線 1 本 = ε L として、格子線の数え方ごとの ε の範囲")
schemes = [
    ("A 中心から延びる線（4 本 / 6 本）", lambda nH: 2),
    ("B 外側の辺（テトラ 6 本 / 八面体 12 本）", lambda nH: 6),
    ("E A＋B（10 本 / 18 本）", lambda nH: 8),
    ("C H が占めていない軸（0 本 / 6−n_H 本）", lambda nH: 6 - nH),
]
for label, dn in schemes:
    lo, hi = window(lambda L, nH, dn=dn: L * dn(nH))
    ok = f"成立（幅 {100 * (hi - lo) / lo:.0f}%）" if hi > lo else "不成立"
    print(f"   {label:<34} ε = {lo:.5f} 〜 {hi:.5f}  {ok}")
lo0, hi0 = window(lambda L, nH: 6 - nH)
print(f"   （C を p = 0 で試すと {'成立' if hi0 > lo0 else '不成立'}）")
print("   A・B・E は本数の差が分子によらず一定なので ε の目盛りが変わるだけで区別できない。")
print("   C は XH2 と XH3 で差が違うため実際の検証になり、p = 1 でのみ狭い範囲で成立。")
print("=" * 80)
print(" 未確定")
print("  ・どの格子線が結合角に効くか（A・B・C のどれか）")
print("  ・格子線 1 本のエネルギーと臂 1 本の整列エネルギーの比 ε（約 0.4〜2% が必要）")
print("    ε が TTT から決まれば、相の境界が当てはめなしで予測できる")
print("=" * 80)
