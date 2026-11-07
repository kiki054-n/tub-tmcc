#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
09k：双錐の均衡 ― 中心点と頂点の長さの釣り合いから LP の臂の長さを決める

09b〜09j では全ての臂（H・LP）を同じ長さ L に置いていた。
双錐の均衡状態では格子線は同じ長さにならない（本人の指摘）。
ここでは「頂点の重心が中心点に一致する」を均衡条件として LP の臂の長さ ρL を決める：
    L Σ_H e_i + ρ L Σ_LP w·l_k = 0   →   ρ = −(Σe·Σwl) / |Σwl|²
  四面体格子 : 4 本等長で自動的に均衡（ρ = 1）
  直交 XH2   : LP は 09i の φ = 61.87°（LP–X–H = 109.47°）→ ρ = 1.5
  直交 XH3   : 案(a) LP を C3 軸に 1 個 → ρ = √3、案(b) 3 か所に 1/3 ずつ → ρ = 1.937
エネルギー（09i と同じ形、LP–H 斥力は実距離で評価）
    U = −L Σ cos δ_i + A Σ_HH 1/d_HH + B Σ_H,LP w/d_HL,   d_HL = |L h − ρL l|
定数 A, B は H2O・H2S で決定（σ なし）。
後半で 09j と同じ方法で相の選択（四面体格子と直交格子のエネルギー比較）を調べる。
※ 「頂点の重心＝中心点」は均衡の一つの読み方（解釈）。LP の臂の長さが二つの格子で
   違うため、エネルギー比較の基準の取り方に結果が依存する。
"""

import importlib.util, io, contextlib, math, os
import numpy as np
from scipy.optimize import minimize, fsolve

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("m09i", os.path.join(HERE, "09i-lp-tetra-angle-rule-fbi-frame.py"))
m = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(m)

NAMES = ["NH3", "H2O", "PH3", "H2S", "AsH3", "H2Se", "SbH3", "H2Te"]
TETRA_SET = {"NH3", "H2O"}
PRED_09I_B = {"NH3": 0.17, "PH3": 0.05, "AsH3": 0.50, "H2Se": 0.30, "SbH3": 0.08, "H2Te": 0.68}


def rho(nH, ortho, var):
    if not ortho:
        return 1.0
    axes, lps, w = m.geometry(nH, ortho, var)
    s, t = sum(axes), w * sum(lps)
    return -float(np.dot(s, t)) / float(np.dot(t, t))


def solve(name, A, B, var, ortho=None):
    _, el, nH, _ = m.MOL[name]
    L = m.arm(el)
    if ortho is None:
        ortho = L >= m.L_CRITICAL
    axes, lps, w = m.geometry(nH, ortho, var)
    r = rho(nH, ortho, var) * L

    def U(x):
        H = [m.unit(x[3*i:3*i+3]) for i in range(nH)]
        u = -L * sum(float(np.dot(H[i], axes[i])) for i in range(nH))
        u += A * sum(1 / (2 * L * math.sin(m.gam(H[i], H[j]) / 2)) for i in range(nH) for j in range(i + 1, nH))
        u += B * w * sum(1 / np.linalg.norm(L * h - r * lp) for h in H for lp in lps)
        return u

    res = minimize(U, np.concatenate(axes), method="BFGS", options={"gtol": 1e-11, "maxiter": 20000})
    H = [m.unit(res.x[3*i:3*i+3]) for i in range(nH)]
    ang = float(np.mean([math.degrees(m.gam(H[i], H[j])) for i in range(nH) for j in range(i + 1, nH)]))
    return ang, res.fun, L


mae = lambda xs: sum(abs(x) for x in xs) / len(xs)
print("=" * 84)
print(" 09k：頂点の重心＝中心点の均衡から LP の臂の長さを決める（σ なし・定数 2 個）")
print("=" * 84)
for var in ("a", "b"):
    print(f" 案({var})  ρ(XH2) = {rho(2, True, var):.4f},  ρ(XH3) = {rho(3, True, var):.4f},  ρ(四面体) = 1")
print("-" * 84)
fits, preds = {}, {}
for var in ("a", "b"):
    A, B = fsolve(lambda p: [solve(n, p[0], p[1], var)[0] - m.MOL[n][3] for n in m.REF], [278.0, 545.0], xtol=1e-10)
    fits[var] = (A, B)
    preds[var] = {n: solve(n, A, B, var)[0] for n in m.TEST}
    print(f" 案({var})  A = {A:.1f}, B = {B:.1f}")
print("-" * 84)
print(f"{'分子':<6}{'実験値':>9}{'09i(b)誤差':>12}{'均衡(a)':>10}{'誤差':>8}{'均衡(b)':>10}{'誤差':>8}")
for n in m.TEST:
    e = m.MOL[n][3]
    pa, pb = preds["a"][n], preds["b"][n]
    print(f"{n:<6}{e:>8.2f}°{PRED_09I_B[n]:>+11.2f}{pa:>9.2f}°{pa - e:>+8.2f}{pb:>9.2f}°{pb - e:>+8.2f}")
print("-" * 84)
print(f" 検証 6 点 MAE：09i(b) {mae(PRED_09I_B.values()):.2f}° ／ 均衡(a) "
      f"{mae([preds['a'][n] - m.MOL[n][3] for n in m.TEST]):.2f}° ／ 均衡(b) "
      f"{mae([preds['b'][n] - m.MOL[n][3] for n in m.TEST]):.2f}°   （09e：定数2個＋σ で 0.24°）")
print("=" * 84)

print(" 相の選択：ΔU = U_min(直交格子) − U_min(四面体格子)")
for var in ("a", "b"):
    A, B = fits[var]
    rows = []
    for n in NAMES:
        _, ut, L = solve(n, A, B, var, False)
        _, uo, _ = solve(n, A, B, var, True)
        rows.append((n, L, uo - ut))
    print(f" 案({var}) " + "  ".join(f"{n} {d:+.2f}" for n, L, d in rows))
    print("   全分子で ΔU < 0：長さの均衡を入れると直交格子が常に低い → 四面体を保つ格子の項 λL^p が必要")
    for p in (-1, -1.5, -2, -3):
        t = max(-d / L ** p for n, L, d in rows if n in TETRA_SET)
        o = min(-d / L ** p for n, L, d in rows if n not in TETRA_SET)
        ok = f"成立（幅 {100 * (o - t) / t:.0f}%）" if o > t else "不成立"
        print(f"     p = {p:>5}: λ の範囲 {t:12.1f} 〜 {o:12.1f}  {ok}")
print("   四面体を保つ項は 1/L² 程度かそれより急に弱まる必要がある")
print("   （1/L² は関係項 F_R ∝ 1/r³ のポテンシャル ∝ 1/r² と同じ形。ただし要素数 n を掛けない形）")
print("=" * 84)
print(" 未確定")
print("  ・「頂点の重心＝中心点」は均衡条件の一つの解釈（F の大きさの釣り合い等、別の読み方もあり得る）")
print("  ・LP の臂の長さが二つの格子で違うため、エネルギー比較の基準の取り方に p の条件が依存する")
print("  ・09j の格子線の数え方（等長前提）は見直しが必要")
print("=" * 84)
