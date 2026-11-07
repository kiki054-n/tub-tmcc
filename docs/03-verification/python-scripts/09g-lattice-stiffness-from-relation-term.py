#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
09g：格子剛性を TTT の関係項 F_R から導き、一つのエネルギー関数で全 9 分子を解く

出発点（three-field-gravity-equation.md、第3.4節・第5節）
    F_R = −G_R (E1 + E2) / r³ · f(θ, φ) · n̂
    f(θ, φ) = cos θ_Ce1 · cos θ_Ce2 · cos(φ1 − φ2)
    n̂ は角度方向（r に垂直）

(1) F_R は角度方向の力なので、ポテンシャル U_R から  F_R = −(1/r) ∂U_R/∂θ  と書ける。
        U_R = −G_R (E_X + E_H) / r² · cos θ_Ce1        （H の Ce は結合軸にそろう：θ_Ce2 = 0）
    θ_Ce1 は中心原子 X の格子軸と X–H 結合のずれ角 δ_i。
    → F_R は X と H の 2 体の項なので、**回復力は H ごと（格子軸への引き戻し）**にかかる。
       H–X–H の「角ごと」の回復力は TTT の式からは出てこない（09f の問いへの答え）。

(2) 剛性のスケーリング
        ∂²U_R/∂δ² = G_R (E_X + E_H) / L²  ≈  G_R · E_X / L²
    質量＝格子線エネルギーの蓄積（同文書 第9節）なので E_X ∝ n_X（中心原子の要素数）。
        剛性 ∝ n / L²       ← molecular-bond-angles.md 1.3 節の「回復力 ∝ n/L²」と一致
    臂長 L = n^(1/3) なので  剛性 ∝ L。

(3) 臂どうしの斥力（1/r² の力、ポテンシャル ∝ 1/d）
        d = 2 L sin(γ/2)    γ：2 本の臂のなす角
    角度方向の力のモーメント（トルク）∝ L · 1/L² = 1/L
    → 変位 ∝ (1/L) / L = 1/L²   ← 09b〜09f で経験的に使ってきた 1/L² が導かれる。

エネルギー関数（全分子共通、変わるのは格子軸だけ）
    U = − κ L Σ_H cos δ_i  +  (A/L) Σ_{H–H} 1/(2 sin(γ/2))  +  (B/L) Σ_{H–LP} 1/(2 sin(γ/2))
    κ = 1 に規格化。A：H–H 斥力、B：LP–H 斥力。
    格子軸：四面体相は正四面体の頂点、直交回帰相は ±XYZ（6 軸）。

09e との違い
  ・σ（四面体相で H–H 項を 0 にする仮定）が不要。四面体相でも H–H 斥力をそのまま入れ、
    CH4 は対称性だけで 109.47° になる。
  ・定数 A, B は両相で共通。H2O と H2S の 2 点で決め、残り 7 点（CH4 含む）を予測する。
  ・LP の向きは固定（四面体相：正四面体の頂点、XH3：C3 軸、XH2：傾き φ を走査）。
"""

import math
import warnings
import numpy as np
from scipy.optimize import minimize, fsolve

warnings.filterwarnings("ignore", category=RuntimeWarning)

E_ELECTRON, E_PROTON, E_NEUTRON = 500, 1500, 2000
ELEMENTS = {
    "C": (6, 6), "N": (7, 7), "O": (8, 8),
    "P": (15, 16), "S": (16, 16),
    "As": (33, 42), "Se": (34, 46),
    "Sb": (51, 70), "Te": (52, 78),
}
MOLECULES = [  # formula, central, n_H, exp（出典は 09b）
    ("CH4", "C", 4, 109.47), ("NH3", "N", 3, 106.7), ("H2O", "O", 2, 104.5),
    ("PH3", "P", 3, 93.5), ("H2S", "S", 2, 92.1),
    ("AsH3", "As", 3, 91.8), ("H2Se", "Se", 2, 91.0),
    ("SbH3", "Sb", 3, 91.7), ("H2Te", "Te", 2, 90.3),
]
MOL = {m[0]: m for m in MOLECULES}
L_CRITICAL = 34.0
PRED_09E = {"CH4": 109.47, "NH3": 106.75, "H2O": 104.50, "PH3": 93.13, "H2S": 92.10,
            "AsH3": 91.88, "H2Se": 91.14, "SbH3": 91.40, "H2Te": 90.82}
REF = ("H2O", "H2S")


def arm(s):
    Z, N = ELEMENTS[s]
    return (Z * E_PROTON + N * E_NEUTRON + Z * E_ELECTRON) ** (1.0 / 3.0)


def unit(v):
    return v / np.linalg.norm(v)


def gam(a, b):
    return math.acos(max(-1.0, min(1.0, float(np.dot(a, b)))))


TETRA = np.array([unit(np.array(v, float)) for v in
                  [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]])
X, Y, Z_ = np.eye(3)


def geometry(nH, ortho, phi):
    if not ortho:
        return list(TETRA[:nH]), list(TETRA[nH:])
    if nH == 3:
        return [X, Y, Z_], [unit(-np.ones(3))]
    c, s = math.cos(math.radians(phi)), math.sin(math.radians(phi))
    r2 = math.sqrt(2)
    return [X, Y], [np.array([-c / r2, -c / r2, s]), np.array([-c / r2, -c / r2, -s])]


def solve_angle(nH, L, A, B, ortho, phi):
    axes, lps = geometry(nH, ortho, phi)

    def energy(x):
        H = [unit(x[3 * i:3 * i + 3]) for i in range(nH)]
        U = -L * sum(float(np.dot(H[i], axes[i])) for i in range(nH))
        U += (A / L) * sum(1 / (2 * math.sin(gam(H[i], H[j]) / 2))
                           for i in range(nH) for j in range(i + 1, nH))
        U += (B / L) * sum(1 / (2 * math.sin(gam(h, lp) / 2)) for h in H for lp in lps)
        return U

    res = minimize(energy, np.concatenate(axes), method="BFGS",
                   options={"gtol": 1e-11, "maxiter": 20000})
    H = [unit(res.x[3 * i:3 * i + 3]) for i in range(nH)]
    return float(np.mean([math.degrees(gam(H[i], H[j]))
                          for i in range(nH) for j in range(i + 1, nH)]))


def predict(name, A, B, phi):
    _, c, nH, _ = MOL[name]
    L = arm(c)
    return solve_angle(nH, L, A, B, L >= L_CRITICAL, phi)


def fit(phi):
    def resid(p):
        return [predict(m, p[0], p[1], phi) - MOL[m][3] for m in REF]
    # 初期値：線形近似からの概算
    return fsolve(resid, [3000.0, 3500.0], xtol=1e-10)


mae = lambda xs: sum(abs(x) for x in xs) / len(xs)
TEST = [m[0] for m in MOLECULES if m[0] not in REF]

print("=" * 86)
print(" 09g：関係項 F_R から導いた格子剛性（H ごと、∝ n/L²）による全 9 分子の統一計算")
print("    U = −L Σcos δ_i + (A/L) Σ_HH 1/(2sin(γ/2)) + (B/L) Σ_HLP 1/(2sin(γ/2))")
print("=" * 86)

phi0 = math.degrees(math.acos(math.sqrt(2.0 / 3.0)))
results = {}
for phi in (phi0, 45.0, 54.74):
    A, B = fit(phi)
    preds = {m[0]: predict(m[0], A, B, phi) for m in MOLECULES}
    results[phi] = (A, B, preds)

A, B, preds = results[phi0]
print(f" XH2 の LP 傾き φ = {phi0:.2f}°（LP–X–LP = 70.5°）")
print(f" 定数（両相共通、H2O・H2S で決定）：A = {A:.1f}, B = {B:.1f}, B/A = {B / A:.3f}")
print("-" * 86)
print(f"{'分子':<6}{'相':<8}{'臂長L':>7}{'09g':>9}{'誤差':>8}{'09e':>9}{'誤差':>8}{'実験値':>9}  区分")
print("-" * 86)
for name, c, nH, e in MOLECULES:
    L = arm(c)
    ph = "四面体" if L < L_CRITICAL else "直交回帰"
    p = preds[name]
    kind = "基準" if name in REF else "検証"
    print(f"{name:<6}{ph:<8}{L:>7.2f}{p:>8.2f}°{p - e:>+8.2f}{PRED_09E[name]:>8.2f}°"
          f"{PRED_09E[name] - e:>+8.2f}{e:>8.2f}°  {kind}")
print("-" * 86)
t7 = [preds[n] - MOL[n][3] for n in TEST]
print(f" 検証 7 点（CH4 含む）MAE：09g {mae(t7):.2f}° ／ 09e {mae([PRED_09E[n] - MOL[n][3] for n in TEST]):.2f}°")
print(f"   ※ 09e は H2O・H2S の 2 点で定数 2 個、ただし σ（四面体相の H–H 項 0）を仮定。")
print(f"     09g は同じ 2 点・定数 2 個で σ なし、両相で A・B 共通。")
print("=" * 86)
print(" 感度：XH2 の LP 傾き φ（定数はそのつど H2O・H2S で再決定）")
print(f"{'φ':>8}{'A':>9}{'B':>9}" + "".join(f"{n:>8}" for n in TEST) + f"{'MAE':>7}")
for phi, (A_, B_, pr) in results.items():
    errs = [pr[n] - MOL[n][3] for n in TEST]
    print(f"{phi:>7.2f}°{A_:>9.1f}{B_:>9.1f}" + "".join(f"{x:>+8.2f}" for x in errs)
          + f"{mae(errs):>7.2f}")
print("=" * 86)
print(" [2] LP 傾き φ も定数として決める（A, B, φ の 3 個を H2O・H2S・PH3 で決定）")
REF3 = ("H2O", "H2S", "PH3")
A3, B3, phi3 = fsolve(lambda p: [predict(m, p[0], p[1], p[2]) - MOL[m][3] for m in REF3],
                      [378.0, 665.0, 65.0], xtol=1e-10)
print(f"     A = {A3:.1f}, B = {B3:.1f}, φ = {phi3:.2f}°（LP–X–LP = {2 * phi3:.1f}°）, B/A = {B3 / A3:.3f}")
T3 = [m[0] for m in MOLECULES if m[0] not in REF3]
print(f"{'分子':<6}{'09g[2]':>9}{'誤差':>8}{'09e':>9}{'誤差':>8}{'実験値':>9}")
e3, e3e = [], []
for n in T3:
    p_ = predict(n, A3, B3, phi3)
    e = MOL[n][3]
    print(f"{n:<6}{p_:>8.2f}°{p_ - e:>+8.2f}{PRED_09E[n]:>8.2f}°{PRED_09E[n] - e:>+8.2f}{e:>8.2f}°")
    if n != "CH4":
        e3.append(p_ - e)
        e3e.append(PRED_09E[n] - e)
print(f" 検証 5 点（NH3, AsH3, H2Se, SbH3, H2Te）MAE：09g[2] {mae(e3):.2f}°（定数3個、σ なし）"
      f" ／ 09e {mae(e3e):.2f}°（定数2個 ＋ σ）")
print("   ※ 09e は PH3 を当てはめに使っていないので、PH3 を含めれば 09e の方が検証点は 1 つ多い")
print("=" * 86)
print(" 導出で入れた仮定")
print("  ・H の Ce は結合軸にそろう（θ_Ce2 = 0）、E_X + E_H ≈ E_X ∝ n_X")
print("  ・臂どうしの斥力は 1/r²（ポテンシャル ∝ 1/d）、LP も臂長 L の位置にあるとみなす")
print("  ・LP の向きは固定。相の切り替え（L_CRITICAL）は依然として任意")
print("=" * 86)
