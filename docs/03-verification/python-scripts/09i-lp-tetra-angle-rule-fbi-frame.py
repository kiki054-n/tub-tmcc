#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
09i：LP の位置を「LP–X–H = 109.47° を保つ」規則で決める（FBI 枠の仮説）

背景（掛け算モデル・3階テンソル）
    E = (xX + yY + zZ) ⊗ (ω_x, ω_y, ω_z) ⊗ (f R + b S + i T) + …
    R, S, T は各臂の局所的な F・B・I 軸（R = F：臂に沿う直線、S = B・T = I：臂に垂直な曲線）。
仮説
    H はへその緒で中心原子とつながり、XYZ 格子軸 (xX+yY+zZ) に従う（直交回帰相で 90° 付近）。
    LP はへその緒を持たず、F 線（R 軸）はテトラ構造の角度を保つ：LP–X–H = 109.47°。

エネルギー関数は 09g と同じ（σ なし、両相共通の A, B、格子剛性 ∝ n/L²）：
    U = −L Σ cos δ_i + (A/L) Σ_HH 1/(2 sin(γ/2)) + (B/L) Σ_H,LP w_LP/(2 sin(γ/2))

LP の置き方
  四面体相 : 正四面体の頂点（LP–X–H = 109.47° は自動的に成立）
  直交 XH2 : (−c/√2, −c/√2, ±s)。LP·x = LP·y = −1/3 から c/√2 = 1/3
             → φ = 61.87°、LP–X–LP = 123.75°（当てはめなし）
  直交 XH3 : 案(a) C3 軸 −(1,1,1)/√3 に 1 個（LP–X–H = 125.26°、規則は満たせない）
             案(b) 2 本の H と 109.47° をなす (−1/3, −1/3, −√7/3) とその C3 回転の 3 か所に
                   1/3 ずつ（LP が ω で 3 か所を巡り、時間平均で C3 対称）
定数 A, B は H2O・H2S で決定。検証 6 点：NH3, PH3, AsH3, H2Se, SbH3, H2Te
"""

import math
import warnings
import numpy as np
from scipy.optimize import minimize, fsolve

warnings.filterwarnings("ignore", category=RuntimeWarning)

E_E, E_P, E_N = 500, 1500, 2000
ELEMENTS = {"C": (6, 6), "N": (7, 7), "O": (8, 8), "P": (15, 16), "S": (16, 16),
            "As": (33, 42), "Se": (34, 46), "Sb": (51, 70), "Te": (52, 78)}
MOLECULES = [("CH4", "C", 4, 109.47), ("NH3", "N", 3, 106.7), ("H2O", "O", 2, 104.5),
             ("PH3", "P", 3, 93.5), ("H2S", "S", 2, 92.1), ("AsH3", "As", 3, 91.8),
             ("H2Se", "Se", 2, 91.0), ("SbH3", "Sb", 3, 91.7), ("H2Te", "Te", 2, 90.3)]
MOL = {m[0]: m for m in MOLECULES}
REF = ("H2O", "H2S")
TEST = ["NH3", "PH3", "AsH3", "H2Se", "SbH3", "H2Te"]
L_CRITICAL = 34.0
PRED_09E = {"NH3": 106.75, "PH3": 93.13, "AsH3": 91.88, "H2Se": 91.14, "SbH3": 91.40, "H2Te": 90.82}


def arm(s):
    Z, N = ELEMENTS[s]
    return (Z * (E_P + E_E) + N * E_N) ** (1 / 3)


def unit(v):
    return v / np.linalg.norm(v)


def gam(a, b):
    return math.acos(max(-1.0, min(1.0, float(np.dot(a, b)))))


TETRA = [unit(np.array(v, float)) for v in [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]]
X, Y, Zv = np.eye(3)
PHI = math.degrees(math.acos(math.sqrt(2) / 3))           # 61.87°
c, s = math.cos(math.radians(PHI)), math.sin(math.radians(PHI))
LP_XH2 = [np.array([-c / math.sqrt(2), -c / math.sqrt(2), s]),
          np.array([-c / math.sqrt(2), -c / math.sqrt(2), -s])]
r7 = math.sqrt(7) / 3
LP_XH3 = {"a": ([unit(-np.ones(3))], 1.0),
          "b": ([np.array([-1/3, -1/3, -r7]), np.array([-1/3, -r7, -1/3]),
                 np.array([-r7, -1/3, -1/3])], 1/3)}


def geometry(nH, ortho, variant):
    if not ortho:
        return TETRA[:nH], TETRA[nH:], 1.0
    if nH == 2:
        return [X, Y], LP_XH2, 1.0
    lps, w = LP_XH3[variant]
    return [X, Y, Zv], lps, w


def solve(name, A, B, variant):
    _, el, nH, _ = MOL[name]
    L = arm(el)
    axes, lps, w = geometry(nH, L >= L_CRITICAL, variant)

    def U(x):
        H = [unit(x[3*i:3*i+3]) for i in range(nH)]
        u = -L * sum(float(np.dot(H[i], axes[i])) for i in range(nH))
        u += (A / L) * sum(1 / (2 * math.sin(gam(H[i], H[j]) / 2))
                           for i in range(nH) for j in range(i + 1, nH))
        u += (B / L) * w * sum(1 / (2 * math.sin(gam(h, lp) / 2)) for h in H for lp in lps)
        return u

    res = minimize(U, np.concatenate(axes), method="BFGS", options={"gtol": 1e-11, "maxiter": 20000})
    H = [unit(res.x[3*i:3*i+3]) for i in range(nH)]
    return float(np.mean([math.degrees(gam(H[i], H[j])) for i in range(nH) for j in range(i + 1, nH)]))


mae = lambda xs: sum(abs(x) for x in xs) / len(xs)

print("=" * 84)
print(" 09i：LP–X–H = 109.47° を保つ規則（FBI 枠の仮説）、σ なし・定数 2 個")
print("=" * 84)
print(f" XH2 の LP：φ = {PHI:.2f}°, LP–X–LP = {2 * PHI:.2f}°（規則から決定、当てはめなし）")
print(" XH3 の LP：")
for k, (lps, w) in LP_XH3.items():
    angs = sorted(round(math.degrees(gam(lp, e)), 2) for lp in lps[:1] for e in (X, Y, Zv))
    print(f"   案({k}) {len(lps)} か所 × 重み {w:.3f}、1 か所目と 3 本の H 軸のなす角 {angs}")
print("   ※ 3 本の H すべてと 109.47° をなす方向は存在しない（その条件は H–X–H = 109.47° と同値）")
print("-" * 84)

res = {}
for k in ("a", "b"):
    A, B = fsolve(lambda p: [solve(n, p[0], p[1], k) - MOL[n][3] for n in REF], [426.0, 722.0], xtol=1e-10)
    res[k] = (A, B, {n: solve(n, A, B, k) for n in [m[0] for m in MOLECULES]})
    print(f" 案({k}) A = {A:.1f}, B = {B:.1f}")
print("-" * 84)
print(f"{'分子':<6}{'実験値':>9}{'案(a)':>10}{'誤差':>8}{'案(b)':>10}{'誤差':>8}{'09e':>10}{'誤差':>8}  区分")
for name, el, nH, e in MOLECULES:
    pa, pb = res["a"][2][name], res["b"][2][name]
    pe = PRED_09E.get(name)
    kind = "基準" if name in REF else ("当然" if name == "CH4" else "検証")
    pe_s = f"{pe:>9.2f}°{pe - e:>+8.2f}" if pe else f"{'—':>10}{'':>8}"
    print(f"{name:<6}{e:>8.2f}°{pa:>9.2f}°{pa - e:>+8.2f}{pb:>9.2f}°{pb - e:>+8.2f}{pe_s}  {kind}")
print("-" * 84)
for k in ("a", "b"):
    print(f" 案({k}) 検証 6 点 MAE = {mae([res[k][2][n] - MOL[n][3] for n in TEST]):.2f}°")
print(f" 09e（定数 2 個 ＋ σ）検証 6 点 MAE = {mae([PRED_09E[n] - MOL[n][3] for n in TEST]):.2f}°")
print("=" * 84)
print(" 未確定")
print("  ・「LP の F 線は H と 109.47° を保つ」「XH3 の LP は ω で 3 か所を巡る」は仮説。")
print("    テンソル E からの導出は未了。")
print("  ・案(b) の LP は残る 1 本の H と 151.87° をなす。この角の意味は未検討。")
print("  ・最大誤差 H2Te（+0.68°）は実験値の一次文献確認が必要。")
print("=" * 84)
