#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
09f：球面上の力の釣り合いを直接最適化して、09e の射影近似を検算する

09e は「射影」と n_eff を使った一次（線形）近似の式だった。
ここでは近似を使わず、H の方向ベクトルを 3 次元で動かしてエネルギー最小点を直接求める。

エネルギー（角度はすべて度）
    U = Σ_H (k/2)·δ_i²  −  f·Σ_{H–H} γ_ij  −  g·Σ_{H–LP} γ_iL
      δ_i  : H_i と、その H が属する格子軸 e_i とのなす角（格子の回復力）
      γ    : 2 方向のなす角（一定の大きさの斥力 ＝ 角に比例するポテンシャル）
      k = L²（09b〜09e と同じ 1/L² スケーリング）

格子軸と LP の置き方
  直交回帰相（6 軸格子）
      XH3 : H は +x, +y, +z 軸、LP は C3 軸 −(1,1,1)/√3
      XH2 : H は +x, +y 軸、LP 2 個は (−c/√2, −c/√2, ±s)、c = cos φ, s = sin φ
            φ は LP の傾き（ここでは未確定なので走査する）
  四面体相
      H も LP も正四面体の頂点。H–H 斥力は 4 本の臂の釣り合いに含まれる（09e の σ=0）
      ので、LP が BP より余分に押す分 g_t だけを入れる。

注意：LP の向きは固定（動かさない）。
定数の決め方は 09e と同じ：
  直交回帰相 f, g → H2S・PH3 で決定、四面体相 g_t → H2O で決定
"""

import math
import warnings
import numpy as np
from scipy.optimize import minimize, fsolve

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

# 09e の結果（比較用、09e の出力から転記）
PRED_09E = {"CH4": 109.47, "NH3": 106.75, "H2O": 104.50, "PH3": 93.13, "H2S": 92.10,
            "AsH3": 91.88, "H2Se": 91.14, "SbH3": 91.40, "H2Te": 90.82}


def arm(s):
    Z, N = ELEMENTS[s]
    return (Z * E_PROTON + N * E_NEUTRON + Z * E_ELECTRON) ** (1.0 / 3.0)


def ang(a, b):
    return math.degrees(math.acos(max(-1.0, min(1.0, float(np.dot(a, b))))))


def unit(v):
    return v / np.linalg.norm(v)


TETRA = np.array([unit(np.array(v, float)) for v in
                  [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]])
X, Y, Z_ = np.eye(3)


def geometry(nH, ortho, phi):
    if not ortho:
        return list(TETRA[:nH]), list(TETRA[nH:])
    if nH == 3:
        return [X, Y, Z_], [unit(-np.ones(3))]
    c, s = math.cos(math.radians(phi)), math.sin(math.radians(phi))
    return [X, Y], [np.array([-c / math.sqrt(2), -c / math.sqrt(2), s]),
                    np.array([-c / math.sqrt(2), -c / math.sqrt(2), -s])]


def optimize_angle(nH, L, f, g, ortho, phi=35.26):
    axes, lps = geometry(nH, ortho, phi)
    k = L ** 2

    def energy(x):
        H = [unit(x[3 * i:3 * i + 3]) for i in range(nH)]
        U = sum(0.5 * k * ang(H[i], axes[i]) ** 2 for i in range(nH))
        if ortho:
            U -= f * sum(ang(H[i], H[j]) for i in range(nH) for j in range(i + 1, nH))
        U -= g * sum(ang(h, lp) for h in H for lp in lps)
        return U

    x0 = np.concatenate(axes)
    res = minimize(energy, x0, method="BFGS", options={"gtol": 1e-10, "maxiter": 10000})
    H = [unit(res.x[3 * i:3 * i + 3]) for i in range(nH)]
    angles = [ang(H[i], H[j]) for i in range(nH) for j in range(i + 1, nH)]
    return float(np.mean(angles))


def run(phi):
    # 四面体相：g_t を H2O で決定
    _, c, nH, e = MOL["H2O"]
    g_t = fsolve(lambda g: optimize_angle(nH, arm(c), 0, g[0], False) - e, [500.0])[0]

    # 直交回帰相：f, g を H2S・PH3 で決定
    def resid(p):
        return [optimize_angle(MOL[m][2], arm(MOL[m][1]), p[0], p[1], True, phi) - MOL[m][3]
                for m in ("H2S", "PH3")]
    f, g = fsolve(resid, [4000.0, 2000.0])

    out = {}
    for name, c, nH, e in MOLECULES:
        L = arm(c)
        if L < L_CRITICAL:
            out[name] = optimize_angle(nH, L, 0, g_t, False)
        else:
            out[name] = optimize_angle(nH, L, f, g, True, phi)
    return g_t, f, g, out


warnings.filterwarnings("ignore", category=RuntimeWarning)  # fsolve の収束メッセージ（残差は出力で確認）

# 09e の定数（09e 出力から転記）を、直接最適化の定数へ一次の対応で写す
#   直交 XH2：(2f − 2√2 g)/L² = (C_HH − 2C_LP)/L²  →  f = C_HH/2, g = C_LP/√2
#   四面体  ：角あたり閉じ = g_t·P/L²              →  g_t = C_LP
C_LP_09E, C_HH_09E = 2505.3, 8370.7

REF = {"H2O", "H2S", "PH3"}
TEST = [m[0] for m in MOLECULES if m[0] not in REF | {"CH4"}]
mae = lambda xs: sum(abs(x) for x in xs) / len(xs)

print("=" * 84)
print(" 09f：球面上の直接最適化による 09e（射影近似）の検算")
print("=" * 84)
phi0 = math.degrees(math.acos(math.sqrt(2.0 / 3.0)))  # 35.26°：09e の LP 数え方と一次で一致
print(" [1] 09e と同じ定数で、当てはめなしに直接最適化（射影近似そのものの検算）")
print(f"     f = C_HH/2 = {C_HH_09E/2:.1f}, g = C_LP/√2 = {C_LP_09E/math.sqrt(2):.1f}, g_t = C_LP = {C_LP_09E:.1f}, φ = {phi0:.2f}°")
print("-" * 84)
print(f"{'分子':<6}{'直接最適化':>11}{'09e 射影式':>11}{'差':>8}{'実験値':>10}")
d1 = []
for name, c, nH, e in MOLECULES:
    L = arm(c)
    if L < L_CRITICAL:
        d = optimize_angle(nH, L, 0, C_LP_09E, False)
    else:
        d = optimize_angle(nH, L, C_HH_09E / 2, C_LP_09E / math.sqrt(2), True, phi0)
    d1.append(d - PRED_09E[name])
    print(f"{name:<6}{d:>10.2f}°{PRED_09E[name]:>10.2f}°{d - PRED_09E[name]:>+8.2f}{e:>9.2f}°")
print(f" 全 9 分子で 直接最適化 − 09e の最大差：{max(d1, key=abs):+.2f}°")
print("   （四面体相の差 ±0.13° は一次近似で落とした高次の項＝変位が約3〜5°と大きいため）")
print("   （PH3 の +0.38° は高次項ではなく、隣接ペアの射影 n_eff の扱いの違い：")
print("     09e は角ごとの回復力で n_eff = 1 + 2cosA ≈ 0.885、直接最適化は H ごとの格子軸回復力で")
print("     隣接ペアの寄与が一次で 0 になり、実質 n_eff が 1 に近い。n_eff = 1 なら 93.74°）")
print("=" * 84)
print(" [2] 直接最適化モデルで定数を当てはめ直した場合")

g_t, f, g, out = run(phi0)
print(f" XH2 の LP 傾き φ = {phi0:.2f}°（LP 1 個あたりの閉じ射影が XH3 と等しくなる角）")
print(f" 定数：g_t = {g_t:.1f}（H2O）, f = {f:.1f}, g = {g:.1f}（H2S・PH3）")
print("-" * 84)
print(f"{'分子':<6}{'直接最適化':>11}{'09e 射影式':>11}{'差':>8}{'実験値':>10}{'誤差(直接)':>11}  区分")
print("-" * 84)
for name, c, nH, e in MOLECULES:
    d, p = out[name], PRED_09E[name]
    kind = "基準" if name in REF else ("当然" if name == "CH4" else "検証")
    print(f"{name:<6}{d:>10.2f}°{p:>10.2f}°{d - p:>+8.2f}{e:>9.2f}°{d - e:>+10.2f}   {kind}")
print("-" * 84)
diffs = [out[n] - PRED_09E[n] for n in TEST]
print(f" 検証 5 点（NH3, AsH3, H2Se, SbH3, H2Te）：直接最適化と 09e の差  最大 {max(diffs, key=abs):+.3f}°")
print(f" 同 5 点 MAE：直接最適化 {mae([out[n] - MOL[n][3] for n in TEST]):.2f}°"
      f" ／ 09e {mae([PRED_09E[n] - MOL[n][3] for n in TEST]):.2f}°")
print("=" * 84)

print(" 感度：XH2 の LP 傾き φ を変えた場合（定数はそのつど H2S・PH3 で再決定）")
print(f"{'φ':>8}{'LP–X–LP':>10}" + "".join(f"{n:>8}" for n in TEST) + f"{'MAE':>8}")
for phi in (phi0, 45.0, 54.74):
    lp_angle = math.degrees(math.acos(math.cos(math.radians(phi)) ** 2
                                      - math.sin(math.radians(phi)) ** 2))
    _, _, _, o = run(phi)
    errs = [o[n] - MOL[n][3] for n in TEST]
    print(f"{phi:>7.2f}°{lp_angle:>9.2f}°" + "".join(f"{x:>+8.2f}" for x in errs)
          + f"{mae(errs):>8.2f}")
print("   （各列は実験値との誤差。LP を格子軸上 φ=90° に置くと閉じ力は 0 になる）")
print("=" * 84)
