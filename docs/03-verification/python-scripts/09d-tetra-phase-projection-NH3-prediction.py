#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
09d：案B（球面射影 ＋ LP 閉じ項）を四面体相へ拡張し、NH3 を当てはめに使わずに予測する

考え方（09c 案B と同じ構造）
  正四面体（109.47°）では 4 本の臂（BP・LP）の押し合いが釣り合っている。
  BP を LP に置き換えると、その LP は余分の強さで H を押す。
  H–X–H 角 θ を閉じる量は、各 LP が H1・H2 をそれぞれ押す力の
  「角を開く方向」への射影の和（符号反転）で決まる：

      P = − Σ_LP Σ_{i=1,2}  F(LP→H_i) · u_i
          F(LP→H_i) : H_i の球面接線方向で LP から離れる向きの単位ベクトル
          u_i       : H_i の球面接線方向で相手の H から離れる向き（角を開く向き）

      θ = 109.4712° − K · P / L²

  P は正四面体の幾何から計算する（手で与えない）。
  直交回帰相（09c 案B）と同じく、変位は 1/L² でスケールし、剛性の分子依存は入れない。
  定数は K のみ（H2O で決定）。NH3 は検証点。CH4 は P=0 で 109.47°（構造上の当然）。

09b 四面体相との違い
  09b：δθ ∝ n_LP / (n_pairs + 1)   → H2O : NH3 = 1 : 1/4
  09d：δθ ∝ P / L²                  → H2O : NH3 = 2 : 1（射影から）× L² 補正
"""

import math
import numpy as np

E_ELECTRON, E_PROTON, E_NEUTRON = 500, 1500, 2000
THETA_TETRA = math.degrees(math.acos(-1.0 / 3.0))

ELEMENTS = {"C": (6, 6), "N": (7, 7), "O": (8, 8)}
# formula, central, n_H, n_LP, exp（出典は 09b）
MOLECULES = [
    ("CH4", "C", 4, 0, 109.47),
    ("NH3", "N", 3, 1, 106.7),
    ("H2O", "O", 2, 2, 104.5),
]
REF = "H2O"


def arm(s):
    Z, N = ELEMENTS[s]
    return (Z * E_PROTON + N * E_NEUTRON + Z * E_ELECTRON) ** (1.0 / 3.0)


# 正四面体の 4 頂点方向
V = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], float)
V /= np.linalg.norm(V, axis=1)[:, None]


def tangent(p, q):
    """点 p（単位ベクトル）で、q へ向かう球面接線方向の単位ベクトル"""
    t = q - np.dot(p, q) * p
    return t / np.linalg.norm(t)


def closing_projection(n_H, n_LP):
    """H0–X–H1 角を閉じる射影量 P（頂点 0,1 を H、残りから n_LP 個を LP とする）"""
    if n_LP == 0:
        return 0.0
    H = [V[0], V[1]]
    LP = [V[k] for k in range(4 - n_LP, 4)]
    P = 0.0
    for i in (0, 1):
        u = -tangent(H[i], H[1 - i])          # 角を開く向き
        for lp in LP:
            F = -tangent(H[i], lp)            # LP から離れる向き
            P -= np.dot(F, u)
    return P


# ---- 定数 K を H2O で決定 ----
mol = {m[0]: m for m in MOLECULES}
_, c, nH, nLP, e = mol[REF]
K = (THETA_TETRA - e) * arm(c) ** 2 / closing_projection(nH, nLP)

# 09b 四面体相（比較用）
C_TETRA_09B = (THETA_TETRA - e) / (nLP / (nH * (nH - 1) / 2 + 1))

print("=" * 78)
print(" 09d：四面体相への射影モデル拡張（定数 K 1個、H2O で決定）")
print("=" * 78)
print(f" K = {K:.1f}   （09c 案B の直交回帰相：C_LP = 3459.3, C_HH = 10278.6）")
print("-" * 78)
print(f"{'分子':<6}{'臂長L':>7}{'射影P':>8}{'09b予測':>10}{'誤差':>8}"
      f"{'09d予測':>10}{'誤差':>8}{'実験値':>9}  区分")
print("-" * 78)
for f, c, nH, nLP, e in MOLECULES:
    L = arm(c)
    P = closing_projection(nH, nLP)
    p09d = THETA_TETRA - K * P / L ** 2
    p09b = THETA_TETRA - C_TETRA_09B * nLP / (nH * (nH - 1) / 2 + 1)
    kind = "基準" if f == REF else ("当然（P=0）" if P == 0 else "検証")
    print(f"{f:<6}{L:>7.2f}{P:>8.3f}{p09b:>9.2f}°{p09b - e:>+8.2f}"
          f"{p09d:>9.2f}°{p09d - e:>+8.2f}{e:>8.2f}°  {kind}")
print("=" * 78)

# ---- 感度：L² 補正を外した場合、実験値の揺れ ----
_, c, nH, nLP, e = mol["NH3"]
P = closing_projection(nH, nLP)
no_L = THETA_TETRA - (THETA_TETRA - mol[REF][4]) * P / closing_projection(2, 2)
print(" 感度確認")
print(f"  NH3（L² 補正なし、射影比のみ）: {no_L:.2f}°（誤差 {no_L - e:+.2f}°）")
for e_h2o in (104.45, 104.5, 104.52):
    Kx = (THETA_TETRA - e_h2o) * arm("O") ** 2 / closing_projection(2, 2)
    th = THETA_TETRA - Kx * P / arm("N") ** 2
    print(f"  H2O 基準を {e_h2o:.2f}° にした場合の NH3: {th:.2f}°（106.7° との差 {th - e:+.2f}°）")
print("=" * 78)
print(" 注意")
print("  ・射影 P は理想正四面体の幾何で評価した一次近似。角が縮んだ後の")
print("    H–H 隣接ペアの射影変化（09c 案B の n_eff に相当）は未導入。")
print("  ・四面体相の K と直交回帰相の C_LP は別々に決めており、相互関係は未確定。")
print("  ・検証点は NH3 の 1 点のみ。")
print("=" * 78)
