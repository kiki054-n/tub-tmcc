#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
09c：3ペア系（XH3）の実効ペア数 補正の比較
     09b（正直版・定数2個）の未解決課題に対する2つの補正案を並べて検証する。

09b の残差
  直交回帰相の3ペア系は式の C(3,2)=3 ペアでは開きすぎる。
  実験値を再現するには実効ペア数が 3 の 0.50〜0.65 倍（PH3 1.63, AsH3 1.51, SbH3 1.96）。

案A：暫定的な経験的補正（理論的根拠なし）
  n_eff / n_pairs = 1 − γ·(n_H − 2)·(L_ref / L) + β·(L − L_ref)/L_ref
  γ, β は PH3・AsH3・SbH3 の3点に最小二乗で当てはめる → この3点の誤差は小さくて当然。
  当てはめに使わない検証点は NH3（四面体相）のみ。

案B：幾何学的補正（球面上の射影）＋孤立電子対（LP）の閉じ項
  (1) 各 H–X–H 角を開くのは、その角をはさむ H–H ペアの斥力が主である。
      XH3 の他の2ペアの斥力は、球面三角形の頂角 A の余弦で射影される：
          cos A = cos θ / (1 + cos θ)        （3つの角が等しい三角錐）
      → n_eff = 1 + 2 cos A   （XH2 では n_eff = 1）
      θ ≈ 90° では cos A ≈ 0 なので、n_eff ≈ 1（式の 3 ではない）。
  (2) LP は結合対を閉じる向きに押す：  − C_LP · n_LP / L²
      XH2 は LP 2個、XH3 は LP 1個 → XH3 の方が閉じる力が弱い。
  θ = 90° + (C_HH · n_eff(θ) − C_LP · n_LP) / L²  を自己無撞着に解く。
  C_HH, C_LP は H2S・PH3 の2点で決める。検証点は AsH3・H2Se・SbH3・H2Te の4点。
  ※ 案Bは直交回帰相のみ。四面体相（NH3）には未適用。

出力は実行したそのままを 09c-hydride-3pair-corrections.output.txt に記録。
"""

import math
import numpy as np

E_ELECTRON, E_PROTON, E_NEUTRON = 500, 1500, 2000
THETA_TETRA = math.degrees(math.acos(-1.0 / 3.0))
THETA_ORTHO = 90.0

ELEMENTS = {
    "C": (6, 6), "N": (7, 7), "O": (8, 8),
    "P": (15, 16), "S": (16, 16),
    "As": (33, 42), "Se": (34, 46),
    "Sb": (51, 70), "Te": (52, 78),
}

# formula, central, n_H, n_LP, exp_angle（出典は 09b を参照）
MOLECULES = [
    ("CH4", "C", 4, 0, 109.47),
    ("NH3", "N", 3, 1, 106.7),
    ("H2O", "O", 2, 2, 104.5),
    ("PH3", "P", 3, 1, 93.5),
    ("H2S", "S", 2, 2, 92.1),
    ("AsH3", "As", 3, 1, 91.8),
    ("H2Se", "Se", 2, 2, 91.0),
    ("SbH3", "Sb", 3, 1, 91.7),
    ("H2Te", "Te", 2, 2, 90.3),
]
MOL = {m[0]: m for m in MOLECULES}
L_CRITICAL = 34.0


def n_elements(s):
    Z, N = ELEMENTS[s]
    return Z * E_PROTON + N * E_NEUTRON + Z * E_ELECTRON


def arm(s):
    return n_elements(s) ** (1.0 / 3.0)


def n_pairs(nH):
    return nH * (nH - 1) / 2.0


L_REF = arm("O")

# ---- 09b の定数（H2O, H2S で決定） ----
C_TETRA = (THETA_TETRA - MOL["H2O"][4]) / (MOL["H2O"][3] / (n_pairs(2) + 1))
C_HH_09B = (MOL["H2S"][4] - THETA_ORTHO) * arm("S") ** 2 / n_pairs(2)


def pred_09b(m):
    f, c, nH, nLP, e = m
    L = arm(c)
    if L < L_CRITICAL:
        return THETA_TETRA - C_TETRA * nLP / (n_pairs(nH) + 1)
    return THETA_ORTHO + C_HH_09B * n_pairs(nH) / L ** 2


# ---- 案A：経験的補正 ----
def ratio_A(nH, L, g, b):
    return 1 - g * (nH - 2) * (L_REF / L) + b * (L - L_REF) / L_REF


fitA = ["PH3", "AsH3", "SbH3"]
rows, y = [], []
for f in fitA:
    _, c, nH, nLP, e = MOL[f]
    L = arm(c)
    need = (e - THETA_ORTHO) * L ** 2 / C_HH_09B / n_pairs(nH)
    rows.append([-(nH - 2) * L_REF / L, (L - L_REF) / L_REF])
    y.append(need - 1)
GAMMA, BETA = np.linalg.lstsq(np.array(rows), np.array(y), rcond=None)[0]


def pred_A(m):
    f, c, nH, nLP, e = m
    L = arm(c)
    neff = n_pairs(nH) * ratio_A(nH, L, GAMMA, BETA) if nH == 3 else n_pairs(nH)
    if L < L_CRITICAL:
        return THETA_TETRA - C_TETRA * nLP / (neff + 1)
    return THETA_ORTHO + C_HH_09B * neff / L ** 2


# ---- 案B：幾何学的射影 ＋ LP 閉じ項 ----
def neff_B(nH, theta):
    if nH == 2:
        return 1.0
    c = math.cos(math.radians(theta))
    return 1 + 2 * c / (1 + c)


fitB = ["H2S", "PH3"]
A = [[neff_B(MOL[f][2], MOL[f][4]), -MOL[f][3]] for f in fitB]
bvec = [(MOL[f][4] - THETA_ORTHO) * arm(MOL[f][1]) ** 2 for f in fitB]
C_HH_B, C_LP_B = np.linalg.solve(np.array(A), np.array(bvec))


def pred_B(m):
    f, c, nH, nLP, e = m
    L = arm(c)
    if L < L_CRITICAL:
        return None
    th = THETA_ORTHO
    for _ in range(500):
        th = THETA_ORTHO + (C_HH_B * neff_B(nH, th) - C_LP_B * nLP) / L ** 2
    return th


# ---- 出力 ----
def tag(f, fit):
    return "基準" if f in fit else "検証"


print("=" * 92)
print(" 09c：3ペア系の補正比較（09b／案A 経験的／案B 幾何学的）")
print("=" * 92)
print(f" 09b : C_TETRA = {C_TETRA:.4f}（H2O）, C_HH = {C_HH_09B:.1f}（H2S）")
print(f" 案A : γ = {GAMMA:.4f}, β = {BETA:+.4f}（PH3・AsH3・SbH3 に最小二乗）")
print(f"       ※ β < 0：臂長が伸びると実効ペア数はさらに減る向き")
print(f" 案B : C_HH = {C_HH_B:.1f}, C_LP = {C_LP_B:.1f}（H2S・PH3 で決定）")
print("-" * 92)
print(f"{'分子':<6}{'実験値':>8} | {'09b':>8}{'誤差':>7} | {'案A':>8}{'誤差':>7}{'':>5}"
      f" | {'案B':>8}{'誤差':>7}{'':>5}{'n_eff':>7}")
print("-" * 92)

refs_09b = ["H2O", "H2S"]
fitA_all = refs_09b + fitA
fitB_all = fitB
eb, ea, ebB = [], [], []
for m in MOLECULES:
    f, c, nH, nLP, e = m
    p0, pa, pb = pred_09b(m), pred_A(m), pred_B(m)
    line = f"{f:<6}{e:>7.2f}° | {p0:>7.2f}°{p0 - e:>+7.2f}"
    line += f" | {pa:>7.2f}°{pa - e:>+7.2f} {tag(f, fitA_all):<4}"
    if f not in refs_09b:
        eb.append((f, p0 - e))
    if f not in fitA_all:
        ea.append((f, pa - e))
    if pb is None:
        line += f" | {'—':>8}{'':>7}{'':>5}{'':>7}"
    else:
        line += f" | {pb:>7.2f}°{pb - e:>+7.2f} {tag(f, fitB_all):<4}{neff_B(nH, pb):>7.3f}"
        if f not in fitB_all:
            ebB.append((f, pb - e))
    print(line)

print("=" * 92)


def mae(lst):
    return sum(abs(x) for _, x in lst) / len(lst)


print(f" 09b  検証 {len(eb)} 点 MAE = {mae(eb):.2f}°  ({', '.join(f for f, _ in eb)})")
print(f" 案A  検証 {len(ea)} 点 MAE = {mae(ea):.2f}°  ({', '.join(f for f, _ in ea)})")
print(f" 案B  検証 {len(ebB)} 点 MAE = {mae(ebB):.2f}°  ({', '.join(f for f, _ in ebB)})")
same = [(f, x) for f, x in eb if f in dict(ebB)]
print(f"   同じ4点での 09b MAE = {mae(same):.2f}°  → 案B {mae(ebB):.2f}°")
print("=" * 92)
print(" 注意")
print("  ・案A は3点に2パラメータを当てはめた経験式。NH3 1点の改善のみが検証になる。")
print("  ・案B の n_eff は 0.87〜0.94 で、式の 3 ではなく「各角は自分のペアで開く」に近い。")
print("    XH2 と XH3 の差は、主に LP の数（2個／1個）の差で説明される。")
print("  ・案B は C_HH と C_LP が大きく打ち消し合う（差 ≈ 09b の C_HH）。感度の確認が必要。")
print("  ・案B は四面体相（NH3）に未適用。")
print("=" * 92)
