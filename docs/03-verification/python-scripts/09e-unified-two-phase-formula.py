#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
09e：四面体相と直交回帰相を一つの式で表す試み（共通 LP 定数の検証）

統一式
    θ = θ_base + [ σ · C_HH · n_eff(θ) − C_LP · n_LP ] / L²

    θ_base = 109.4712°, σ = 0   （四面体相：L < L_CRITICAL）
    θ_base =  90°     , σ = 1   （直交回帰相）

    n_eff(θ) = 1 + 2 cos A,  cos A = cos θ / (1 + cos θ)   （XH3、09c 案B の球面射影）
    n_eff    = 1                                        （XH2）

各項の意味
  ・LP 閉じ項 −C_LP·n_LP/L²：両相で同じ形・同じ定数。
      四面体相では 09d の射影 P が n_LP と一致する（H2O: P=2, NH3: P=1）。
  ・H–H 開き項：四面体相では 4 本の臂の釣り合い（109.47°）に含まれているので σ=0。
      直交回帰相では 90° を保つのは格子であり、臂どうしの斥力ではないので
      H–H 斥力が正味の開きとして現れる（σ=1）。

仮説：C_LP は相によらず一つ（09c 案B の 3459 と 09d の K=2505 を一本化する）
定数
  C_LP : H2O で決定（四面体相）
  C_HH : H2S で決定（直交回帰相、C_LP は固定）
  → 定数 2 個で 9 分子。検証点は CH4 を除く 6 分子（NH3, PH3, AsH3, H2Se, SbH3, H2Te）
"""

import math

E_ELECTRON, E_PROTON, E_NEUTRON = 500, 1500, 2000
THETA_TETRA = math.degrees(math.acos(-1.0 / 3.0))
THETA_ORTHO = 90.0
L_CRITICAL = 34.0

ELEMENTS = {
    "C": (6, 6), "N": (7, 7), "O": (8, 8),
    "P": (15, 16), "S": (16, 16),
    "As": (33, 42), "Se": (34, 46),
    "Sb": (51, 70), "Te": (52, 78),
}
# formula, central, n_H, n_LP, exp（出典は 09b）
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


def arm(s):
    Z, N = ELEMENTS[s]
    return (Z * E_PROTON + N * E_NEUTRON + Z * E_ELECTRON) ** (1.0 / 3.0)


def n_eff(nH, theta):
    if nH != 3:
        return 1.0
    c = math.cos(math.radians(theta))
    return 1 + 2 * c / (1 + c)


def phase(L):
    return (THETA_TETRA, 0) if L < L_CRITICAL else (THETA_ORTHO, 1)


def solve(m, C_HH, C_LP):
    f, c, nH, nLP, e = m
    L = arm(c)
    base, sigma = phase(L)
    th = base
    for _ in range(500):
        th = base + (sigma * C_HH * n_eff(nH, th) - C_LP * nLP) / L ** 2
    return th


# ---- 定数決定 ----
f, c, nH, nLP, e = MOL["H2O"]
C_LP = (THETA_TETRA - e) * arm(c) ** 2 / nLP
f, c, nH, nLP, e = MOL["H2S"]
C_HH = ((e - THETA_ORTHO) * arm(c) ** 2 + C_LP * nLP) / n_eff(nH, e)

# 比較用：09b（定数2個）
C_T09B = (THETA_TETRA - MOL["H2O"][4]) / (2 / 2)
C_H09B = (MOL["H2S"][4] - 90) * arm("S") ** 2


def pred_09b(m):
    f, c, nH, nLP, e = m
    L = arm(c)
    npair = nH * (nH - 1) / 2
    if L < L_CRITICAL:
        return THETA_TETRA - C_T09B * nLP / (npair + 1)
    return 90 + C_H09B * npair / L ** 2


print("=" * 80)
print(" 09e：二相統一式  θ = θ_base + [σ·C_HH·n_eff(θ) − C_LP·n_LP] / L²")
print("=" * 80)
print(f" C_LP = {C_LP:.1f}（H2O で決定、両相共通）")
print(f" C_HH = {C_HH:.1f}（H2S で決定）")
print(f" 参考：別々に決めた場合  09d K = 2505.3 / 09c 案B C_LP = 3459.3, C_HH = 10278.6")
print("-" * 80)
print(f"{'分子':<6}{'相':<8}{'臂長L':>7}{'n_eff':>7}{'09b':>9}{'誤差':>7}"
      f"{'09e':>9}{'誤差':>7}{'実験値':>9}  区分")
print("-" * 80)
err09b, err09e = [], []
for m in MOLECULES:
    f, c, nH, nLP, e = m
    L = arm(c)
    base, sigma = phase(L)
    t = solve(m, C_HH, C_LP)
    t0 = pred_09b(m)
    if f in ("H2O", "H2S"):
        kind = "基準"
    elif f == "CH4":
        kind = "当然（n_LP=0）"
    else:
        kind = "検証"
        err09b.append(t0 - e)
        err09e.append(t - e)
    ph = "四面体" if sigma == 0 else "直交回帰"
    ne = n_eff(nH, t) if sigma == 1 else float("nan")
    ne_s = f"{ne:7.3f}" if sigma == 1 else f"{'—':>7}"
    print(f"{f:<6}{ph:<8}{L:>7.2f}{ne_s}{t0:>8.2f}°{t0 - e:>+7.2f}"
          f"{t:>8.2f}°{t - e:>+7.2f}{e:>8.2f}°  {kind}")
print("=" * 80)
mae = lambda xs: sum(abs(x) for x in xs) / len(xs)
print(f" 検証 6 点（NH3, PH3, AsH3, H2Se, SbH3, H2Te）")
print(f"   09b（定数2個） MAE = {mae(err09b):.2f}°  最大 {max(err09b, key=abs):+.2f}°")
print(f"   09e（定数2個） MAE = {mae(err09e):.2f}°  最大 {max(err09e, key=abs):+.2f}°")
print("=" * 80)

# ---- 感度：C_LP を共通にしない場合との差 ----
print(" 感度確認：C_LP を 2505〜3459 の間で動かす（C_HH は H2S で再決定）")
for clp in (2505.3, 2800.0, 3100.0, 3459.3):
    chh = ((MOL["H2S"][4] - 90) * arm("S") ** 2 + clp * 2) / 1.0
    es = []
    for m in MOLECULES:
        if m[0] in ("CH4", "H2O", "H2S"):
            continue
        if arm(m[1]) < L_CRITICAL:
            t = THETA_TETRA - clp * m[3] / arm(m[1]) ** 2
        else:
            t = solve(m, chh, clp)
        es.append(t - m[4])
    print(f"   C_LP = {clp:7.1f}, C_HH = {chh:8.1f} → 検証6点 MAE = {mae(es):.2f}°")
print("=" * 80)
print(" 注意")
print("  ・σ（四面体相で H–H 項を 0 とする）は、基底角の切り替えと同時に入れた仮定。")
print("  ・L_CRITICAL は依然として 31.75〜39.58 の間で任意。")
print("  ・XH2 の直交回帰相予測は C_HH − 2C_LP（= 09b の C_HH）だけで決まるため、")
print("    H2Se・H2Te は 09b と同じ値になる。統一の効果が出るのは XH3 と NH3。")
print("=" * 80)
