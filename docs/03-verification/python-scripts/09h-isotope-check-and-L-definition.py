#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
09h：同位体チェックと臂長 L の定義の比較

問題
  09b〜09g では臂長を L = n^(1/3)、n = Z·(陽子 1500 + 電子 500) + N·2000（中性子を含む）
  としてきた。この定義では中性子数が変わると L が変わり、結合角の予測も変わる。
  実際の分子では、中心原子の同位体を替えても結合角はほとんど変わらない
  （¹³CH4 は対称性で 109.47°、¹⁵NH3・H2¹⁸O も主同位体とほぼ同じ）。

検証内容
  [1] 外部から提案された相判定 Λ 条件  n_elec / (n_nuc·L) > Λ（Λ は H2O で決定）を同位体に適用
  [2] 09e 統一式を 2 通りの L で計算
        L_n : 中性子を含む要素数の立方根（09b〜09g の定義）
        L_Z : 中性子を含まない要素数 Z·2000 の立方根
              （電子だけの Z·500 でも Z^(1/3) に比例するので、定数が変わるだけで予測は同じ）
      定数はどちらも C_LP を H2O、C_HH を H2S で決める（09e と同じ手順）。
  [3] 中心原子の同位体を替えたときの予測角の変化
"""

import math

E_E, E_P, E_N = 500, 1500, 2000
THETA_T = math.degrees(math.acos(-1 / 3))

# (Z, N) 主同位体
MAIN = {"C": (6, 6), "N": (7, 7), "O": (8, 8), "P": (15, 16), "S": (16, 16),
        "As": (33, 42), "Se": (34, 46), "Sb": (51, 70), "Te": (52, 78)}
MOLS = [("CH4", "C", 3 + 1, 0, 109.47), ("NH3", "N", 3, 1, 106.7), ("H2O", "O", 2, 2, 104.5),
        ("PH3", "P", 3, 1, 93.5), ("H2S", "S", 2, 2, 92.1), ("AsH3", "As", 3, 1, 91.8),
        ("H2Se", "Se", 2, 2, 91.0), ("SbH3", "Sb", 3, 1, 91.7), ("H2Te", "Te", 2, 2, 90.3)]
MOLS = [(f, c, (4 if f == "CH4" else nH), nLP, e) for f, c, nH, nLP, e in MOLS]
PERIOD = {"C": 2, "N": 2, "O": 2, "P": 3, "S": 3, "As": 4, "Se": 4, "Sb": 5, "Te": 5}
# 同位体（中心原子）
ISOTOPES = [("¹³CH4", "C", 6, 7, 4, 0, "CH4"), ("¹⁵NH3", "N", 7, 8, 3, 1, "NH3"),
            ("H2¹⁸O", "O", 8, 10, 2, 2, "H2O"), ("H2³⁴S", "S", 16, 18, 2, 2, "H2S")]


def L_n(Z, N):
    return (Z * (E_P + E_E) + N * E_N) ** (1 / 3)


def L_Z(Z, N):
    return (Z * (E_P + E_E)) ** (1 / 3)


def n_eff(nH, th):
    if nH != 3:
        return 1.0
    c = math.cos(math.radians(th))
    return 1 + 2 * c / (1 + c)


def make_model(Lfun):
    # 相の境界：周期2の最大 L（O）と周期3の最小 L（P）の幾何平均（範囲内ならどこでも同じ）
    Lc = math.sqrt(Lfun(*MAIN["O"]) * Lfun(*MAIN["P"]))
    C_LP = (THETA_T - 104.5) * Lfun(*MAIN["O"]) ** 2 / 2
    C_HH = (92.1 - 90) * Lfun(*MAIN["S"]) ** 2 + 2 * C_LP

    def pred(Z, N, nH, nLP):
        L = Lfun(Z, N)
        if L < Lc:
            return THETA_T - C_LP * nLP / L ** 2
        th = 90.0
        for _ in range(500):
            th = 90 + (C_HH * n_eff(nH, th) - C_LP * nLP) / L ** 2
        return th
    return pred, Lc, C_LP, C_HH


mae = lambda xs: sum(abs(x) for x in xs) / len(xs)

print("=" * 82)
print(" 09h：同位体チェックと臂長 L の定義")
print("=" * 82)

# ---- [1] Λ 条件 ----
def K(Z, N):
    ne, nn = Z * E_E, Z * E_P + N * E_N
    return ne / (nn * L_n(Z, N))

Lam = K(*MAIN["O"])
print(" [1] 相判定 Λ 条件  n_elec/(n_nuc·L) > Λ ⇒ 四面体相（Λ = H2O の値 = %.5f）" % Lam)
print(f"{'分子':<8}{'Z':>4}{'N':>4}{'K':>10}{'判定':>10}   実際")
rows = [("CH4", "C", 6, 6), ("¹³CH4", "C", 6, 7), ("NH3", "N", 7, 7), ("¹⁵NH3", "N", 7, 8),
        ("H2O", "O", 8, 8), ("H2¹⁸O", "O", 8, 10), ("PH3", "P", 15, 16)]
for name, el, Z, N in rows:
    k = K(Z, N)
    j = "四面体" if k >= Lam * (1 - 1e-12) else "直交回帰"
    actual = "直交回帰" if PERIOD[el] >= 3 else "四面体"
    mark = "" if j == actual else "  ← 誤判定"
    print(f"{name:<8}{Z:>4}{N:>4}{k:>10.5f}{j:>10}   {actual}{mark}")
print("   中性子 1〜2 個の追加で周期2の分子が直交回帰相（約90°）と判定される → 条件として不成立")
print("   また Λ を H2O で決めているため、周期2では L < L_O と同値（L_CRITICAL = L_O は循環）")
print("=" * 82)

# ---- [2] 09e を 2 通りの L で ----
models = {"L_n（中性子込み）": make_model(L_n), "L_Z（中性子なし）": make_model(L_Z)}
print(" [2] 09e 統一式：L の定義による違い（C_LP は H2O、C_HH は H2S で決定）")
for k_, (p, Lc, clp, chh) in models.items():
    print(f"   {k_}: L_CRITICAL = {Lc:.2f}, C_LP = {clp:.1f}, C_HH = {chh:.1f}")
names = list(models)
print(f"{'分子':<6}{'実験値':>9}" + "".join(f"{n[:3]+'予測':>10}{'誤差':>8}" for n in names) + "  区分")
errs = {n: [] for n in names}
for f, c, nH, nLP, e in MOLS:
    Z, N = MAIN[c]
    line = f"{f:<6}{e:>8.2f}°"
    for n in names:
        th = models[n][0](Z, N, nH, nLP)
        line += f"{th:>9.2f}°{th - e:>+8.2f}"
        if f not in ("H2O", "H2S", "CH4"):
            errs[n].append(th - e)
    kind = "基準" if f in ("H2O", "H2S") else ("当然" if f == "CH4" else "検証")
    print(line + "  " + kind)
for n in names:
    print(f"   {n} 検証 6 点 MAE = {mae(errs[n]):.2f}°  最大 {max(errs[n], key=abs):+.2f}°")
print("=" * 82)

# ---- [3] 同位体 ----
print(" [3] 中心原子の同位体を替えたときの予測角の変化（主同位体との差）")
print(f"{'分子':<8}" + "".join(f"{n[:3]+'予測':>10}{'変化':>8}" for n in names) + "   実際")
for name, el, Z, N, nH, nLP, base in ISOTOPES:
    line = f"{name:<8}"
    for n in names:
        p = models[n][0]
        th = p(Z, N, nH, nLP)
        th0 = p(*MAIN[el], nH, nLP)
        line += f"{th:>9.2f}°{th - th0:>+8.2f}"
    print(line + "   主同位体とほぼ同じ")
print("=" * 82)
print(" まとめ")
print("  ・L_n では同位体だけで NH3 が +0.12°、H2O が +0.38° 動く。L_Z では 0（定義上）。")
print("  ・L_Z にしても 検証 6 点 MAE は 0.24° → 0.28° とほぼ同じ。")
print("  ・したがって結合角を決める臂長には、中性子を含まない L_Z の方が整合的。")
print("    中性子の要素数（2000）は、結合角の幾何には寄与しないとみなすことになる。")
print("=" * 82)
