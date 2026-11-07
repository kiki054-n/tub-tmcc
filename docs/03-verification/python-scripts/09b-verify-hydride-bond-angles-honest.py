#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TTT 2状態モデル：第2〜第5周期 水素化物の結合角 検証スクリプト（正直版）

方針
  1. 質量引力項（C_MASS_PULL × 要素数）は外す。
     旧版では第3周期以降で 5〜22° の減算になり、全体を破綻させていた。
     また H2S と PH3 で逆符号の補正が必要になるため、
     要素数に比例する単一の正の項では表現できない。
  2. 実験値は出典付きの値に改める（H2Te 89.5° → 約90°、SbH3 91.5° → 91.7°）。
  3. 定数は「どの分子で決めたか」を明示する。決めた分子の誤差 0 は当然であり、
     予測の評価には含めない。
  4. 出力は実行したそのままを記録する。

モデル
  四面体相（L < L_CRITICAL）:
      θ = 109.4712° − C_TETRA × N_LP / (N_pair + 1)
  直交回帰相（L ≥ L_CRITICAL）:
      θ = 90° + C_HH × N_pair / L²
  N_pair = C(N_H, 2)、L = (中心原子の要素数)^(1/3)
"""

import math

# ------------------------------------------------------------------
# 1. TTT 要素数
# ------------------------------------------------------------------
E_ELECTRON = 500
E_PROTON = 1500
E_NEUTRON = 2000

THETA_TETRA = math.degrees(math.acos(-1.0 / 3.0))  # 109.4712°
THETA_ORTHO = 90.0

# 中心原子 (Z, N)：最多同位体の中性子数
ELEMENTS = {
    "C": (6, 6), "N": (7, 7), "O": (8, 8),
    "P": (15, 16), "S": (16, 16),
    "As": (33, 42), "Se": (34, 46),
    "Sb": (51, 70), "Te": (52, 78),
}

# ------------------------------------------------------------------
# 2. 実験値（気相）と出典
#    ※ 個々の値は一次文献（分光測定）での再確認を推奨
# ------------------------------------------------------------------
MOLECULES = [
    # formula, central, period, n_H, n_LP, exp_angle, source
    ("CH4",  "C",  2, 4, 0, 109.47, "正四面体（定義値）"),
    ("NH3",  "N",  2, 3, 1, 106.7,  "気相分光値（教科書値 107°）"),
    ("H2O",  "O",  2, 2, 2, 104.5,  "気相分光値 104.45°"),
    ("PH3",  "P",  3, 3, 1, 93.5,   "Wikipedia: Phosphine"),
    ("H2S",  "S",  3, 2, 2, 92.1,   "AUS-e-TUTE Group 16 hydrides"),
    ("AsH3", "As", 4, 3, 1, 91.8,   "Wikipedia: Arsine"),
    ("H2Se", "Se", 4, 2, 2, 91.0,   "AUS-e-TUTE Group 16 hydrides"),
    ("SbH3", "Sb", 5, 3, 1, 91.7,   "Wikipedia: Stibine"),
    ("H2Te", "Te", 5, 2, 2, 90.3,   "教科書値 90°／分光値 約90.25°（要一次文献確認）"),
]

# 定数を決める基準分子
REF_TETRA = "H2O"
REF_ORTHO = "H2S"

# 相転移の臨界臂長：第2周期の最大 L（O: 31.75）と第3周期の最小 L（P: 39.58）の
# 間ならどこでも結果は同じ。現時点では理論から決まっていない。
L_CRITICAL = 34.0


def n_elements(symbol):
    Z, N = ELEMENTS[symbol]
    return Z * E_PROTON + N * E_NEUTRON + Z * E_ELECTRON


def arm_length(symbol):
    return n_elements(symbol) ** (1.0 / 3.0)


def n_pairs(n_H):
    return n_H * (n_H - 1) / 2.0


def tetra_factor(n_H, n_LP):
    return n_LP / (n_pairs(n_H) + 1.0)


def ortho_factor(n_H, L):
    return n_pairs(n_H) / L ** 2


# ------------------------------------------------------------------
# 3. 定数の決定（各相 1 個ずつ、計 2 個の自由パラメータ）
# ------------------------------------------------------------------
mol = {m[0]: m for m in MOLECULES}

_, c, _, nH, nLP, exp, _ = mol[REF_TETRA]
C_TETRA = (THETA_TETRA - exp) / tetra_factor(nH, nLP)

_, c, _, nH, nLP, exp, _ = mol[REF_ORTHO]
C_HH = (exp - THETA_ORTHO) / ortho_factor(nH, arm_length(c))


def predict(m):
    formula, c, period, nH, nLP, exp, src = m
    L = arm_length(c)
    if L < L_CRITICAL:
        return L, "四面体相", THETA_TETRA - C_TETRA * tetra_factor(nH, nLP)
    return L, "直交回帰相", THETA_ORTHO + C_HH * ortho_factor(nH, L)


# ------------------------------------------------------------------
# 4. 出力
# ------------------------------------------------------------------
print("=" * 86)
print(" TTT 2状態モデル 水素化物結合角 検証（質量引力項なし・定数2個）")
print("=" * 86)
print(f" C_TETRA = {C_TETRA:.4f}°  （{REF_TETRA} で決定）")
print(f" C_HH    = {C_HH:.2f}    （{REF_ORTHO} で決定）")
print(f" L_CRITICAL = {L_CRITICAL}（31.75〜39.58 の間で任意。理論的根拠は未確定）")
print("-" * 86)
print(f"{'分子':<6}{'周期':>4}{'要素数':>10}{'臂長L':>8}  {'相':<8}"
      f"{'予測角':>9}{'実験値':>9}{'誤差':>8}  区分")
print("-" * 86)

pred_errors = []
for m in MOLECULES:
    formula, c, period, nH, nLP, exp, src = m
    L, phase, theta = predict(m)
    err = theta - exp
    is_ref = formula in (REF_TETRA, REF_ORTHO)
    tag = "基準（誤差0は当然）" if is_ref else "予測"
    if not is_ref:
        pred_errors.append((formula, err))
    print(f"{formula:<6}{period:>4}{n_elements(c):>10,}{L:>8.2f}  {phase:<8}"
          f"{theta:>8.2f}°{exp:>8.2f}°{err:>+7.2f}°  {tag}")

print("=" * 86)
mae = sum(abs(e) for _, e in pred_errors) / len(pred_errors)
worst = max(pred_errors, key=lambda x: abs(x[1]))
print(f" 予測分子 {len(pred_errors)} 個の平均絶対誤差 (MAE): {mae:.2f}°")
print(f" 最大誤差: {worst[1]:+.2f}° ({worst[0]})")

two = [e for f, e in pred_errors if f in ("CH4", "H2Se", "H2Te")]
three = [e for f, e in pred_errors if f in ("NH3", "PH3", "AsH3", "SbH3")]
print(f"   2ペア以下の系（CH4・H2Se・H2Te）MAE: {sum(map(abs, two)) / len(two):.2f}°")
print(f"   3ペア系（NH3・PH3・AsH3・SbH3）  MAE: {sum(map(abs, three)) / len(three):.2f}°")
print("=" * 86)

# ------------------------------------------------------------------
# 5. 未解決課題の定量化：3ペア系に必要な「実効ペア数」
# ------------------------------------------------------------------
print(" 直交回帰相の3ペア系：実験値を再現するのに必要な実効ペア数")
print(f"   （式では C(3,2) = 3。H2S 基準の C_HH を固定した場合）")
for formula in ("PH3", "AsH3", "SbH3"):
    _, c, _, nH, nLP, exp, _ = mol[formula]
    L = arm_length(c)
    need = (exp - THETA_ORTHO) * L ** 2 / C_HH
    print(f"   {formula:<5} 必要 {need:.2f} ペア（式の 3 に対し {need / 3:.2f} 倍）")
print("=" * 86)
print(" 出典")
for m in MOLECULES:
    print(f"   {m[0]:<5} {m[5]:>7.2f}°  {m[6]}")
print("=" * 86)
