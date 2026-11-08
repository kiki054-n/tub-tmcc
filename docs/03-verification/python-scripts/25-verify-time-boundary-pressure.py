"""
TUB-TMCC 25：時間軸（過去・現在・未来）・境界・誕生の段・境界の圧力34 の検証（2026-11-08）

この記録の数値はすべてこのスクリプトで出す。
  python3 25-verify-time-boundary-pressure.py

1. 三時刻の式（ベルレ法）で力と量の循環が切れること、総和（エネルギー）が保たれること
2. 誕生の段 T1〜T5：N_k = 2·3^(k-1), V_k = k^3, 隙間 G_k = V_k − N_k
3. 3^3+4^3+5^3 = 6^3 が成り立つのは T6 だけであること
4. 今の時刻 n、境界の圧力 P = 34/(6L^2) と宇宙定数の観測値の比較
5. 記号の関係（F±, B±, I±, π）と外積 B×I = F
6. 1/L^2 で薄まる圧力と p = −ρ（FBI = −1）を両立させるのに必要な出入り
"""
import itertools
import math

import numpy as np

print("=" * 72)
print("1. 三時刻の式（過去・現在・未来）")
print("=" * 72)


def two_body(scheme, steps=20000, dt=0.01):
    """2体の引力（G=1, m=1）。円軌道から始め、エネルギーのずれを返す。
    scheme='euler'   : 現在の力で現在の速度を更新（一時刻の式）
    scheme='verlet'  : x(n+1) = 2x(n) − x(n−1) + a(n) dt^2（三時刻の式）
    """
    m = 1.0
    r0 = 1.0
    v0 = math.sqrt(m * 2 / (4 * r0)) * 1.0  # 相対距離 2r0、換算質量で円軌道
    x = np.array([[r0, 0.0], [-r0, 0.0]])
    v = np.array([[0.0, v0], [0.0, -v0]])

    def acc(x):
        d = x[1] - x[0]
        r = np.linalg.norm(d)
        f = m * m / r**2 * d / r  # 粒子0が受ける力（粒子1の向き）
        return np.array([f / m, -f / m])  # F21 = −F12（双極の F）

    def energy(x, v):
        r = np.linalg.norm(x[1] - x[0])
        return 0.5 * m * (v**2).sum() - m * m / r

    e0 = energy(x, v)
    if scheme == "euler":
        for _ in range(steps):
            a = acc(x)
            x = x + v * dt
            v = v + a * dt
        return energy(x, v) - e0
    # verlet：最初の「過去」x(−1) は誕生（初期条件）から渡す
    x_prev = x - v * dt + 0.5 * acc(x) * dt**2
    for _ in range(steps):
        x_next = 2 * x - x_prev + acc(x) * dt**2
        x_prev, x = x, x_next
    v = (x - x_prev) / dt  # 半時刻の速度（過去の量）
    return energy(x, v) - e0


de_euler = two_body("euler")
de_verlet = two_body("verlet")
print(f"  一時刻の式（オイラー）   エネルギーのずれ = {de_euler:+.3e}")
print(f"  三時刻の式（ベルレ）     エネルギーのずれ = {de_verlet:+.3e}")
print("  → 三時刻の式は時間を逆に回しても同じ形で、総和を長期間崩さない")

print()
print("=" * 72)
print("2. 誕生の段 T1〜T5")
print("=" * 72)
print(f"  {'段':>3} {'粒子N':>6} {'広さV':>6} {'隙間G':>6}  状態")
for k in range(1, 6):
    N = 2 * 3 ** (k - 1)
    V = k**3
    G = V - N
    state = "拘束（N > V）" if G < 0 else "隙間あり"
    print(f"  T{k:<2} {N:>6} {V:>6} {G:>6}  {state}")
print("  記録との照合：T3 の隙間 9、T4 の隙間 10、T5 で 37 個があふれる → 一致")
print("  T5 の体積拘束で3倍の増え方は止まり、以後は境界が光速で一辺1ずつ広がる")

print()
print("=" * 72)
print("3. 3^3 + 4^3 + 5^3 = 6^3 は T6 だけか")
print("=" * 72)
hits = [(a, a + 1, a + 2, a + 3) for a in range(1, 10000)
        if a**3 + (a + 1) ** 3 + (a + 2) ** 3 == (a + 3) ** 3]
print(f"  連続する3つの立方数の和が次の立方数になる a（1〜9999）：{hits}")
print(f"  6^3 − 5^3 = {6**3 - 5**3} = 3^3 + 4^3 = {3**3 + 4**3}")

print()
print("=" * 72)
print("4. 今の時刻と境界の圧力")
print("=" * 72)
t_P = 5.391247e-44            # プランク時間 s（CODATA 2018）
l_P = 1.616255e-35            # プランク長 m
age = 13.8e9 * 365.25 * 86400  # 宇宙年齢 s（Planck 2018 で 13.8 Gyr）
n_now = age / t_P
L_now = 5 + n_now
print(f"  今の時刻 n = {n_now:.3e}  （一辺 L = 5 + n）")

Lam = 1.1056e-52              # 宇宙定数 m^-2（Planck 2018, ΛCDM）
Lam_lP2 = Lam * l_P**2
rho_ratio = Lam_lP2 / (8 * math.pi)
print(f"  観測：Λ·ℓP^2 = {Lam_lP2:.3e}   ρΛ/ρP = Λ·ℓP^2/(8π) = {rho_ratio:.3e}")
print()
print(f"  {'読み方':<26} {'式':<12} {'今の値':>11} {'÷Λ·ℓP²':>9} {'÷ρΛ/ρP':>9}")
rows = [
    ("どこでも同じ強さ（定数）", "34", 34.0),
    ("全体積に分ける", "34/L^3", 34 / L_now**3),
    ("境界の力÷境界の面（採用）", "34/(6L^2)", 34 / (6 * L_now**2)),
    ("  同、元と新の半分", "17/(6L^2)", 17 / (6 * L_now**2)),
]
for name, f, val in rows:
    print(f"  {name:<24} {f:<12} {val:>11.3e} {val / Lam_lP2:>9.2e} {val / rho_ratio:>9.2e}")
print(f"  起点 n=0 の圧力 34/150 = {34 / 150:.3f}")

print()
print("  取り下げた読み方（34 を毎ステップ積み上がる数とした場合）")
print(f"    34n/L^3 = {34 * n_now / L_now**3:.3e}   (34/2000)n/L^3 = {34 / 2000 * n_now / L_now**3:.3e}")

print()
print("=" * 72)
print("5. 記号の関係")
print("=" * 72)
ok = []
for Fp, Fm, BI, pi, XYZ in itertools.product([1, -1], repeat=5):
    conds = [
        Fp + Fm == 0,          # F は必ず ± の組で、足すと 0
        Fp * Fm == -1,         # Fi × Fi = −1
        Fp * BI == -1,         # FBI = −1
        Fp * Fm * BI == 1,     # FFBI = 1
        Fp * BI * Fm == 1,     # FBIF = 1（電荷）
        Fp * BI * pi == 1,     # FBIπ = 1
        XYZ * pi == 1,         # XYZπ = 1
        Fp + BI == 0,          # B と I の回転で F と釣り合う
    ]
    if all(conds):
        ok.append(dict(Fp=Fp, Fm=Fm, BI=BI, pi=pi, XYZ=XYZ))
print(f"  すべての関係を満たす割り当て：{ok}")
B = np.array([1, 0, 0])
I = np.array([0, 1, 0])
F = np.cross(B, I)
print(f"  B×I = {F}（F 軸）  I×B = {np.cross(I, B)}  F + I×B = {F + np.cross(I, B)}")
print(f"  三重積 F·(B×I) = {F @ np.cross(B, I)}   F·(I×B) = {F @ np.cross(I, B)}")

print()
print("=" * 72)
print("6. 1/L^2 で薄まる圧力と p = −ρ の両立")
print("=" * 72)
print("  境界が光速で一定に進む：a ∝ n、H = 1/t（加速も減速もしない）")
print("  ρ ∝ L^-2 ∝ a^-2 を、出入りのない流体として扱うと w = −1 + 2/3 = −1/3")
print("  p = −ρ（w = −1）を保つには、連続の式 dρ/dt + 3H(ρ + p) = Q から")
print("  Q = dρ/dt = −2ρ/t：境界（UM=1 の誕生の側）へ、毎時刻 2ρ/t の割合で出ていく必要がある")
print("  JA0XUY-BPFF verify/bpff_boundary_w.py のケースC（粒子地平線）は、出入りなしで")
print("  w0 = +0.217 となり加速しない。出入り Q を TTT から導けるかが未解決")
