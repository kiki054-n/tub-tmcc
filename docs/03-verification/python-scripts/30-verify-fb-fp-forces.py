"""
TUB-TMCC 30：Ff＝Fb＋Fp ― 2粒子間の力の向きと大きさ、出たエネルギーの行き先（2026-11-10）

  python3 30-verify-fb-fp-forces.py

Fb：磁石（磁気モーメント）どうしの力。Fp：電荷の力（p は電流の記号。ここでは電荷の力として計算）
値：CODATA 2018、AME2020、文献値（本文の参考文献）
"""
import math
import numpy as np

k_e = 8.9875517923e9          # N m² / C²
e = 1.602176634e-19
k_m = 1e-7                    # μ0/4π
muN = 5.0507837461e-27
mu_p = 2.79284734463 * muN
mu_n = -1.91304273 * muN
h_eV = 4.135667696e-15
Eh = 27.211386245988          # eV
a0 = 0.529177210903e-10       # m
U = 931.49410242
line = "=" * 74


def f_dipole(m1, m2, r):
    """磁石1が磁石2に及ぼす力（r は 1→2 のベクトル）"""
    R = np.linalg.norm(r); n = r / R
    return 3 * k_m / R**4 * ((m1 @ n) * m2 + (m2 @ n) * m1 + (m1 @ m2) * n - 5 * (m1 @ n) * (m2 @ n) * n)


z = np.array([0.0, 0.0, 1.0]); x = np.array([1.0, 0.0, 0.0])
ARR = [("縦並び・同じ向き", z, z), ("縦並び・逆向き", z, -z), ("横並び・同じ向き", x, x), ("横並び・逆向き", x, -x)]

print(line)
print("1. 2粒子間の Fp と Fb（＋は反発、−は引き合う。軸方向の成分）")
for name, R, ma, mb, qq in [("水素分子 陽子−陽子 0.7414 Å", 0.7414e-10, mu_p, mu_p, 1),
                            ("ヘリウム核 陽子−陽子 2 fm", 2e-15, mu_p, mu_p, 1),
                            ("ヘリウム核 陽子−中性子 2 fm（スピンの向きで表示）", 2e-15, mu_p, mu_n, 0)]:
    Fp = k_e * qq * e * e / R**2
    print(f"  {name}")
    print(f"    Fp = {Fp:+.3e} N" + ("（反発）" if qq else "（中性子は電荷 0）"))
    for lab, a, b in ARR:
        Fz = f_dipole(ma * a, mb * b, R * z) @ z
        ratio = f"  Fb/Fp = {abs(Fz) / Fp:.2e}" if qq else ""
        print(f"    Fb {lab:10s} {Fz:+.3e} N{ratio}")

print(line)
print("2. 水素原子：押す（動き）と引く（電荷）")
F_a0 = k_e * e * e / a0**2
print(f"  ボーア半径での Fp = {F_a0:.4e} N、F×r = {F_a0 * a0 / e:.4f} eV")
print(f"  運動 +{Eh / 2:.3f} eV、電荷 −{Eh:.3f} eV、和 −{Eh / 2:.3f} eV（押す＝引くの半分）")

print(line)
print("3. 水素原子2つと水素分子（核を止めた計算、Kolos-Wolniewicz の E = −1.174475 Eh）")
E_H2 = -1.1744757 * Eh
E_2H = -1.0 * Eh
for nm, E in [("水素原子2つ", E_2H), ("水素分子", E_H2)]:
    print(f"  {nm:8s} 運動 +{-E:.2f} eV  電荷 {2 * E:.2f} eV  和 {E:.2f} eV")
print(f"  差 {E_H2 - E_2H:.3f} eV、振動の最低エネルギー 0.27 eV を引くと 解離 4.478 eV")
w = 4401.21 * 1.23984198e-4             # eV
mu_red = 0.5 * 1.007276 * U * 1e6 * e / (2.99792458e8**2)   # kg（換算質量）
omega = w / 6.582119569e-16
x_rms = math.sqrt(1.054571817e-34 / (2 * mu_red * omega))
print(f"  振動 {w:.4f} eV、振れ幅 ±{x_rms * 1e10:.3f} Å（距離の {x_rms / 0.7414e-10 * 100:.0f}%）")
F_pp = k_e * e * e / (0.7414e-10)**2
print(f"  陽子にかかる力：水素原子 {F_a0:.3e} N、水素分子（陽子どうし）{F_pp:.3e} N、比 {F_a0 / F_pp:.3f}")

print(line)
print("4. Fb のエネルギー（測定値）")
for nm, f in [("水素原子 21 cm 線", 1420.405751768e6), ("水素分子 陽子どうし（約）", 57.7e3), ("水素分子 回転と陽子（約）", 113.9e3)]:
    print(f"  {nm:22s} {f:.6g} Hz → {f * h_eV:.3e} eV")
print(f"  21 cm 線 ÷ 水素原子の結合 13.6 eV = {1420.405751768e6 * h_eV / 13.598:.2e}")

print(line)
print("5. 核分裂：箱を作るエネルギーと押し合いの解放（ワイツゼッカーの式、U236 を半分に）")
A, Z, a_s, a_c = 236, 92, 17.8, 0.711
Es1, Es2 = a_s * A**(2 / 3), 2 * a_s * (A / 2)**(2 / 3)
Ec1, Ec2 = a_c * Z * Z / A**(1 / 3), 2 * a_c * (Z / 2)**2 / (A / 2)**(1 / 3)
print(f"  表面（箱）  {Es1:.0f} → {Es2:.0f} MeV  差 +{Es2 - Es1:.0f}（入れる）")
print(f"  押し合い Fp {Ec1:.0f} → {Ec2:.0f} MeV  差 −{Ec1 - Ec2:.0f}（出る）")
print(f"  差し引き 約 {(Ec1 - Ec2) - (Es2 - Es1):.0f} MeV が出る")
for a in (236, 118):
    R = 1.2 * a**(1 / 3)
    print(f"  A={a}: 半径 {R:.2f} fm、密度 {a / (4 / 3 * math.pi * R**3):.3f} 個/fm³")
R1, R2 = 1.2 * 92**(1 / 3), 1.2 * 141**(1 / 3)
print(f"  Kr92＋Ba141 が接した瞬間の押し合い {1.44 * 36 * 56 / (R1 + R2):.0f} MeV（実際は変形して離れているので小さい。測定の運動エネルギー 約170 MeV）")

print(line)
print("6. 一つの規則：深い釣り合いへ移ると、深くなった分が出て軽くなる")
rows = [("水素原子 Fb 反転（21 cm）", 1420.405751768e6 * h_eV, 1.007825 * U * 1e6),
        ("水素＋酸素 → 水（1個）", 2.506, 18.0106 * U * 1e6),
        ("水素原子2 → 水素分子", 4.478, 2 * 1.007825 * U * 1e6),
        ("陽子＋電子 → 水素原子", 13.598, 1.007825 * U * 1e6),
        ("陽子＋中性子 → 重水素", 2.2246e6, 2.014102 * U * 1e6),
        ("D＋T → He4＋n", 17.589e6, (2.014102 + 3.016049) * U * 1e6),
        ("陽子2＋中性子2 → He4 核", 28.296e6, 3755.675e6),
        ("U235＋n → かけら2つ", 200e6, 236 * U * 1e6)]
for nm, E, M in rows:
    print(f"  {nm:24s} 出る {E:10.3e} eV   軽くなる割合 {E / M:.1e}")
print("  照合：重さの差と出た光のエネルギー 1 − Δmc²/E = (−1.4 ± 4.4)×10⁻⁷（Rainville 2005）")
print(line)
