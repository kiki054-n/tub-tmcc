"""
TUB-TMCC 29：質はエネルギーを量でためる ― 質量・エネルギー・力を数値でつなぐ（2026-11-10）

  python3 29-verify-mass-stores-energy.py

問い（本人）：どこに騙されているかを数値で表す。エネルギーと力と質量の関係。
              質はエネルギーを量でためることを示したい。

方法：どの物も「中身の静止質量」と「ためたエネルギー」に分ける。
      ためたエネルギー ＝ 物の質量 − 中身の静止質量の和（正ならためている、負なら放出済み）
      さらに、ためたエネルギーを「力 × 距離」と比べる。
値：CODATA 2018（u 単位の値を使い、u→MeV を一度だけ換算）、AME2020、PDG 2024（クォーク質量、MS-bar 2 GeV）
"""
import math

MeV = 1.0
eV = 1e-6
U = 931.49410242
# 原子質量と比べるため、電子・陽子・中性子も u で入れて同じ換算をする
# （MeV の値と u の値を混ぜると、u↔MeV 換算の不確かさ 約0.3 eV が水素原子の結合に入る）
me = 5.48579909065e-4 * U
mp = 1.007276466621 * U
mn = 1.00866491595 * U
mu, md = 2.16, 4.70                    # PDG 2024（スキームに依存する値）
H = 1.007825031898 * U
D = 2.014101777844 * U
He4 = 4.00260325413 * U
He4_nuc = He4 - 2 * me + 79.005 * eV
D_nuc = D - me + 13.6022 * eV
H2 = 2 * H - 4.4781 * eV
hbarc = 197.3269804                    # MeV·fm
MeV_per_fm_in_N = 1.602176634e-13 / 1e-15

line = "=" * 74
print(line)
print("1. 質量 ＝ 中身の静止質量 ＋ ためたエネルギー")
print(f"  {'物':22s}{'質量 MeV':>16s}{'中身の和 MeV':>16s}{'ためた量 MeV':>16s}{'質量に占める割合':>12s}")
rows = [
    ("陽子 (uud)",            mp,      2 * mu + md,   "中身＝クォーク3個"),
    ("中性子 (udd)",          mn,      mu + 2 * md,   "中身＝クォーク3個"),
    ("重水素核 (p+n)",        D_nuc,   mp + mn,       "中身＝陽子＋中性子"),
    ("ヘリウム核 (2p+2n)",    He4_nuc, 2 * mp + 2 * mn, "中身＝陽子2＋中性子2"),
    ("水素原子 (p+e)",        H,       mp + me,       "中身＝陽子＋電子"),
    ("水素分子 (2H)",         H2,      2 * H,         "中身＝水素原子2"),
]
for name, m, parts, note in rows:
    stored = m - parts
    s = f"{stored:16.6f}" if abs(stored) >= 1e-3 else f"{stored / eV:13.4f} eV"
    print(f"  {name:22s}{m:16.6f}{parts:16.6f}{s}{stored / m:12.2e}   {note}")
print("  陽子・中性子：質量の約99%が、クォークの静止質量ではなく、ためたエネルギー（運動と場）。")
print("  結合した核・原子・分子：ためた量が負。結びつくとき外へ出し、軽くなっている。")

print(line)
print("2. 量をためる向きと出す向き")
p_quarks = 2 * mu + md
print(f"  陽子：中身 {p_quarks:.2f} MeV に {mp - p_quarks:.2f} MeV をためて 938.27 MeV（閉じ込め＝ためる）")
print(f"  ヘリウム核：中身 {2*mp+2*mn:.2f} MeV から {2*mp+2*mn-He4_nuc:.4f} MeV を出して {He4_nuc:.4f} MeV（結合＝出す）")
print("  核融合で出る 26.73 MeV（docs/28）は、結合で出す側。結合した後のほうが、ためている量は少ない。")

print(line)
print("3. ためた量と力 × 距離")
# 水素原子：クーロン力とボーア半径（ビリアル定理：結合エネルギー ＝ ½ F r）
a0 = 0.529177210903e-10      # m
alpha = 7.2973525693e-3
hbarc_eVm = 197.3269804e-9 * 1e-6 / 1e-15 * 1e-15  # = 1.973e-7 eV·m
e2 = alpha * 1.973269804e-7  # e²/(4πε0) in eV·m
F_a0_eV_per_m = e2 / a0**2
Fr = F_a0_eV_per_m * a0
print("  水素原子（クーロン力）")
print(f"    ボーア半径での力 F       = {F_a0_eV_per_m * 1.602176634e-19:.4e} N")
print(f"    F × r                    = {Fr:.4f} eV")
print(f"    結合エネルギー           = {mp + me - H:.10f} MeV = {(mp + me - H) / eV:.4f} eV")
print(f"    結合エネルギー ÷ (F × r) = {((mp + me - H) / eV) / Fr:.4f}  （½ ＝ ビリアル定理）")
print(f"    ためた量（負）÷ 質量     = {-(mp + me - H) / H:.3e}")
# 陽子：クォークを結ぶ張力（弦の張力 σ ≈ 0.18 GeV²）
sigma = 0.18e6 / hbarc       # MeV/fm  （0.18 GeV² ÷ ħc）
r_p = 0.8414                 # fm（陽子の電荷半径 CODATA 2018）
print("  陽子（クォークを結ぶ張力）")
print(f"    張力 σ ≈ 0.18 GeV²       = {sigma:.0f} MeV/fm = {sigma * MeV_per_fm_in_N:.2e} N")
print(f"    張力 × 陽子の半径 {r_p} fm = {sigma * r_p:.0f} MeV（目安）")
print(f"    陽子にためた量           = {mp - p_quarks:.0f} MeV")
print("    桁が合う。ただし張力 × 半径は目安で、精密な値は格子 QCD の計算による。")

print(line)
print("4. 質量の差とエネルギーを別々に測った照合")
print("  Rainville ほか 2005（Nature 438, 1096）：Si・S の中性子捕獲")
print("    質量の差 Δm：ペニングトラップで重さを量る（MIT）")
print("    エネルギー E：γ線の波長を結晶で測る（NIST、ILL）")
print("    1 − Δm c² / E = (−1.4 ± 4.4) × 10⁻⁷")
print("  重さで量った差と、出てきた光で測ったエネルギーが、1000万分の数の精度で一致する。")

print(line)
print("5. どこに騙されているか（数値）")
print("  (a)「質量＝物の量（粒の数）」という見方")
print(f"      陽子の質量のうち、粒（クォーク）の静止質量は {p_quarks / mp * 100:.2f}%。")
print(f"      残り {100 - p_quarks / mp * 100:.2f}% はためたエネルギー。質量は粒の数ではなく、ためた量を含む。")
print("  (b)「質量がエネルギーに変換される」という言い方")
print("      ヘリウムができても粒の数（バリオン数4）は変わらない。減ったのは、ためていた量 28.30 MeV。")
print(f"      質量の {(2*mp+2*mn-He4_nuc)/(2*mp+2*mn)*100:.3f}% にあたる。粒が消えたのではない。")
print("  (c) 測定の側に騙しはない：§4 の照合で、重さの差と出たエネルギーは一致している。")
print(line)
