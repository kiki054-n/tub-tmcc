"""
TUB-TMCC 31：Ff＝Fb＋Fe のベクトル解析と材料研究、F のエネルギー変換（2026-11-10）

  python3 31-verify-ff-fb-fe-materials.py

Fb：求心（引き合う）側の力、Fe：電子雲・環境の斥力（本人の「水分子生成反応の内積分解」の記号）
値：標準的な平均結合エネルギー、CODATA、文献値。「要確認」と書いた値は記憶による概略値。
"""
import math

line = "=" * 74
kJ_per_eV = 96.485332          # kJ/mol per eV/個
R_gas = 8.314462618

# ---------------------------------------------------------------------------
print(line)
print("1. 水の生成 H2 + 1/2 O2 -> H2O(g) を結合エネルギーで")
D_HH, D_OO2, D_OH = 436.0, 498.0, 463.0      # kJ/mol（平均結合エネルギー）
react = D_HH + 0.5 * D_OO2
prod = 2 * D_OH
dH = react - prod
print(f"  反応物 Σ D = {D_HH} + {0.5*D_OO2} = {react:.0f} kJ/mol")
print(f"  生成物 Σ D = 2 x {D_OH} = {prod:.0f} kJ/mol")
print(f"  ΔH = {dH:.0f} kJ/mol（測定値 ΔfH H2O(g) = -241.83 kJ/mol）")
print(f"  液体の水なら ΔfH = -285.83 kJ/mol（差 44.0 kJ/mol は蒸発熱）")
print(f"  ΔH と ΔU の差 Δn RT = -0.5 x RT = {-0.5*R_gas*298.15/1000:.2f} kJ/mol")

# ---------------------------------------------------------------------------
print(line)
print("2. Pauling の関係で分解")
chi_O, chi_H = 3.44, 2.20
dchi = chi_O - chi_H
ionic = 96.5 * dchi**2
D_OO1 = 146.0                                # O–O 単結合 kJ/mol
D_OH_est = 0.5 * (D_HH + D_OO1) + ionic
print(f"  Δχ = {dchi:.2f}、極性項 96.5 Δχ² = {ionic:.1f} kJ/mol（1本）")
print(f"  D(O–H) 推定 = ½(436+146) + {ionic:.0f} = {D_OH_est:.0f} kJ/mol（平均値 463）")
part1 = 0.5 * D_OO2 - D_OO1
part2 = -2 * ionic
print(f"  ΔH ≈ [½D(O=O) − D(O–O)] + [−2×極性項] = {part1:+.0f} {part2:+.0f} = {part1+part2:+.0f} kJ/mol")

# ---------------------------------------------------------------------------
print(line)
print("3. 斥力ベクトル2本の合力（同じ大きさ 1、なす角 θ）")
print("   |F_total|² = 1 + 1 + 2cosθ  →  |F_total| = 2|cos(θ/2)|")
for th in [90.0, 104.5, 109.47, 120.0, 150.0, 180.0]:
    c = math.cos(math.radians(th))
    print(f"   θ = {th:7.2f}°  cosθ = {c:+.3f}  |F_total| = {math.sqrt(2+2*c):.3f}")
print("   → 合力が最小（0）になるのは 180°（直線形）。折れ線 104.5° は最小ではない。")
print("   → 合力の大きさの二乗はエネルギーではない。エネルギーは ∫F·dr。")

barrier_cm = 11100.0                          # 要確認：水の線形化障壁（概略値）
cm_to_kJ = 0.01196266                          # 1 cm⁻¹ = 0.01196266 kJ/mol
print(f"  水の線形化障壁 約 {barrier_cm:.0f} cm⁻¹ = {barrier_cm*cm_to_kJ:.0f} kJ/mol（要確認）")

# ---------------------------------------------------------------------------
print(line)
print("4. 階層の分離：原子核と化学結合")
eV_bond = D_OH / kJ_per_eV
E_exc = 6.049e6                               # 酸素16 最初の励起準位 0+ 6.049 MeV
print(f"  O–H 1本 = {D_OH} kJ/mol = {eV_bond:.2f} eV")
print(f"  酸素16 の最初の励起準位 = 6.049 MeV、比 = {E_exc/eV_bond:.2e}")
r_ch = 2.699e-15                              # 酸素16 電荷半径 m
r_OH = 0.9572e-10
print(f"  酸素16 電荷半径 2.699 fm ÷ O–H 0.9572 Å = {r_ch/r_OH:.2e}")
E_hb = 20.0                                   # 水素結合（概略）
print(f"  水素結合 約 {E_hb:.0f} kJ/mol = {E_hb/kJ_per_eV:.2f} eV（O–H の約 1/{D_OH/E_hb:.0f}）")

# ---------------------------------------------------------------------------
print(line)
print("5. 四面体の網目と体積")
rho_ice, rho_w = 0.9167, 0.9998               # g/cm³（0 °C）
print(f"  氷 Ih {rho_ice}、水 {rho_w} g/cm³ → 凍ると体積 {100*(rho_w/rho_ice-1):+.2f} %")
rho_Si_s, rho_Si_l = 2.33, 2.57
print(f"  Si 固体 {rho_Si_s}、液体 {rho_Si_l} g/cm³ → 固まると体積 {100*(rho_Si_l/rho_Si_s-1):+.1f} %")
tet = math.degrees(math.acos(-1/3))
print(f"  四面体角 {tet:.2f}°、正五角形の内角 108.00°（差 {tet-108:.2f}°）、正六角形 120.00°")

# ---------------------------------------------------------------------------
print(line)
print("6. β-FeSi2 と Si の格子の合い方")
a_Si = 5.431
b, c = 7.791, 7.833                           # β-FeSi2（Cmca）Å
sq = math.sqrt(2) * a_Si
print(f"  √2 × a_Si = {sq:.3f} Å")
print(f"  β の b = {b} Å → ずれ {100*(b-sq)/sq:+.2f} %")
print(f"  β の c = {c} Å → ずれ {100*(c-sq)/sq:+.2f} %")
for name, a in [("γ-FeSi2（蛍石型）", 5.389), ("CoSi2", 5.365), ("NiSi2", 5.406)]:
    print(f"  {name:16s} a = {a} Å → Si とのずれ {100*(a-a_Si)/a_Si:+.2f} %")

# ---------------------------------------------------------------------------
print(line)
print("7. F のエネルギー変換の例")
M = 18.015e-3                                 # kg/mol
dV = M/(rho_ice*1000) - M/(rho_w*1000)        # m³/mol
L = 6010.0                                    # J/mol 融解熱
P0 = 101325.0
print(f"  凍るときの体積増 ΔV = {dV*1e6:.3f} cm³/mol")
print(f"  大気を押しのける仕事 P0ΔV = {P0*dV:.3f} J/mol、融解熱 {L:.0f} J/mol の {P0*dV/L:.1e} 倍")
dPdT = -L / (273.15 * dV)
print(f"  クラペイロン dP/dT = −L/(T ΔV) = {dPdT/1e6:.1f} MPa/K（融点を 1 K 下げるのに要る圧力）")
print(f"  融解熱 1 分子あたり = {L/1000/kJ_per_eV:.3f} eV")

print("  熱電変換の効率（ZT と温度から。FeSi2 の ZT 値ではない）")
Th, Tc = 900.0, 300.0
carnot = 1 - Tc/Th
for ZT in [0.2, 0.5, 1.0]:
    m = math.sqrt(1 + ZT)
    eta = carnot * (m - 1) / (m + Tc/Th)
    print(f"    Th={Th:.0f} K, Tc={Tc:.0f} K, ZT={ZT}: η = {100*eta:.1f} %（カルノー {100*carnot:.1f} %）")
print(line)
