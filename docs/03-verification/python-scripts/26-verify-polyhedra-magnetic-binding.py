"""
TUB-TMCC 26：多面体の要素数・双錘のねじれ・磁気モーメント・見えない力の検証（2026-11-08 夜）

この記録の数値はすべてこのスクリプトで出す。
  python3 26-verify-polyhedra-magnetic-binding.py

1. 要素数 = 面の角の合計 = 頂点の線の合計 = 2E（プラトン5立体と n角双錘）
2. 中心からの距離（正三角形の双錘 / 球に内接させた双錘）
3. 3本（120°）と4本（正四面体）の向きの和、正8面体を対角線から見た120°
4. 陽子・中性子の磁気モーメント：回転軸と4つ目の符号の8通り
5. 臍の寄与と重水素の磁気モーメント
6. 結合エネルギー（重水素・ヘリウム4・水素分子）と1本あたりの比
7. 質と量を分けた分母（読み方A / 読み方B）
8. 見えない力：水素分子のビリアル定理による内訳
"""
import itertools
import math

import numpy as np

line = "=" * 72

# ---------------------------------------------------------------------------
print(line)
print("1. 要素数 = 2E")
print(line)
solids = {
    # name: (faces as list of polygon sizes, vertex degrees as list)
    "正四面体": ([3] * 4, [3] * 4),
    "正6面体": ([4] * 6, [3] * 8),
    "正8面体": ([3] * 8, [4] * 6),
    "正12面体": ([5] * 12, [3] * 20),
    "正20面体": ([3] * 20, [5] * 12),
}
for n in (3, 4, 5):
    # n角双錘：面 2n（三角形）、極点2（次数 n）、赤道 n点（次数 4）
    solids[f"{n}角双錘"] = ([3] * (2 * n), [n] * 2 + [4] * n)
for name, (faces, degs) in solids.items():
    V, F = len(degs), len(faces)
    s_face, s_vert = sum(faces), sum(degs)
    E = s_face // 2
    assert s_face == s_vert, name
    assert V - E + F == 2, name
    print(f"{name:6s} V={V:2d} E={E:2d} F={F:2d}  面で数える={s_face:2d}  点で数える={s_vert:2d}")
print("接続点を3本と数えた3角双錘：", 3 * 2 + 3 * 3, "（面の18と合わない）")

# ---------------------------------------------------------------------------
print()
print(line)
print("2. 中心からの距離（辺 a = 1 の正三角形の双錘）")
print(line)
for n in (3, 4, 5):
    R = 1 / (2 * math.sin(math.pi / n))  # 赤道の外接半径
    h = math.sqrt(1 - R * R)  # 極点の高さ
    print(f"{n}角双錘  極点 {h:.3f}  接続点 {R:.3f}")
print("球（半径1）に内接させた場合の辺の長さ：")
for n in (3, 5):
    print(f"{n}角双錘  赤道の辺 {2 * math.sin(math.pi / n):.3f}  極点への辺 {math.sqrt(2):.3f}")

# ---------------------------------------------------------------------------
print()
print(line)
print("3. 3本のY・4本の正四面体・正8面体の対角線方向")
print(line)
Y = [np.array([math.cos(t), math.sin(t), 0.0]) for t in (0, 2 * math.pi / 3, 4 * math.pi / 3)]
r = math.sqrt(8) / 3
T = [np.array([0, 0, 1.0])] + [np.array([r * math.cos(t), r * math.sin(t), -1 / 3]) for t in (0, 2 * math.pi / 3, 4 * math.pi / 3)]
print("Y（120°）の和 :", np.round(sum(Y), 12))
print("正四面体4本の和:", np.round(sum(T), 12), " 隣り合う角度", round(math.degrees(math.acos(-1 / 3)), 2), "°")
d = np.array([1, 1, 1]) / math.sqrt(3)
proj = [e - e.dot(d) * d for e in np.eye(3)]
ang = math.degrees(math.acos(proj[0].dot(proj[1]) / (np.linalg.norm(proj[0]) * np.linalg.norm(proj[1]))))
print("+X,+Y,+Z を対角線方向から見たときの角度:", round(ang, 6), "°")

# ---------------------------------------------------------------------------
print()
print(line)
print("4. 磁気モーメント：8通り")
print(line)
mu_p_obs, mu_n_obs = 2.79284734, -1.91304273
obs = mu_n_obs / mu_p_obs
z = np.array([0, 0, 1.0])
print("観測の比", round(obs, 4))
for pax, nax, apex in itertools.product(["共通軸", "外向き"], ["共通軸", "外向き"], [1, -1]):
    mp = sum((z if pax == "共通軸" else v) for v in Y)
    w = [apex, 1, 1, 1]
    mn = sum(wi * (z if nax == "共通軸" else v) for wi, v in zip(w, T))
    P, N = mp[2], mn[2]
    ratio = f"{N / P:+.4f}" if abs(P) > 1e-9 else "求まらない"
    print(f"陽子={pax} 中性子={nax} 4つ目={apex:+d}  陽子={P:+.3f} 中性子={N:+.3f} 横成分={np.hypot(mn[0], mn[1]):.0e} 比={ratio}")
m = mu_p_obs / 3
print(f"陽子に合わせた1束 m = {m:.4f}  中性子モデル = {-2 * m:.4f}（観測 {mu_n_obs}）  比のずれ {(-2 / 3 - obs) / obs * 100:+.2f}%")

# ---------------------------------------------------------------------------
print()
print(line)
print("5. 臍の寄与と重水素")
print(line)
m2 = -mu_n_obs / 2
navel = mu_p_obs - 3 * m2
print(f"1束 m = {m2:.4f}  陽子の3束 = {3 * m2:.4f}  臍 = {navel:+.4f}（陽子の {navel / mu_p_obs * 100:+.2f}%）")
mu_d_obs = 0.8574382
print(f"重水素：モデル {mu_p_obs + mu_n_obs:.4f}  観測 {mu_d_obs}  差 {mu_d_obs - (mu_p_obs + mu_n_obs):+.4f}")

# ---------------------------------------------------------------------------
print()
print(line)
print("6. 結合エネルギー")
print(line)
mp_, mn_, me = 938.27208816, 939.56542052, 0.51099895  # MeV
md, ma = 1875.61294257, 3727.3794066
Bd = mp_ + mn_ - md
Ba = 2 * mp_ + 2 * mn_ - ma
mH = mp_ + me - 13.598434e-6
D0, De = 4.4781e-6, 4.7474e-6  # MeV
print(f"重水素   {Bd:.4f} MeV  割合 {Bd / (mp_ + mn_):.4e}")
print(f"ヘリウム4 {Ba:.4f} MeV  割合 {Ba / (2 * mp_ + 2 * mn_):.4e}")
print(f"水素分子 {D0 * 1e6:.4f} eV  割合 {D0 / (2 * mH):.4e}")
print(f"ヘリウム4/重水素 = {Ba / Bd:.3f}  1本あたり {Bd:.3f} / {Ba / 6:.3f} MeV  比 {Ba / 6 / Bd:.4f}  （3/√2 = {3 / math.sqrt(2):.4f}）")
print(f"割合の比 重水素/水素分子 = {(Bd / (mp_ + mn_)) / (D0 / (2 * mH)):.3e}")

# ---------------------------------------------------------------------------
print()
print(line)
print("7. 質と量を分けた分母")
print(line)
unit = me / 500  # 読み方A：1要素 = 電子の質量の1/500
q_p, q_n = 1500 * unit, 2000 * unit
print(f"読み方A：1要素 {unit * 1e3:.3f} keV  陽子 質 {q_p:.3f} 量 {mp_ - q_p:.3f}  中性子 質 {q_n:.3f} 量 {mn_ - q_n:.3f} MeV")
print(f"  重水素   結合/量 {Bd / (mp_ - q_p + mn_ - q_n):.4e}")
print(f"  ヘリウム4 結合/量 {Ba / (2 * (mp_ - q_p) + 2 * (mn_ - q_n)):.4e}")
print(f"読み方B：要素数を電子の質量の倍率と読むと 中性子の量 = {mn_ / me:.2f} − 2000 = {mn_ / me - 2000:.2f}（マイナス）")

# ---------------------------------------------------------------------------
print()
print(line)
print("8. 見えない力（水素分子のビリアル定理）")
print(line)
Eh = 27.211386  # eV
E_H2 = -1.17447 * Eh  # 平衡距離での全エネルギー（核の振動を除く）
print(f"全エネルギー {E_H2:.2f} eV → 運動 {-E_H2:+.2f} eV、位置 {2 * E_H2:+.2f} eV")
print(f"原子2個（{-Eh:.2f} eV）からの差 {E_H2 + Eh:+.3f} eV = De")
print("重水素の内訳（核力模型 AV18 の計算値、参考）：運動 約 +19.8 MeV、位置 約 −22.0 MeV、差 −2.22 MeV")
