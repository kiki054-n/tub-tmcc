# TUB-TMCC リポジトリ ディレクトリ構造

> トリテトラ理論（TTT）3要素方程式による結合・物質状態・組成設計の記録リポジトリ

---

## 推奨ディレクトリ構造

```
TUB-TMCC/
├── README.md                          # リポジトリの概要・核心原理
├── LICENSE                            # CC BY-NC-SA 4.0
├── CITATION.cff                       # 引用情報
│
├── docs/                              # 理論ドキュメント（核心・基礎）
│   ├─ 00-philosophy/                 # 哲学・原理
│   │   └── origin-of-zero.md          # 0の起源：双極と3つの力
│   │
│   ├── 01-foundations/                # 基礎方程式・定義
│   │   ├── three-field-equation.md    # 3要素の万有引力の法則
│   │   ├── mass-and-energy.md         # 質（M）と量（E）の定義
│   │   ├── cm-and-ce.md             # 質の中心と量の中心
│   │   └── constants.md               # 定数一覧（G_M, G_E, G_ME, G_R）
│   │
│   ├── 02-mathematical-framework/     # 数学的枠組み
│   │   ├── vector-analysis.md         # ベクトル分解（大きさ・向き・方向）
│   │   ├── geometry.md              # 幾何学的構造（Tri-Tetra）
│   │   └── 3-4-5-6-relation.md      # 3³+4³+5³=6³ の宇宙原理
│   │
│   └── 03-verification/               # 検証方法論
│       ├── methodology.md             # 検証の手法（仮説実証→失敗→評価）
│       └── python-scripts/            # 検証用Pythonスクリプト
│           ├── calculate_forces.py
│           ├── thermal_expansion.py
│           └── bonding_analysis.py
│
├── elements/                          # 周期律表：元素ごとの記録
│   ├── README.md                      # 元素記録のガイドライン
│   │
│   ├── period-1/                      # 第1周期（1段目）
│   │   ├── 001-H/                     # 水素 Hydrogen
│   │   │   ├── element.md             # 元素の基本データ（質・量・Cm・Ce）
│   │   │   ├── atomic-structure.md    # 原子構造の3要素解析
│   │   │   ├── bonding/             # 結合の記録
│   │   │   │   ├── H-H.md           # H₂分子
│   │   │   │   ├── H-O.md           # OH（水酸基）
│   │   │   │   └── H-C.md           # CH（メタン等）
│   │   │   └── states/              # 物質状態の記録
│   │   │       ├── solid.md           # 固体水素
│   │   │       ├── liquid.md          # 液体水素
│   │   │       └── plasma.md        # 水素プラズマ
│   │   │
│   │   └── 002-He/                    # ヘリウム Helium
│   │       ├── element.md
│   │       ├── atomic-structure.md
│   │       ├── bonding/
│   │       │   └── He-He.md         # He₂（不安定）
│   │       └── states/
│   │           ├── solid.md
│   │           ├── liquid.md
│   │           └── plasma.md
│   │
│   ├── period-2/                      # 第2周期（2段目）
│   │   ├── 003-Li/                    # リチウム Lithium
│   │   │   ├── element.md
│   │   │   ├── atomic-structure.md
│   │   │   ├── bonding/
│   │   │   │   ├── Li-Li.md
│   │   │   │   ├── Li-H.md
│   │   │   │   ├── Li-O.md
│   │   │   │   └── Li-F.md
│   │   │   └── states/
│   │   │       ├── solid.md
│   │   │       ── liquid.md
│   │   │       └── plasma.md
│   │   │
│   │   ├── 004-Be/                    # ベリリウム Beryllium
│   │   │   ├── element.md
│   │   │   ├── atomic-structure.md
│   │   │   ├── bonding/
│   │   │   │   ├── Be-Be.md
│   │   │   │   ├── Be-O.md
│   │   │   │   └── Be-F.md
│   │   │   └── states/
│   │   │       ├── solid.md
│   │   │       ├── liquid.md
│   │   │       └── plasma.md
│   │   │
│   │   ├── 005-B/                     # ホウ素 Boron
│   │   │   ├── element.md
│   │   │   ├── atomic-structure.md
│   │   │   ├── bonding/
│   │   │   │   ├── B-B.md
│   │   │   │   ├── B-O.md
│   │   │   │   └── B-N.md
│   │   │   └── states/
│   │   │       └── solid.md
│   │   │
│   │   ├── 006-C/                     # 炭素 Carbon
│   │   │   ├── element.md
│   │   │   ├── atomic-structure.md
│   │   │   ├── bonding/
│   │   │   │   ├── C-C.md           # 共有結合
│   │   │   │   ├── C-H.md           # CH（有機化合物の基礎）
│   │   │   │   ├── C-O.md           # CO, CO₂
│   │   │   │   └── C-N.md           # CN（シアン基）
│   │   │   └── states/
│   │   │       ├── diamond.md         # ダイヤモンド（3次元共有結合）
│   │   │       ├── graphite.md        # 黒鉛（2次元共有結合）
│   │   │       ├── graphene.md        # グラフェン（単層）
│   │   │       └── plasma.md
│   │   │
│   │   ├── 007-N/                     # 窒素 Nitrogen
│   │   │   ├── element.md
│   │   │   ├── atomic-structure.md
│   │   │   ├── bonding/
│   │   │   │   ├── N-N.md           # N₂（三重結合）
│   │   │   │   ├── N-H.md           # NH₃（アンモニア）
│   │   │   │   ├── N-O.md           # NO, NO₂, N₂O
│   │   │   │   └─ N-C.md           # シアン化物
│   │   │   └── states/
│   │   │       ├── solid.md
│   │   │       ├── liquid.md
│   │   │       └── plasma.md
│   │   │
│   │   ├── 008-O/                     # 酸素 Oxygen
│   │   │   ├── element.md
│   │   │   ├── atomic-structure.md
│   │   │   ├── bonding/
│   │   │   │   ├── O-O.md           # O₂（二重結合）
│   │   │   │   ├── O-H.md           # OH（水酸基）
│   │   │   │   ├── H2O.md           # 水（特別扱い）
│   │   │   │   ├── O-C.md           # CO, CO₂
│   │   │   │   └── O-N.md           # NO, NO₂
│   │   │   └── states/
│   │   │       ├── ice.md             # 氷（四面体構造）
│   │   │       ├── water.md           # 水（水素結合）
│   │   │       └── plasma.md
│   │   │
│   │   ├── 009-F/                     # フッ素 Fluorine
│   │   │   ├── element.md
│   │   │   ├── atomic-structure.md
│   │   │   ├── bonding/
│   │   │   │   ├── F-F.md
│   │   │   │   ├── F-H.md           # HF
│   │   │   │   ├── F-Li.md          # LiF（イオン結合）
│   │   │   │   └── F-C.md           # CF₄（フロン）
│   │   │   └── states/
│   │   │       ├── solid.md
│   │   │       ├── liquid.md
│   │   │       └── plasma.md
│   │   │
│   │   └── 010-Ne/                    # ネオン Neon
│   │       ├── element.md
│   │       ├── atomic-structure.md
│   │       ├── bonding/
│   │       │   └── Ne-...md         # ほぼ反応しない（ファンデルワールス）
│   │       └── states/
│   │           ├── solid.md
│   │           ├── liquid.md
│   │           └── plasma.md
│   │
│   ├── period-3/                      # 第3周期（3段目）
│   │   ├── 011-Na/                    # ナトリウム
│   │   ├── 012-Mg/                    # マグネシウム
│   │   ├── 013-Al/                    # アルミニウム
│   │   ├── 014-Si/                    # ケイ素
│   │   ├── 015-P/                     # リン
│   │   ├── 016-S/                     # 硫黄
│   │   ├── 017-Cl/                    # 塩素
│   │   └── 018-Ar/                    # アルゴン
│   │
│   └── period-4-plus/                 # 第4周期以降（必要に応じて追加）
│       └── ...
│
├── compounds/                         # 化合物：重要な化合物の詳細解析
│   ├── README.md                      # 化合物記録のガイドライン
│   │
│   ├── inorganic/                     # 無機化合物
│   │   ├── oxides/                    # 酸化物
│   │   │   ├── H2O.md                 # 水（最重要）
│   │   │   ├── CO2.md                 # 二酸化炭素
│   │   │   ├── SiO2.md                # 二酸化ケイ素（石英
│   │   │   └── Al2O3.md             # 酸化アルミニウム
│   │   │
│   │   ├── salts/                     # 塩
│   │   │   ├── NaCl.md                # 塩化ナトリウム
│   │   │   ├── LiF.md                 # フッ化リチウム
│   │   │   └── CaCO3.md             # 炭酸カルシウム
│   │   │
│   │   └── superconductors/           # 超電導物質
│   │       ├── YBCO.md                # YBa₂Cu₃O₇（ペロブスカイト）
│   │       ├── BSCCO.md             # Bi-Sr-Ca-Cu-O系
│   │       ├── LK99.md                # 鉛-アパタイト系（検証中）
│   │       └── perovskite-ABO3.md   # ペロブスカイト一般式
│   │
│   └── organic/                       # 有機化合物
│       ├── hydrocarbons/              # 炭化水素
│       │   ├── methane.md             # メタン CH₄
│       │   ├── ethane.md              # エタン C₂H₆
│       │   └── benzene.md           # ベンゼン C₆H₆
│       │
│       ├── alcohols/                  # アルコール
│       │   └── ethanol.md           # エタノール C₂H₅OH
│       │
│       └── biomolecules/              # 生体分子
│           ├── amino-acids.md       # アミノ酸
│           └── DNA.md               # DNA（水素結合の集大成）
│
├── materials/                         # 材料：物質状態・相転移・特性
│   ├── README.md                      # 材料記録のガイドライン
│   │
│   ├── phases/                        # 相と相転移
│   │   ├── solid.md                   # 固体：3つの力が強く釣り合う
│   │   ├── liquid.md                  # 液体：3つの力が揺らぐ
│   │   ├── gas.md                     # 気体：3つの力が弱い
│   │   ├── plasma.md                  # プラズマ：3つの場が分離
│   │   └── phase-transitions.md     # 相転移の統一的説明
│   │
│   ├── bonding-types/                 # 結合の種類
│   │   ├── covalent.md                # 共有結合（d ≈ 0）
│   │   ├── ionic.md                   # イオン結合（d > 0）
│   │   ├── metallic.md                # 金属結合（F_Eの集団共鳴）
│   │   ├── hydrogen.md                # 水素結合（F_R支配）
│   │   └── vanderwaals.md           # ファンデルワールス力
│   │
│   ├── thermal-properties/            # 熱的性質
│   │   ├── thermal-expansion.md     # 熱膨張（温度振動と距離）
│   │   ├── specific-heat.md         # 比熱（量の蓄積）
│   │   └── ice-volume-expansion.md  # 氷の体積膨張（対称性変化）
│   │
│   └── superconductivity/           # 超電導
│       ├── theory.md                  # 3要素方程式による超電導理論
│       ├── composition-design.md    # 組成設計の指針
│       └── verification-data.md     # 実験データとの照合
│
├── phenomena/                         # 現象：自然現象の3要素解析
│   ├── nuclear-fusion.md            # 核融合（質と量の変）
│   ├── lightning.md                 # 雷（プラズマの創発）
│   ├── solar-physics.md             # 太陽（プラズマの3つの場）
│   └── water-anomalies.md           # 水の異常性（4°Cで密度最大等）
│
├── calculations/                      # 計算記録
│   ├── force-calculations/          # 力の計算
│   │   ├── H-atom/                  # 水素原子
│   │   ├── H2-molecule/           # 水素分子
│   │   ├── He-atom/               # ヘリウム原子
│   │   └── ...
│   │
│   ├── temperature-effects/         # 温度効果
│   │   ├── Cu-thermal-expansion/  # 銅の熱膨張
│   │   └── ice-formation/         # 氷の形成
│   │
│   └── python-scripts/              # 計算用スクリプト
│       ├── constants.py             # 物理定数・TTT定数
│       ├── force_calculator.py    # 3要素力計算機
│       ├── bonding_analyzer.py    # 結合解析ツール
│       └── thermal_simulator.py   # 温度シミュータ
│
├── experiments/                       # 実験記録
│   ├── proposed-experiments.md      # 提案実験
│   └── data-comparison/           # 実験データとの比較
│
├── discussions/                       # 議論記録
│   ├── open-questions.md            # 未解決問題
│   ├── research-log.md              # 研究ログ（日記形式）
│   └── collaborators/             # 協力者からの貢献
│
└── references/                        # 参考文献・関連理論
    ├── TTT-core.md                  # TTTコアリポジトリへのリンク
    ├── standard-physics.md          # 標準物理との対応
    └── bibliography.md              # 文献一覧
```

---

## ファイ命名規則

### 元素ファイル

```
elements/period-{N}/{XXX}-{Symbol}/
├── element.md           # 元素の基本データ
├── atomic-structure.md  # 原子構造の3要素解析
├── bonding/             # 結合の記録
│   ├── {A}-{B}.md       # A-B結合の解析
│   └── ...
└── states/              # 物質状態の記録
    ├── solid.md
    ├── liquid.md
    └── plasma.md
```

### 化合物ファイル

```
compounds/{category}/{formula}.md
```

例：
- `compounds/inorganic/oxides/H2O.md`
- `compounds/inorganic/superconductors/YBCO.md`
- `compounds/organic/hydrocarbons/methane.md`

### 計算ファイル

```
calculations/{type}/{subject}/
├── calculation.md       # 計算の記録
├── parameters.md        # パラメータ一覧
├── results.md           # 結果
└── scripts/
    ├── calculate.py
    └── plot.py
```

---

## 各MDファイルのテンプレート構造

### element.md のテンプレート

```markdown
# {元素名} ({Symbol}) - 原子番号 {Z}

## 基本データ

| 項目 | 値 | 説明 |
|------|-----|------|
| 原子番号 | {Z} | |
| 原子量 | {mass} u | |
| 電子配置 | {config} | |
| TTT質（M） | {M} | 要素数 |

## 3要素解析

### 質の中心（Cm）

### 量の中心（Ce）

### 中心ずれ（d）

## 結合特性

## 物質状態

## 計算結果

## 未解決・次の課題
```

### {A}-{B}.md（結合ファイル）のテンプレート

```markdown
# {A}-{B} 結合の3要素解析

## 結合の種類

## パラメータ

| パラメータ | A側 | B側 |
|-----------|-----|-----|
| 質（M） | | |
| 量（E） | | |
| Cm | | |
| Ce | | |
| d | | |

## 力の計算

### F_M（質的引力）

### F_E（量的引力）

### F_ME（交差引力）

### F_R（関係項）

## 釣り合い距離

## 実験値との比較

## 結論
```

---

## ライセンス

**CC BY-NC-SA 4.0**

本リポジトリのすべてのファイルは、表示・非営利・継承の条件で共有・改変できます。
商用利用は禁止されています。

---

## 著者

**川上真潔（Kawakami Naoyuki）**
- ORCID: 0009-0009-2972-6511
- Company: 友朋共生

---

> 草案段階。ディレクトリ構造は必要に応じて変更・拡張してください。
> みんなで直していけば、それがいい。
