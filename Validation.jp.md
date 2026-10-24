# 検証状況

本家ISODISTORT の出力との比較と、独立した数学的検査を用いて検証しています。
本書の数値は **2026 年 10 月 23 日時点**のスナップショットです。

## 有効母集団の定義

- 表示された完全なmode basisの有理数上の線形独立性
- 表示mode fieldの報告部分群作用に対する不変性

の2つの独立検査を突破した本家出力のみを有効母集団とします（現時点で**215件**が除外されています。大半はパラメータが1/2を超えるK由来です。）

## 判定の定義

大きく3種に区別しています。

- **Identical**: OPD の表示 setting、子構造、mode の種類・本数・原子対応・正規化係数・vector が、順序等のBYU固有の規約まで含め本家出力と完全に一致します。
- **Equivalent**: 明示的な結晶学的・mode 空間の同値（mode の順序、符号、同じ mode 空間内の基底、等価な setting）によって、子構造と張られる mode 空間が 本家出力と同値であると証明されたものです。Identical はすべて Equivalent に含まれます。
- **Unproven**: 比較器が同値を証明できなかったものです。物理的に異なる出力のほか、比較器が同値を判定しきれない場合です。APOSTRUCTの出力が誤りであるという判定ではありません。

## 検証母集団

### 1. MP pool: 広域ストレス母集団

Materials Project 由来の構造を基礎に、固定・parametric K、multi-K、一般 OPD を広く含む母集団です。

#### 構造 pool の選定

CIF pool は Materials Project から作成しました。
同じ空間群で、占有される Wyckoff site の multiplicity・letter の組が同じ構造は同じ topology とみなし、`(space group, occupied Wyckoff pattern)` ごとに1件だけを取得。
元素の置換だけが異なる構造や、同じ Wyckoff pattern の自由座標だけが異なる構造は、この段階では重複として除きました。
この方法で得た pool は **19,086 CIF** です。
以下の母集団で使用されたのは、現時点で**5,483 CIF** です。

#### 結果

| 項目 | 件数 |
|---|---:|
| All | 9,634 |
| Identical | 8,844 |
| Equivalent | 9,630 |
| Unproven | 4 |

- **Equivalent: 9,630 / 9,634 = 99.96%**
- **Identical: 8,844 / 9,634 = 91.80%**

#### K-signature 別

`F` は固定 K、`P` は parametric K を表します。

| K | Validation母集団 | Equivalent | Identical |
|---|---:|---:|---:|
| F | 6,233 | 6,233 (100.00%) | 5,913 (94.87%) |
| FF | 846 | 846 (100.00%) | 742 (87.71%) |
| FFF | 189 | 189 (100.00%) | 167 (88.36%) |
| FFFF | 202 | 202 (100.00%) | 186 (92.08%) |
| P | 1,287 | 1,287 (100.00%) | 1,169 (90.83%) |
| PP | 20 | 18 (90.00%) | 7 (35.00%) |
| FP | 126 | 126 (100.00%) | 108 (85.71%) |
| FFP | 453 | 453 (100.00%) | 364 (80.35%) |
| FFFP | 246 | 246 (100.00%) | 171 (69.51%) |
| FPP | 31 | 29 (93.55%) | 17 (54.84%) |
| FFPP | 1 | 1 (100.00%) | 0 (0.00%) |

### MAGNDATA: practical 母集団

MAGNDATA に収録された実在磁気構造を基にした母集団です。
実用的な入力に対する信頼性を見ることを目的とします。

#### 構造と入力の選定

MAGNDATAのmcifをLocalへ直接入力するのではなく、同じparent構造の通常CIFをMaterials Project、次いでCrystallography Open Database（COD）から取得しました。

親CIFの対応は、MAGNDATAのparent SGを必須とし、次のような一意性が確認できた場合だけを採用しました。

- MAGNDATA/mcifの化学式と、元素順を無視して一致する候補が1件だけ: MP 495件
- MAGNDATA top pageの化学式と一意に一致: MP 153件
- 記録されたICSD IDにより一意に対応: MP 344件
- 化学式、DOI、citation、元素集合をparent SGと組み合わせて一意に対応: COD 59件

以上で **1,051件**のparent CIFを対応できました。

入力生成では、mcifから直接読める情報だけを使います。

- primary IRをSourceの一意なIR labelへ対応させ、IRに属するKをslotへ入れる
- child magnetic space groupはmcif記載のBNS番号をそのまま使う
- ある元素のsiteが1つでも磁気momentを持てば、その元素全体をmagnetic onにする
- primary IRにordinary成分があれば全元素をdisplacive onにし、strainはoffにする

primary IRやBNSの欠落、Source labelの非一意性、または同一のinput・OPDとしてWebへ投入できない不正なmcif記録を除き、最終的な母集団は**951件**です。
内訳は、one-K 862、two-K 88、three-K 1 です。

| 項目 | 件数 |
|---|---:|
| All | 951 |
| Identical | 914 |
| Equivalent | 951 |
| Unproven | 0 |

- **Equivalent: 951 / 951 = 100.00%**
- **Identical: 914 / 951 = 96.11%**


#### K-signature 別


| K | 母数 | Equivalent | Identical |
|---|---:|---:|---:|
| F | 821 | 821 (100.00%) | 788 (95.98%) |
| FF | 84 | 84 (100.00%) | 81 (96.43%) |
| P | 41 | 41 (100.00%) | 40 (97.56%) |
| PP | 3 | 3 (100.00%) | 3 (100.00%) |
| FP | 1 | 1 (100.00%) | 1 (100.00%) |
| FPP | 1 | 1 (100.00%) | 1 (100.00%) |

## 再現性

`branches/`には、各branchのcanonical input、選択OPD、比較結果、数学検査statusを収録しています。
`verification/`には上記結果を再現する検証ツールが含まれます。
