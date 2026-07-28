# Distortropy の検証状況

Distortropy は、公開 Web 版 ISODISTORT の出力との比較と、Web から独立した数学的検査の両方を用いて検証しています。
本書の数値は **2026 年 7 月 27 日時点**のスナップショットです。

## 判定

検証単位は、1つの入力構造と1つの order-parameter direction（OPD）の組です。

- **Strict pass**: OPD の表示 setting、子構造、mode の種類・本数・原子対応・正規化係数・vector が、許容された数値精度の範囲で Web 出力と完全に一致します。
- **Physical pass**: 表示上の mode 順序、符号、同じ mode 空間内の基底選択、等価な crystallographic setting などに違いがあっても、子構造と張られる物理的 mode 空間が等価であり、規約差のみを含むものです。
  Strict pass はすべて Physical pass に含まれます。
- **Physical fail**: OPD が存在しない、子構造が物理的に異なる、必要な mode が欠ける、または mode vector が同じ物理空間を張らない、本家とは物理的に異なる場合です（後述するように、本家の一部の出力は数学的要請を満たさない場合があるため、厳密には少し異なります）。

## 二つの検証母集団

### MP pool: 広域ストレス母集団

Materials Project 由来の構造を基礎に、固定・parametric K、multi-K、一般 OPD を広く含む母集団です。

#### 構造 pool の選定

CIF pool は Materials Project から作成しました。
同じ空間群で、占有される Wyckoff site の multiplicity・letter の組が同じ構造は同じ topology とみなし、`(space group, occupied Wyckoff pattern)` ごとに1件だけを取得。
元素の置換だけが異なる構造や、同じ Wyckoff pattern の自由座標だけが異なる構造は、この段階では重複として除きました。
この方法で得た pool は **19,086 CIF** です。
以下の Web 母集団で実際に使われたのは、現時点で2,665 CIF です。

#### 入力のサンプリング

現在の母集団は、1回の単純 random sampling ではなく、複数の再現可能な campaign の和集合です。
各 campaign は seed、filter、目標件数、生成された順序付き input ID を保存しています。

1. CIF を空間群、結晶系、centering、site数、atom数などでfilterします。
2. campaign に応じて、空間群・結晶系・centering の bucket を先に一様選択し、その中の CIF を一様選択します。
   比較用に、bucketを使わず全候補から選ぶ campaign も含みます。
3. strain と、各元素の displacive/magnetic 有効・無効を選びます。
   random指定時は各元素を独立に選びますが、ordinary-only、magnetic-only、mixed などを固定した層も含みます。
4. K slot数を1から4の範囲で選び、fixed/parametric の個数と配置を選びます。
   初期母集団で薄かったmulti-parametric、parametric の位置違い、多次元parametric K、multi-Kは後続campaignで明示的に補いました。
5. 各slotについて、その空間群、選択元素のWyckoff rows、mode kindに対してSource上有効なK/IRだけを列挙し、その候補から選びます。
6. parametric Kの値はSource domain内かつ分母6以下の既約有理数を列挙し、canonical値を重複除去した集合から選びます。
7. WebでOPD一覧を取得した後、campaignが明示した1-based ordinalを収集します。
   1 inputから複数OPDを収集した場合は、それぞれを別branchとして数えます。

下記のpass率は、randomサンプリングを中心に、不足している難しい形状を意図的に厚くした**層別ストレス試験の達成率**です。

| 項目 | 件数 |
|---|---:|
| Web 母集団 | 3,206 branches |
| 入力 / CIF | 3,094 / 2,665 |
| Validation 済み | 3,195 |
| Judged | 3,148 |
| Strict pass | 1,142 |
| Physical-only pass | 1,807 |
| Physical fail | 199 |
| Unjudged | 47 |
| 未検証 | 11 |

- **Physical pass: 2,949 / 3,148 = 93.68%**
- **Strict pass: 1,142 / 3,148 = 36.28%**
- Unjudged も失敗として保守的に数える場合の Physical pass は **2,949 / 3,195 = 92.30%** です。

#### K-signature 別

`F` は固定 K、`P` は parametric K を表します。

| K構成（順序不問） | 母数 | Judged | Physical pass | Strict pass | 未判定 |
|---|---:|---:|---:|---:|---:|
| F | 788 | 788 | 782 (99.24%) | 512 (64.97%) | 0 |
| P | 761 | 748 | 685 (91.58%) | 212 (28.34%) | 13 |
| FF | 297 | 294 | 293 (99.66%) | 186 (63.27%) | 3 |
| FP | 141 | 140 | 123 (87.86%) | 32 (22.86%) | 1 |
| PP | 24 | 22 | 16 (72.73%) | 2 (9.09%) | 2 |
| FFF | 189 | 185 | 185 (100.00%) | 66 (35.68%) | 4 |
| FFP | 490 | 477 | 423 (88.68%) | 50 (10.48%) | 13 |
| FPP | 40 | 36 | 21 (58.33%) | 1 (2.78%) | 4 |
| FFFF | 202 | 200 | 200 (100.00%) | 70 (35.00%) | 2 |
| FFFP | 273 | 258 | 221 (85.66%) | 11 (4.26%) | 15 |

### MAGNDATA: practical 母集団

MAGNDATA に収録された実在磁気構造を基にした母集団です。
実用的な入力に対する信頼性を見ることを目的とします。

#### 構造と入力の選定

MAGNDATAのmcifをLocalへ直接入力するのではなく、同じparent構造の通常CIFをMaterials Project、次いでCrystallography Open Database（COD）から取得しました。

親CIFの対応は、MAGNDATAのparent SGを必須とし、次のような一意性が確認できた場合だけを採用しました。
曖昧な候補から任意の1件を選んではいません。

- MAGNDATA/mcifの化学式と、元素順を無視して一致する候補が1件だけ: MP 495件
- MAGNDATA top pageの化学式と一意に一致: MP 153件
- 記録されたICSD IDにより一意に対応: MP 344件
- 化学式、DOI、citation、元素集合をparent SGと組み合わせて一意に対応: COD 59件

以上で **1,051件**のparent CIFを対応できました。

入力生成では、mcifから推測ではなく直接読める情報だけを使います。

- primary IRをSourceの一意なIR labelへ対応させ、IRに属するKをslotへ入れる
- child magnetic space groupはmcif記載のBNS番号をそのまま使う
- ある元素のsiteが1つでも磁気momentを持てば、その元素全体をmagnetic onにする
- primary IRにordinary成分があれば全元素をdisplacive onにし、strainはoffにする

primary IRの欠落・Source labelの非一意性・BNS欠落・Web側との不整合などを除き、最終的なValidation母集団は**927件**です。
内訳は、zero-K 458、nonzero Type I/III 28、nonzero Type IV 412、two-K 19、three-or-more-K 10 です。

| 項目 | 件数 |
|---|---:|
| 母集団 / 入力 / CIF | 927 / 927 / 927 |
| Validation 済み・Judged | 927 |
| Strict pass | 783 |
| Physical-only pass | 144 |
| Physical fail | 0 |
| Unjudged | 0 |

- **Physical pass: 927 / 927 = 100.00%**
- **Strict pass: 783 / 927 = 84.47%**

#### K-signature 別

MP poolと同様に、slot順ではなくF/Pの個数で集約しています。

| K構成（順序不問） | 母数 | Physical pass | Strict pass |
|---|---:|---:|---:|
| F | 798 | 798 (100.00%) | 679 (85.09%) |
| P | 40 | 40 (100.00%) | 30 (75.00%) |
| FF | 84 | 84 (100.00%) | 69 (82.14%) |
| FP | 1 | 1 (100.00%) | 1 (100.00%) |
| PP | 3 | 3 (100.00%) | 3 (100.00%) |
| FPP | 1 | 1 (100.00%) | 1 (100.00%) |

ストレス母集団と比較して明確に良い結果ですが、これは同一K-signatureでも、ストレス母集団には重い複雑なケースを意図的に多数含めているためです。

## Web 比較と数学的検査

Web 一致率は互換性を測る重要な指標ですが、Web を唯一の真理とは扱っていません。
Distortropy の検証では、Web と Local の両方を対象に、表示された mode basis の exact rational rank など、相手方の出力を正解として仮定しない必要条件も検査します。

外部の blind A/B 監査では、Web と Local のどちらにも、表示された部分群や mode basis の数学的必要条件を満たさない出力が存在することが確認されました。
このため、残ったlong-tailに対しては、次の二軸を分けて扱う予定です。

1. **Web compatibility**: Web と Strict または Physical に一致するか。
2. **Independent mathematics**: 群作用に対する不変性や basis の独立性などを満たすか。

## 数値の再現性と注意事項

このスナップショットの完了済み authority は次の通りです。

- MP pool: runtime `6b9aa754cf70`、comparator `ae1ad1486dbb`
- MAGNDATA: runtime `000ce2ae1f98`、comparator `validation.v7+e276b8fbd237`

また、計算時間制限、Web 側の取得失敗、Local 側の timeout、malformed input は、科学的な Physical/Strict 判定とは別の運用上の失敗として記録します。
特に重い multi-K・parametric ケースでは、結果の正しさだけでなく計算時間にも long tail があります。
