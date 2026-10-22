# 検証状況

公開 Web 版 ISODISTORT の出力との比較と、Web から独立した数学的検査を用いて検証しています。
本書の数値は **2026 年 10 月 15 日時点**のスナップショットです。
各 Web branch について main へ land 済みの local runtime が生成した最新の Validation 出力を参照しています。

## 判定

検証単位は、1つの入力構造と1つの order-parameter direction（OPD）の組です。

- **Identical**: OPD の表示 setting、子構造、mode の種類・本数・原子対応・正規化係数・vector が、順序まで含め Web 出力と完全に一致します。
- **Equivalent**: 比較器が持つ明示的な結晶学的・mode 空間の同値（mode の順序、符号、同じ mode 空間内の基底、等価な setting）によって、子構造と張られる mode 空間が Web 出力と同値であると証明できたものです。Identical はすべて Equivalent に含まれます。
- **Unproven**: 比較器が同値を証明できなかったものです。物理的に異なる出力のほか、比較器が同値を判定しきれない出力や、Web 出力そのものが数学検査で反証される出力も含みます。ローカル出力が誤りであるという判定ではありません。

## 有効母集団の定義

各出力に二つの必要条件検査を行います。

- **Mode basis:** 表示された完全なmode basisの有理数上の線形独立性
- **Group invariance:** verifierが保持するsettingにおける、表示mode fieldの報告部分群作用に対する不変性

Validation母集団は、WebとLocalがともに計算済みで、Web出力がいずれの検査でもrefutedされなかったbranchで構成します。

## 再現性

公開checkoutの`branches/`には、各branchのcanonical input、選択OPD、比較結果、数学検査statusを収録します。
独立した`verification/`には、private monorepo Validation運用と同じparser、数学検査、Strict/Physical comparator、result projectionを収録します。
第三者はbranchを選び、APOSTRUCTを実行し、同じ入力を公開ISODISTORT Webへ自分で投入することで、privateなsweep/store schemaに依存せず両出力を検査できます。

## 二つの検証母集団

### MP pool: 広域ストレス母集団

Materials Project 由来の構造を基礎に、固定・parametric K、multi-K、一般 OPD を広く含む母集団です。

#### 構造 pool の選定

CIF pool は Materials Project から作成しました。
同じ空間群で、占有される Wyckoff site の multiplicity・letter の組が同じ構造は同じ topology とみなし、`(space group, occupied Wyckoff pattern)` ごとに1件だけを取得。
元素の置換だけが異なる構造や、同じ Wyckoff pattern の自由座標だけが異なる構造は、この段階では重複として除きました。
この方法で得た pool は **19,086 CIF** です。
以下の Web 母集団で実際に使われたのは、現時点で**2,557 CIF** です。

#### 入力のサンプリング

現在の母集団は、1回の単純 random sampling ではなく、複数の再現可能な campaign の和集合です。
各 campaign は seed、filter、目標件数、生成された順序付き input ID を保存しています。

1. CIF を空間群、結晶系、centering、site数、atom数などでfilterします。
2. campaign に応じて、空間群・結晶系・centering の bucket を先に一様選択し、その中の CIF を一様選択します。比較用に、bucketを使わず全候補から選ぶ campaign も含みます。
3. strain と、各元素の displacive/magnetic 有効・無効を選びます。random指定時は各元素を独立に選びますが、ordinary-only、magnetic-only、mixed などを固定した層も含みます。
4. K slot数を1から4の範囲で選び、fixed/parametric の個数と配置を選びます。初期母集団で薄かったmulti-parametric、parametric の位置違い、多次元parametric K、multi-Kは後続campaignで明示的に補いました。
5. 各slotについて、その空間群、選択元素のWyckoff rows、mode kindに対してSource上有効なK/IRだけを列挙し、その候補から選びます。
6. parametric Kの値はSource domain内かつ分母6以下の既約有理数を列挙し、canonical値を重複除去した集合から選びます。
7. WebでOPD一覧を取得した後、campaignが明示した1-based ordinalを収集します。1 inputから複数OPDを収集した場合は、それぞれを別branchとして数えます。

下記のpass率は、randomサンプリングを中心に、不足している難しい形状を意図的に厚くした**層別ストレス試験の達成率**です。

#### 結果

| 項目 | 件数 |
|---|---:|
| 有効Validation母集団 | 3,040 |
| Identical | 2,599 |
| Equivalent | 3,035 |
| Unproven | 5 |

- **Equivalent: 3,035 / 3,040 = 99.84%**
- **Identical: 2,599 / 3,040 = 85.49%**

#### K-signature 別

`F` は固定 K、`P` は parametric K を表します。
すべての列で、上で定義したValidation母集団を使います。

| K | Validation母集団 | Equivalent | Identical |
|---|---:|---:|---:|
| F | 788 | 788 (100.00%) | 734 (93.15%) |
| FF | 297 | 297 (100.00%) | 252 (84.85%) |
| FFF | 189 | 189 (100.00%) | 167 (88.36%) |
| FFFF | 202 | 202 (100.00%) | 186 (92.08%) |
| P | 687 | 686 (99.85%) | 599 (87.19%) |
| PP | 20 | 18 (90.00%) | 7 (35.00%) |
| FP | 126 | 126 (100.00%) | 106 (84.13%) |
| FFP | 453 | 453 (100.00%) | 363 (80.13%) |
| FFFP | 246 | 246 (100.00%) | 168 (68.29%) |
| FPP | 31 | 29 (93.55%) | 17 (54.84%) |
| FFPP | 1 | 1 (100.00%) | 0 (0.00%) |

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

primary IRやBNSの欠落、Source labelの非一意性、または同一のinput・OPDとしてWebへ投入できない不正なmcif記録を除き、最終的なValidation母集団は**927件**です。
内訳は、zero-K 458、nonzero Type I/III 28、nonzero Type IV 412、two-K 19、three-or-more-K 10 です。

| 項目 | 件数 |
|---|---:|
| 母集団 / 入力 / CIF | 927 / 927 / 927 |
| Validation 済み・Judged | 927 |
| Identical | 892 |
| Equivalent | 927 |
| Unproven | 0 |

- **Equivalent: 927 / 927 = 100.00%**
- **Identical: 892 / 927 = 96.22%**


#### K-signature 別


| K | 母数 | Equivalent | Identical |
|---|---:|---:|---:|
| F | 798 | 798 (100.00%) | 767 (96.12%) |
| FF | 84 | 84 (100.00%) | 81 (96.43%) |
| P | 40 | 40 (100.00%) | 39 (97.50%) |
| PP | 3 | 3 (100.00%) | 3 (100.00%) |
| FP | 1 | 1 (100.00%) | 1 (100.00%) |
| FPP | 1 | 1 (100.00%) | 1 (100.00%) |


## 検証限界

successful Web branch 4,141件のうち、Local比較が記録されたbranchは4,139件です（残る2件はLocal側のoperational failure）。
このうち現行数学検査で**171件**がrefuted、**1件**がindeterminateです。
この172件はすべて有効母集団から除外します。
除外されたWeb出力との一致はcorrectness evidenceには数えません。


## 計算性能

有効Validation母集団のうち、WebとLocalの両方に`wall_s`が記録されているbranchを性能母集団とします。
Web時間は、successful OPD一覧requestと、選択branchのComplete Mode Details requestの和です。
Local時間は、ローカルの実行ファイル`apo`を1 branchにつき1プロセスで呼んだwall timeです（起動を含み、6並列のsweep中に計測）。
実行環境が異なるため、以下は速度の直接比較ではなく参考値です。

| 性能母集団 | Local中央値 | Web中央値 | Local p95 | Web p95 |
|---:|---:|---:|---:|---:|
| 2,081 / 3,967 | 0.03 s | 8.27 s | 1.11 s | 30.89 s |

```text
          0.01        0.1         1           10          100 s
Local         [=│=====]------------|
        p5=0.02  Q1=0.02  median=0.03  Q3=0.11  p95=1.11
Web                                         [│]-----|
        p5=7.31  Q1=7.48  median=8.27  Q3=10.44  p95=30.89
```

### K-signature別性能

| K | 性能母集団 | Local中央値 | Web中央値 | Local p95 | Web p95 |
|---|---:|---:|---:|---:|---:|
| F | 951 | 0.02 s | 7.48 s | 0.04 s | 7.94 s |
| FF | 164 | 0.02 s | 8.80 s | 0.09 s | 10.68 s |
| FFF | 53 | 0.04 s | 9.24 s | 0.31 s | 37.38 s |
| FFFF | 47 | 0.05 s | 9.30 s | 0.30 s | 20.05 s |
| P | 235 | 0.06 s | 8.48 s | 1.39 s | 32.51 s |
| PP | 23 | 0.25 s | 13.47 s | 1.70 s | 37.34 s |
| FP | 127 | 0.11 s | 10.58 s | 1.82 s | 41.10 s |
| FFP | 240 | 0.17 s | 11.54 s | 1.97 s | 47.92 s |
| FFFP | 208 | 0.24 s | 14.03 s | 2.42 s | 50.26 s |
| FPP | 32 | 0.19 s | 13.41 s | 1.94 s | 53.72 s |
| FFPP | 1 | 1.67 s | 66.64 s | 1.67 s | 66.64 s |

parametric Kとmulti-Kを含むstress caseでも、Localのp95は2.5秒以内です。
