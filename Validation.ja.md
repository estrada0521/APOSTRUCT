# 検証状況

公開 Web 版 ISODISTORT の出力との比較と、Web から独立した数学的検査を用いて検証しています。
本書の数値は **2026 年 7 月 30 日時点**のスナップショットです。
各 Web branch について main へ land 済みの local runtime が生成した最新の Validation 出力を参照しています。

## 判定

検証単位は、1つの入力構造と1つの order-parameter direction（OPD）の組です。

- **Strict pass**: OPD の表示 setting、子構造、mode の種類・本数・原子対応・正規化係数・vector が、順序まで含め Web 出力と完全に一致します。
- **Physical pass**: 表示上の mode 順序、符号、同じ mode 空間内の基底選択、等価な crystallographic setting などに違いがあっても、子構造と張られる物理的 mode 空間が等価であり、基底の取り方など、純粋な規約差のみを含むものです。
  Strict pass はすべて Physical pass に含まれます。
- **Physical fail**: OPD が存在しない、子構造が物理的に異なる、必要な mode が欠ける、または mode vector が同じ物理空間を張らない、本家とは物理的に異なる場合です。本家・ローカル共に数学違反がなく、非自明な物理同値を判定しきれていない潜在的な candidate Physical Pass も保守的に fail として扱います。

## 有効母集団の定義

各出力に二つの必要条件検査を行います。

- **Mode basis:** 表示された完全なmode basisの有理数上の線形独立性
- **Group invariance:** verifierが保持するsettingにおける、表示mode fieldの報告部分群作用に対する不変性

Validation母集団は、WebとLocalがともに計算済みで、Web出力がいずれの検査でもrefutedされなかったbranchで構成します。

## 二つの検証母集団

### MP pool: 広域ストレス母集団

Materials Project 由来の構造を基礎に、固定・parametric K、multi-K、一般 OPD を広く含む母集団です。

#### 構造 pool の選定

CIF pool は Materials Project から作成しました。
同じ空間群で、占有される Wyckoff site の multiplicity・letter の組が同じ構造は同じ topology とみなし、`(space group, occupied Wyckoff pattern)` ごとに1件だけを取得。
元素の置換だけが異なる構造や、同じ Wyckoff pattern の自由座標だけが異なる構造は、この段階では重複として除きました。
この方法で得た pool は **19,086 CIF** です。
以下の Web 母集団で実際に使われたのは、現時点で**2,665 CIF** です。

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

#### 結果

| 項目 | 件数 |
|---|---:|
| 有効Validation母集団 | 3,096 |
| Strict pass | 1,421 |
| Physical pass | 3082 |
| Physical fail | 14 |

- **Physical pass: 3,082 / 3,096 = 99.55%**
- **Strict pass: 1,421 / 3,096 = 45.90%**

#### K-signature 別

`F` は固定 K、`P` は parametric K を表します。
すべての列で、上で定義したValidation母集団を使います。

| K構成（順序不問） | Validation母集団 | Physical pass | Strict pass |
|---|---:|---:|---:|
| F | 788 | 788 (100.00%) | 637 (80.84%) |
| FF | 294 | 294 (100.00%) | 226 (76.87%) |
| FFF | 185 | 185 (100.00%) | 97 (52.43%) |
| FFFF | 200 | 200 (100.00%) | 124 (62.00%) |
| P | 724 | 722 (99.72%) | 214 (29.56%) |
| PP | 20 | 18 (90.00%) | 1 (5.00%) |
| FP | 130 | 130 (100.00%) | 39 (30.00%) |
| FFP | 459 | 457 (99.56%) | 59 (12.85%) |
| FFFP | 261 | 261 (100.00%) | 22 (8.43%) |
| FPP | 34 | 26 (76.47%) | 2 (5.88%) |
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
| Strict pass | 801 |
| Physical pass | 927 |
| Physical fail | 0 |

- **Physical pass: 927 / 927 = 100.00%**
- **Strict pass: 801 / 927 = 86.41%**


#### K-signature 別


| K構成（順序不問） | 母数 | Physical pass | Strict pass |
|---|---:|---:|---:|
| F | 798 | 798 (100.00%) | 697 (87.34%) |
| FF | 84 | 84 (100.00%) | 69 (82.14%) |
| P | 40 | 40 (100.00%) | 30 (75.00%) |
| PP | 3 | 3 (100.00%) | 3 (100.00%) |
| FP | 1 | 1 (100.00%) | 1 (100.00%) |
| FPP | 1 | 1 (100.00%) | 1 (100.00%) |


## 検証限界

Webの数学違反により有効母集団から外れたケースは現時点で**65件**を確認済みです。
対応するLocal出力は、内64件が現行の二検査をsatisfied、1件がrefutedでした。
これはLocalを証明するものではないため、今度のより強い数学検査による要検証対象として扱います。


## 計算性能

有効Validation母集団のうち、WebとLocalの両方に`wall_s`が記録されているbranchを性能母集団とします。
Web時間はCIF uploadからComplete Mode Detailsのtext変換までを測ります。
Local時間はローカルのmode-details計算本体を測ります。実行環境が異なるため、以下は速度の直接比較ではなく参考値です。

| 性能母集団 | Local中央値 | Web中央値 | Local p95 | Web p95 |
|---:|---:|---:|---:|---:|
| 2,114 / 4,023 | 1.23 s | 8.37 s | 52.07 s | 31.84 s |

```text
          0.1                 1                  10                100 s
Local     |--------[==========│=============]-----------------|
        p5=0.13  Q1=0.35  median=1.23  Q3=6.40  p95=52.07
Web                                          [│==]--------|
        p5=7.32  Q1=7.49  median=8.37  Q3=10.69  p95=31.84
```

### K-signature別性能

| K構成（順序不問） | 性能母集団 | Local中央値 | Web中央値 | Local p95 | Web p95 |
|---|---:|---:|---:|---:|---:|
| F | 951 | 0.37 s | 7.48 s | 2.57 s | 7.94 s |
| FF | 163 | 0.79 s | 8.79 s | 4.98 s | 10.67 s |
| FFF | 52 | 1.87 s | 9.24 s | 13.03 s | 29.05 s |
| FFFF | 47 | 2.22 s | 9.30 s | 12.50 s | 20.05 s |
| P | 250 | 3.43 s | 8.49 s | 60.47 s | 29.16 s |
| PP | 23 | 14.73 s | 13.47 s | 98.42 s | 61.71 s |
| FP | 131 | 5.23 s | 10.58 s | 85.47 s | 41.10 s |
| FFP | 243 | 8.80 s | 11.84 s | 73.62 s | 47.34 s |
| FFFP | 218 | 13.71 s | 15.00 s | 90.13 s | 62.02 s |
| FPP | 35 | 14.11 s | 14.02 s | 97.78 s | 53.72 s |
| FFPP | 1 | 96.40 s | 66.64 s | 96.40 s | 66.64 s |

parametric Kとmulti-Kを含むstress caseには、Local に明確な long tail が残っています。
