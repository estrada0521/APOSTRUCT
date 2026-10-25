# 検証状況

本家ISODISTORT の出力との比較と、独立した数学的検査を用いて検証しています。

## 数学検査

- 表示された完全なmode basisの有理数上の線形独立性
- 表示mode fieldの報告部分群作用に対する不変性

を通過した本家出力を比較に用います。

## 判定の定義

- **Identical**: OPD の表示 setting、子構造、mode の種類・本数・原子対応・正規化係数・vector が、順序等のBYU固有の規約まで含め本家出力と完全に一致します。
- **Equivalent**: 明示的な結晶学的・mode 空間の同値（mode の順序、符号、同じ mode 空間内の基底、等価な setting）によって、子構造と張られる mode 空間が 本家出力と同値であると証明されたものです。Identical はすべて Equivalent に含まれます。
- **Unproven**: 比較器が同値を証明できなかったものです。物理的に異なる出力のほか、比較器が同値を判定しきれない場合です。APOSTRUCTの出力が誤りであるという判定ではありません。

## 検証母集団

### MP pool

Materials Project 由来の構造を基礎に、固定・parametric K、multi-K、一般 OPD を広く含む母集団です。

#### 構造 pool の選定

CIF pool は Materials Project から作成しました。
同じ空間群で、占有される Wyckoff site の multiplicity・letter の組が同じ構造は同じ topology とみなし、`(space group, occupied Wyckoff pattern)` ごとに1件だけを取得。
元素の置換だけが異なる構造や、同じ Wyckoff pattern の自由座標だけが異なる構造は、この段階では重複として除きました。

### MAGNDATA

MAGNDATAに収録された実在磁気構造を基にした母集団です。親空間群の一致を必須とし、化学式・ICSD ID・引用情報等から一意に対応する親CIFをMaterials Project、次いでCODから取得します。

入力はmcifの主既約表現と子磁気空間群のBNS番号から生成します。主既約表現はSourceの一意なラベルに対応づけ、そのk点を選択します。磁気モーメントを持つサイトがある元素全体で磁気モードを有効にし、主既約表現にordinary成分があれば全元素の変位を有効、strainを無効にします。

必要な情報の欠落、Sourceラベルの非一意性、同一の入力・OPDとして本家に投入できない記録は除きます。

## 再現性

[`branches/`](branches/)には、各branchのcanonical input、選択OPD、比較結果、数学検査statusを収録しています。
[`verification/`](verification/)には結果を再現する検証ツールが含まれます。

入力は[`inputs.jsonl.gz`](branches/inputs.jsonl.gz)、判定は[`results.jsonl.gz`](branches/results.jsonl.gz)、件数と検証器のバージョンは[`manifest.json`](branches/manifest.json)を参照してください。
