# 検証

## 母集団

### MP pool

Materials Projectから取得した一般的なcifを使った広域のストレス母集団

### MAGNDATA

MAGNDATAに記録された磁気構造を再現する実用母集団

## 数学検査
[`verify_local_output.py`](verification/verify_local_output.py)・[`verify_web_text.py`](verification/verify_web_text.py)

APOSTRUCTとISODISTORTの各出力について、次を検査します。

- 表示されたモード基底の有理数上の線形独立性
- 表示されたモード場の報告部分群に対する不変性

## 出力比較
[`verify_comparison.py`](verification/verify_comparison.py)

APOSTRUCTとISODISTORTの出力を比較します。

- **Identical**：表示設定、子構造、モードの種類・本数・原子対応・正規化係数・ベクトルが、順序等の規約も含めて一致。
- **Equivalent**：設定やモード基底の変換によって、子構造とモード空間の同値を確認。名前（空間群、既約表現、k点、Wyckoff文字、OPDのラベル）は一致が必要です。Identicalを含みます。
- **Unproven**：同値を確認できなかったもの。

APOSTRUCTが目指すのはIdenticalではなくEquivalentです。名前はISODISTORTに合わせますが、同値なものの中での選択（どのdomainを取るか、標準設定の中での基底と原点、モードの符号と基底）はAPOSTRUCT自身の短い規則に従い、異なることがあります。ISODISTORTの出力が数学検査に反する場合は、それに合わせません。該当する枝は記録の`web_math`でわかります。

実行失敗は比較判定と分けて記録します。

## 記録

枝（入力とOPDの組）ごとの計算・検査・比較の記録は[`branches/`](branches/)にあります。

- [`inputs.jsonl.gz`](branches/inputs.jsonl.gz)：入力
- [`results.jsonl.gz`](branches/results.jsonl.gz)：枝ごとの判定と数学検査の状態
- [`manifest.json`](branches/manifest.json)：枝数・入力数・ファイルのハッシュ・検証器のバージョン
