# Distortropy

[English](README.md)

**Distortropy** は **ISODISTORT の計算エンジンを内部実装した、オープンソースの
純 Python 実装**です。BYU の Web serviceを呼び出すwrapperではありません。CIFと
BYU ISOTROPYの群論tableから、k vector・既約表現・秩序変数方向・isotropy subgroupを
導出し、displacive・magnetic・strain mode definitionをDistortropy自身が構成します。
全ての計算はlocalかつofflineで実行されます。

<p align="center">
  <img src="media/distortropy.png" alt="Distortropy graphical interface" width="70%">
</p>

## Setup

リポジトリをcloneし、依存packageと`isodistort` commandをeditable installします。

```bash
python -m pip install -e .
```

必要なBYU ISOTROPY群論tableは`Source/`に同封されています。

## How to Use

```bash
# 1. 親空間群の k 点を列挙
isodistort kpoints structure.cif

# 2. 選択した k 点の既約表現を列挙
isodistort irreps structure.cif \
  --k L \
  --k GP 1/3 1/4 2/5 \
  --k B 1/3 2/5 \
  --k W 1/4

# 3. 選択した既約表現の OPD(秩序変数方向)を列挙
isodistort opds structure.cif \
  --k L --irrep L1 \
  --k GP 1/3 1/4 2/5 --irrep GP1GQ1 \
  --k B 1/3 2/5 --irrep mB2BA2 \
  --k W 1/4 --irrep W1WA1 \
  --displacive Sn Fe --magnetic O Fe --strain

# 4. 選択に対するモード詳細を計算
isodistort modes structure.cif \
  --k L --irrep L1 \
  --k GP 1/3 1/4 2/5 --irrep GP1GQ1 \
  --k B 1/3 2/5 --irrep mB2BA2 \
  --k W 1/4 --irrep W1WA1 \
  --displacive Sn Fe --magnetic O Fe --strain \
  --opd 'P1(1)P3(1)C2(1)P2(1)'
```

ケースを JSON に焼いて入力

```json
{ "structure": "structure.cif",
  "k": [
    { "label": "L", "ir": "L1" },
    { "label": "GP", "params": { "a": "1/3", "b": "1/4", "g": "2/5" },
      "ir": "GP1GQ1" },
    { "label": "B", "params": { "a": "1/3", "g": "2/5" }, "ir": "mB2BA2" },
    { "label": "W", "params": { "g": "1/4" }, "ir": "W1WA1" }
  ],
  "sites": {
    "Sn": { "displacive": true,  "magnetic": false },
    "O":  { "displacive": false, "magnetic": true },
    "Fe": { "displacive": true,  "magnetic": true },
    "Ba": { "displacive": false, "magnetic": false }
  },
  "strain": true,
  "opd": "P1(1)P3(1)C2(1)P2(1)"
}
```

```bash
isodistort modes --case case.json
```

ケースを .in に焼いて入力

```text
CIF structure.cif

K L
IR L1
K GP 1/3 1/4 2/5
IR GP1GQ1
K B 1/3 2/5
IR mB2BA2
K W 1/4
IR W1WA1

STRAIN
DISPLACIVE Sn Fe
MAGNETIC O Fe
OPD P1(1)P3(1)C2(1)P2(1)
```

```bash
isodistort modes --case case.in
```

## Output

`kpoints`、`irreps`、`opds` は JSON を標準出力します。`modes` は JSON caseでは
JSON、`.in`または直接指定ではcomplete mode-details textを標準出力します。
`modes --format json|text`で形式を指定でき、`-o PATH`で指定したファイルへ保存できます。
`-o`を指定しない限りファイルは作成しません。

## Server

ローカルサーバーを起動して、本家と同じ使用感で利用可能です。

```bash
isodistort serve --host 127.0.0.1 --port 8300
```

Invariants計算は現時点ではこのGUIから利用できます。

## Validation

検証方法と結果は[Validation.ja.md](Validation.ja.md)を参照してください。

## Notice

元になったISODISTORT・群論 table は Harold T. Stokes, Dorian M. Hatch,
Branton J. Campbell（BYU）の成果です。
研究で用いる場合は謝辞をお願いします（[NOTICE](NOTICE)）。
