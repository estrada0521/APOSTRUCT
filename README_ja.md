# Distortropy

[English](README.md)

**Distortropy** は **ISODISTORT** の純 Python・オフライン実装です。

<p align="center">
  <img src="media/distortropy-1.png" alt="Distortropy structure selection" width="49%">
  <img src="media/distortropy-2.png" alt="Distortropy mode viewer" width="49%">
</p>

## Setup

リポジトリをcloneし、依存packageと`isodistort` commandをeditable installします。

```bash
python -m pip install -e .
```

群論テーブル 9 ファイルをリポジトリ直下の `Source/` に置いてください（BYU の ISOTROPY
Software Suite <https://iso.byu.edu/> に含まれます）。

```
const.dat  data_isotropy.txt  data_irreps.txt  data_images.txt  data_little.txt
data_magnetic.txt  data_space.txt  data_ssgmag.txt  data_wyckoff.txt
```

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

## Validation
現在も検証を続けています。
実用的な入力では本家と物理的に非等価な出力をすることは殆どありませんが、完全一致を保証するものではありません。
（`How to Use`における入力例のような）重い非自明なケースでは差が出ることがあるため、ご注意下さい。


## Notice

元になったISODISTORT・群論 table は Harold T. Stokes, Dorian M. Hatch,
Branton J. Campbell（BYU）の成果です。
研究で用いる場合は謝辞をお願いします（[NOTICE](NOTICE)）。
