# CLI

[English](CLI.md) · [README](README_ja.md)

各commandはstateless。前段のJSONが返すlabelを、次段に渡す。

```text
info -> kpoints -> irreps -> opds -> modes
```

## Commands

| Command | 結果 |
| --- | --- |
| `info` | 正規化したparent cellと、選択できるsite |
| `kpoints` | k点のlabel、kvector、star、parameter名 |
| `irreps` | 選んだsiteとstrainで発現するirrep |
| `opds` | OPD一覧と、各OPDの部分群・basis・origin・domain |
| `modes` | 1つのOPDのmode詳細 |
| `run` | JSONのcaseを指定段まで実行 |
| `serve` | browser画面を起動 |
| `show` | 起動中の画面にcaseを読み込む |

## Parent

CIFのpathを位置引数で渡す。setting、site、占有率、座標はCIFに従う。
座標は小数、誤差付き表記（`0.1234(5)`）、`1/2` のような有理数で書ける。
Wyckoff位置を決められないsiteは `info.sites[].wyckoff_mapping_error` に理由が入り、そのsiteを選ぶと停止する。

## Sites

`--displacive` と `--magnetic` は `info.sites` の `label` か `type` を取る。`label` は1つのsite、`type` はその種類の全siteを選ぶ。
`--strain` はhomogeneous strainを加える。strainだけの `irreps` では `--k` を省略でき、`GM` とみなす。

```sh
apo irreps structure.cif --k GM --displacive O --magnetic Fe
```

## k points

labelとparameter名は `kpoints` の `label` と `parameter_names` に従う。parameterは順に並べるか、名前で指定する。

```sh
apo irreps structure.cif --k DT 1/3 --displacive O
apo irreps structure.cif --k DT b=1/3 --displacive O
```

`miller_love_kvector` は出力の表示用で、入力には使わない。
parametric kがspecial pointに一致する値は拒否する。その場合はfixed pointのlabelを使う。

`--k` と `--irrep` は1から4組まで並べられ、書いた順に対応する。

```sh
apo opds structure.cif --k R --irrep R4- --k M --irrep M3+ --displacive O
```

## OPDs

`--opd` には `opds` が返した `label` をそのまま渡す。directionから組み立て直さない。
複数のkでは、labelにslotごとのdomain番号が付く（`P1(1)P1(2)`）。shellではquoteする。

```sh
apo modes structure.cif --k R --irrep R4- --k M --irrep M3+ --displacive O --opd 'P1(1)P1(2)'
```

`index` は群の指数、`cell_index` はsupercellの倍率。`ferroic_classified` が `false` のときは、`ferroic_properties` が空でも未分類を意味する。

## Modes

各modeは `definition_id`（そのpayload内でのID）と、kind・k点・irrep・direction・site・Wyckoff位置・site irrepからなるidentityを持つ。`role` はprimaryと誘起されたsecondaryを分ける。
`--format text` はmode表を出す（`modes` と `run --upto modes` のみ）。

## Case

選択をJSONにまとめたもの。

```json
{
  "structure": "structure.cif",
  "sites": {"O": {"displacive": true, "magnetic": false}},
  "strain": false,
  "k": [{"label": "R", "ir": "R4-"}],
  "opd": "P1"
}
```

- fieldは `structure`・`cif`・`sites`・`strain`・`k`・`opd` のみ。`structure` と `cif` はどちらか一方。
- `structure` の相対pathは、case fileのdirectoryから解決する。
- `cif` は12文字のcontent IDで、このcheckoutの `assets/cif/` か `assets/magndata/cif/` から探す。
- `k` の各要素は `label`、`params`（名前と値のobject）、`ir` を持つ。`ir` はOPD段から必要。
- fileは拡張子 `.json` が必要。`--case -` で標準入力から読む。

```sh
apo run --case case.json --upto opds
apo modes --case case.json
```

`run` は一つのprocessで全段を計算する。自動処理ではこちらを使う。`--case` は個別の選択optionと併用できない。

## Output

出力はJSONで、`schema`（`APOSTRUCT.cli.info` など）で段を示す。

| Option | 内容 |
| --- | --- |
| `-o PATH` | 標準出力の代わりにfileへ書く |
| `--indent N` | JSONの字下げ |
| `--full-state` | 計算stateをすべて出す |

errorは標準errorに出て、status 2で終わる。見つからないk、irrep、OPDを別のものに置き換えて続けることはない。

## GUI

```sh
apo serve --host 127.0.0.1 --port 8300 --open-browser
apo show --case case.json --server http://127.0.0.1:8300
```

`show` はcaseの最も深い段まで開く。kだけならirrep選択、irrepまでならOPD選択、`opd` まであればmode詳細と構造表示。
`--amplitude DEFINITION_ID=VALUE` でmode sliderの初期値を指定する（`--amplitude magnetic-1=1`）。
CIFの本文は `show` 側で読んでserverに送る。caseと計算結果はserverのmemoryにだけ残る。
