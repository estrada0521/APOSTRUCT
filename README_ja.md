<div align="center">

# APOSTRUCT

Local implementation of ISODISTORT

[CLIガイド](CLI_ja.md) · [検証](Validation.ja.md) · [English](README.md)

</div>

<p align="center">
  <img src="media/APOSTRUCT.png" width="100%" alt="APOSTRUCT graphical interface">
</p>

---

## Overview

親構造のCIFと、変位・磁気・strainの選択、k点と既約表現から、OPD一覧と、選択OPDの部分群・子構造・mode詳細を計算する。

## Setup

C compilerが必要（macOS・Linux）。

```sh
sh build.sh
ln -s "$PWD/bin/apo" ~/.local/bin/apo
```

`bin/apo` と、tableをまとめた `bin/source.blob` が作られる。

## CLI

SrTiO3（Sr(1a), Ti(1b), O(3c)）の例。

```sh
apo info structure.cif
apo kpoints structure.cif
apo irreps structure.cif --k R --displacive O
apo opds structure.cif --k R --irrep R5- --displacive O
apo modes structure.cif --k R --irrep R5- --displacive O --opd P1
```

詳細は[CLIガイド](CLI_ja.md)。

## GUI

```sh
apo serve --open-browser
```

`http://127.0.0.1:8300` で開く。

## Validation

Web版ISODISTORTの出力との比較と、Webから独立した数学検査で検証。詳細は[Validation.ja.md](Validation.ja.md)。

## Layout

| Path | 内容 |
| --- | --- |
| `core/` | 計算本体(C言語) |
| `cli/` | command line |
| `gui/` | `apo serve` とbrowser画面 |
| `source/` | BYU ISOTROPYのtableと、blob生成器 |
| `verification/` | 検証器（Python） |
| `branches/` `assets/` | 判定記録とそのCIF |

## Notice

元のISODISTORTと群論tableは、Harold T. Stokes、Dorian M. Hatch、Branton J. Campbell（Brigham Young University）によるものです。研究で使う場合は謝辞をお願いします。詳細は[NOTICE](NOTICE)。
