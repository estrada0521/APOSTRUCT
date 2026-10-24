<div align="center">

# APOSTRUCT

結晶対称性の解析エンジンです
ISODISTORT相当の計算をローカルで実行します

[CLIガイド](CLI_ja.md) · [検証](Validation.ja.md) · [English](README.md)

</div>

<p align="center">
  <img src="https://media.githubusercontent.com/media/estrada0521/APOSTRUCT/b76ddd9/media/APOSTRUCT.png" width="100%" alt="APOSTRUCT graphical interface">
</p>

---

## Overview

親構造のCIFと、変位・磁気・strainの選択、k点と既約表現から、OPD一覧と、選択OPDの部分群・子構造・mode詳細を計算する。ラベルと、数学的に一意に決まらない選択（既約表現の基底、設定、代表の方向など）は、ISOTROPY Software Suite のデータが記録する規約（標準の結晶学の資料と、ISOTROPY 自身の選択）に従う。

## How to use

C compilerが必要（macOS・Linux）。buildすると `bin/apo` と、tableをまとめた `bin/source.blob` が作られる。WindowsはMinGW-w64でbuildする（MSYS2上またはcross buildで `CC=x86_64-w64-mingw32-gcc sh build.sh`）。`bin/apo.exe` は `bin/source.blob` と同じ場所に置く。

```sh
sh build.sh
ln -s "$PWD/bin/apo" ~/.local/bin/apo
```

CLIの例（SrTiO3、Sr(1a), Ti(1b), O(3c)）。詳細は[CLIガイド](CLI_ja.md)。

```sh
apo info structure.cif
apo kpoints structure.cif
apo irreps structure.cif --k R --displacive O
apo opds structure.cif --k R --irrep R5- --displacive O
apo modes structure.cif --k R --irrep R5- --displacive O --opd P1
```

GUIは `http://127.0.0.1:8300` で開く。

```sh
apo serve --open-browser
```

## Layout

| Path | 内容 |
| --- | --- |
| `core/` | 計算エンジン(C言語) |
| `cli/` | command line |
| `gui/` | `apo serve` とbrowser画面 |
| `source/` | 基本的な静的データとBYU規約の表 |
| `verification/` | 検証器（Python）。本家ISODISTORTの出力との比較と、本家から独立した数学検査。結果は[Validation.ja.md](Validation.ja.md) |
| `branches/` `assets/` | 判定記録とそのCIF |

## Notice

APOSTRUCTの出力はHarold T. Stokes、Dorian M. Hatch、Branton J. Campbell（Brigham Young University: 以下BYU）の規約に従います。研究で使う場合は謝辞をお願いします。詳細は[NOTICE](NOTICE)。
