<div align="center">

# APOSTRUCT

結晶対称性の解析エンジンです。
ISODISTORT相当の計算をローカルで実行します。

[マニュアル](Manual.md) · [検証](Validation.jp.md) ・[NOTICE](NOTICE)

</div>

<p align="center">
  <img src="https://media.githubusercontent.com/media/estrada0521/APOSTRUCT/b76ddd9/media/APOSTRUCT.png" width="100%" alt="APOSTRUCT graphical interface">
</p>

---

## Overview

親構造のCIFから、変位・磁気・占有・ひずみに対応する既約表現、秩序変数の方向（OPD）、部分群、子構造、モードを計算します。

## Layout

| Path | 内容 |
| --- | --- |
| `core/` | 計算エンジン（C言語） |
| `cli/` | コマンドライン |
| `gui/` | `apo serve` とブラウザ画面 |
| `source/` | 基本的な静的データと規約の表 |
| `verification/` | 検証器（Python）。ISODISTORTの出力との比較と、独立した数学検査。結果は[Validation.jp.md](Validation.jp.md) |
| `branches/`・`assets/` | 判定記録とそのCIF |

## Setup

- macOS・Linux：Cコンパイラでビルドします。

  ```sh
  sh build.sh
  ln -s "$PWD/bin/apo" ~/.local/bin/apo
  ```

- Windows：MSYS2上のMinGW-w64、またはクロスビルドを使います。

  ```sh
  CC=x86_64-w64-mingw32-gcc sh build.sh
  ```

ビルドすると `bin/apo`（Windowsでは `bin/apo.exe`）と、テーブルをまとめた `bin/source.blob` が生成されます。実行ファイルと `source.blob` は同じディレクトリに置きます。

## Example

SrTiO3(Sr[1a], Ti[1b], O[3c])の例

```sh
apo info structure.cif  # 親構造とサイトを確認
apo kpoints structure.cif  # k点の一覧を取得
apo irreps structure.cif --k R --displacive O  # R点で酸素の変位に対応する既約表現を取得
apo opds structure.cif --k R --irrep R5- --displacive O  # R5-の秩序変数の方向（OPD）と部分群を取得
apo modes structure.cif --k R --irrep R5- --displacive O --opd P1  # P1方向の子構造と変位モードを計算
```


詳細は[Manual](Manual.md)。