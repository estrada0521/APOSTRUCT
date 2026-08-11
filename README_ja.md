# APOSTRUCT

[English](README.md) | [CLIガイド](CLI_ja.md) | [検証](Validation.ja.md)

**APOSTRUCT** は、結晶の対称性解析をローカルで行う Python パッケージです。
同梱の BYU ISOTROPY テーブルから、k ベクトル、既約表現、秩序変数方向、isotropy 部分群、不変基底、対称性適合モード、Landau不変量を計算します。

<p align="center">
  <img src="media/APOSTRUCT.png" alt="APOSTRUCT graphical interface" width="100%">
</p>

名前について — apo- は「離れて」を意味するギリシャ語の接頭辞で、apomorphy（祖先形質から派生した形質）と同じ語法にあたります。
ISODISTORT の iso- が秩序変数を不変に保つ部分群、すなわち対称性の残る側を名指すのに対し、apo- は親構造から離れる側を名指します。
struct はその対象であり、変位・磁気・歪みを含む構造そのものを指します。

## 入出力と対象範囲

- parent入力: CIF、またはspace groupとordered Wyckoff site
- distortion入力: displacive site、magnetic site、homogeneous strain、1から4個の
  ordered k/irrep factor
- forward pipeline: parent情報、k point、発現irrep、OPD、subgroup、invariant basis、
  atomic/strain mode definition
- embedding query: 指定したparent/subgroupのbasisとoriginに両立するordinaryまたは
  time-odd OP direction
- 出力: compact JSON、optionalなfull pipeline state、text mode table、保存JSON case
- interface: `apo` commandとlocal browser interface

CLIとGUIは同じbackend serviceへの入口です。workflowと表示粒度は異なりますが、科学計算の経路は共通です。

## Install

repositoryをcloneし、packageをinstallします。

```bash
python -m pip install .
```

開発時は`python -m pip install -e .`を使います。どちらも`apo` commandをinstallし、必要な群論tableもpackageに含まれます。

## Quick Start

以下はSrを1a、Tiを1b、Oを3cに置いたSrTiO3 parentの例です。各段が返すidentifierを次段の入力に使います。

```bash
# parentと選択可能なcrystallographic siteを確認
apo info structure.cif

# k pointとirrepを探索
apo kpoints structure.cif
apo irreps structure.cif --k R --displacive O

# OPDを列挙し、返されたexact labelを1つ選択
apo opds structure.cif --k R --irrep R5- --displacive O
apo modes structure.cif \
  --k R --irrep R5- --displacive O --opd P1

# 選択domain上のLandau invariant basisを計算
apo invariants structure.cif \
  --k R --irrep R5- --displacive O --opd P1 \
  --minimum-degree 2 --maximum-degree 6
```

具体的なCIFを必要としない場合はspace-group/Wyckoff経路を使います。

```bash
apo opds \
  --sg 221 --wyckoff 1a 1b 3c \
  --k R --irrep R5- --displacive c
```

この原点選択では反位相八面体回転は`R5-`です。parent原点を
`(1/2,1/2,1/2)`だけ移し、Tiを1a、Oを3dに置く規約では、同じ物理的回転が
Howard-Stokesの慣用記号`R4+`になります。

自由なWyckoff座標はOPD・invariant段まで未指定のまま保持できます。atomic modeの
geometryを構成する`modes`段でのみ具体値が必要です。

[CLIガイド](CLI_ja.md)にcommand workflow、入力規約、coupled selection、parametric k、
saved case、exact embeddingの再利用、secondary invariant factor、machine-readable出力を
まとめています。

## Graphical Interface

```bash
apo serve --host 127.0.0.1 --port 8300 --open-browser
```

既存caseは`apo show --case case.json`で稼働中GUIへ読み込めます。

GUIは通常のCIF workflowに加え、space groupとWyckoff siteを直接指定するsymbolic
workflowを備えています。

## Output

pipeline commandは既定でcompact JSONを標準出力します。`modes`は
`--format text`にも対応し、利用可能なcommandでは`--full-state`で完全な計算stateを
取得できます。`-o PATH`を指定した場合だけfileを作成します。

## Validation

検証方法と結果は[Validation.ja.md](Validation.ja.md)を参照してください。

## Notice

元のISODISTORT softwareと群論tableはHarold T. Stokes、Dorian M. Hatch、
Branton J. Campbell（Brigham Young University）の成果です。研究で使用する際は
謝辞をお願いします。詳細は[NOTICE](NOTICE)を参照してください。
