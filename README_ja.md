# Distortropy

[English](README.md) | [CLIガイド](CLI_ja.md) | [検証](Validation.ja.md)

**Distortropy**は、結晶歪みの対称性解析をlocalで行うPython packageです。同封された
BYU ISOTROPY tableから、k vector、既約表現、秩序変数方向、isotropy subgroup、
invariant basis、symmetry-adapted modeを計算します。通常の計算ではBYU Web serviceや
外部executableを呼びません。

<p align="center">
  <img src="media/distortropy.png" alt="Distortropy graphical interface" width="70%">
</p>

## 入出力と対象範囲

- parent入力: CIF、またはspace groupとordered Wyckoff site
- distortion入力: displacive site、magnetic site、homogeneous strain、1から4個の
  ordered k/irrep factor
- forward pipeline: parent情報、k point、発現irrep、OPD、subgroup、invariant basis、
  atomic/strain mode definition
- embedding query: 指定したparent/subgroupのbasisとoriginに両立するordinaryまたは
  time-odd OP direction
- 出力: compact JSON、optionalなfull pipeline state、text mode table、保存JSON case
- interface: `distortropy` commandとlocal browser interface

CLIとGUIは同じbackend serviceへの入口です。workflowと表示粒度は異なりますが、
科学計算の経路は共通です。

## Install

repositoryをcloneし、packageをinstallします。

```bash
python -m pip install .
```

開発時は`python -m pip install -e .`を使います。どちらも`distortropy` commandを
installし、必要な群論tableもpackageに含まれます。

## Quick Start

各段が返すidentifierを次段の入力に使います。

```bash
# parentと選択可能なcrystallographic siteを確認
distortropy info structure.cif

# k pointとirrepを探索
distortropy kpoints structure.cif
distortropy irreps structure.cif --k R --displacive O

# OPDを列挙し、返されたexact labelを1つ選択
distortropy opds structure.cif --k R --irrep R4- --displacive O
distortropy modes structure.cif \
  --k R --irrep R4- --displacive O --opd P1

# 選択domain上のLandau invariant basisを計算
distortropy invariants structure.cif \
  --k R --irrep R4- --displacive O --opd P1 \
  --minimum-degree 2 --maximum-degree 6
```

具体的なCIFを必要としない場合はspace-group/Wyckoff経路を使います。

```bash
distortropy opds \
  --sg 221 --wyckoff 1a 1b 3c \
  --k R --irrep R4- --displacive c
```

自由なWyckoff座標はOPD・invariant段まで未指定のまま保持できます。atomic modeの
geometryを構成する`modes`段でのみ具体値が必要です。

[CLIガイド](CLI_ja.md)にcommand workflow、入力規約、coupled selection、parametric k、
saved case、exact embeddingの再利用、secondary invariant factor、machine-readable出力を
まとめています。

## Graphical Interface

```bash
distortropy serve --host 127.0.0.1 --port 8300 --open-browser
```

GUIは通常のCIF workflowに加え、space groupとWyckoff siteを直接指定するsymbolic
workflowを備えています。

## Output

pipeline commandは既定でcompact JSONを標準出力します。`modes`は
`--format text`にも対応し、利用可能なcommandでは`--full-state`で完全な内部stateを
取得できます。`-o PATH`を指定した場合だけfileを作成します。

## Validation

検証方法と結果は[Validation.ja.md](Validation.ja.md)を参照してください。

## Notice

元のISODISTORT softwareと群論tableはHarold T. Stokes、Dorian M. Hatch、
Branton J. Campbell（Brigham Young University）の成果です。研究で使用する際は
謝辞をお願いします。詳細は[NOTICE](NOTICE)を参照してください。
