# APOSTRUCT CLIガイド

[English](CLI.md) | [README](README_ja.md)

CLIはlocal GUIと同じ計算serviceを公開します。
各commandはstatelessです。
前段の JSONからlabelやselectorを選び、それを次のcommandへ渡します。

```text
parent -> k points -> irreps -> OPDs -> modes
```

`apo --version`はinstall済みversionを表示し、Git checkout内ではcommitも併記します。
全optionは`apo COMMAND --help`で確認できます。
この文書では各commandに共通する科学的・machine-readable契約を説明します。

## Command一覧

| Command | 結果 |
|---|---|
| `info` | 正規化されたparent cellと選択可能なcrystallographic site |
| `kpoints` | k-point label、式、star、parameter名 |
| `irreps` | 選択したmode kindで発現するordinary/magnetic irrep |
| `opds` | OPD、subgroup、basis、origin、domain、ferroic metadata |
| `modes` | 1つのOPDに対するatomic/strain mode definition |
| `run` | 保存したJSON caseを指定段まで実行 |
| `serve` | local GUIを起動 |
| `show` | 保存JSON caseを稼働中のGUIへ読み込む |

## Parent入力

CIF pathを位置引数へ渡します。

```bash
apo info structure.cif
apo kpoints structure.cif
```

parent setting、site、occupancy、座標はCIFをauthorityとします。
数値座標には通常の decimal、estimated-standard-deviation付き表記、`1/2`のようなrational表記を使えます。
Wyckoff位置を解決できないsiteは`info.sites[].wyckoff_mapping_error`に理由を保持します。
そのsiteをdisplaciveまたはmagnetic計算に選ぶと、label、座標、理由を示して停止します。

## SiteとModeの選択

利用可能なselectorは`info.sites`から取得します。
固有の`label`は1つの crystallographic siteを選び、`type`は同じtypeのsiteを全て選びます。

```bash
apo info structure.cif | jq '.sites[] | {label, type, wyckoff}'
apo irreps structure.cif --k GM --displacive O --magnetic Fe
```

`--strain`はhomogeneous strainを追加します。
strain-onlyの`irreps`は`GM`を推論します。
atomic modeを含む通常のselectionではsiteとk pointを明示します。

## K pointとParameter

選択可能なlabelとparameter名のauthorityは`kpoints`です。

```bash
apo kpoints structure.cif | \
  jq '.kpoints[] | {label, kvector, parameter_names, miller_love_kvector}'
```

CLI入力はGUIと同じ`kvector`と`parameter_names`の規約を使います。
parameterはその順序で位置指定するか、名前付きで指定します。

```bash
apo irreps structure.cif --k DT 1/3 --displacive O
apo irreps structure.cif --k DT b=1/3 --displacive O
```

`miller_love_kvector`と解決後の`miller_love_parameters`は対応するMiller-Love形式を示します。
CLI入力には常に`kvector`と`parameter_names`の規約を使います。
parametric kが special k pointへ一致する値は拒否されます。
その場合は対応するfixed-point labelを選びます。

1から4個のk/irrep factorをcoupleできます。
繰り返した`--k`と`--irrep`は指定順で対応します。

```bash
apo opds structure.cif \
  --k R --irrep R4- \
  --k M --irrep M3+ \
  --displacive O
```

## OPDとDomain

`opds`が返すlabelが後段で受理されるexact labelです。
表示されたdirectionから labelを再構成しないでください。

```bash
apo opds structure.cif \
  --k R --irrep R4- --displacive O \
  -o opds.json

jq '.opds[] | {label, opd, subgroup, basis, origin, index, cell_index}' \
  opds.json
```

coupled selectionのlabelにはfactorごとのdomain番号が含まれます。
例えば `P1(1)P1(2)`です。
丸括弧を含むlabelはshell quoteしてください。

```bash
apo modes structure.cif \
  --k R --irrep R4- \
  --k M --irrep M3+ \
  --displacive O --opd 'P1(1)P1(2)'
```

`index`はgroup index、`cell_index`はsupercell factorです。
分類済みだが空の `ferroic_properties`と、未分類の結果は`ferroic_classified`で区別できます。

## Mode Details

`modes`は既定でcompact JSONを返します。

```bash
apo modes structure.cif \
  --k R --irrep R4- --displacive O --opd P1 \
  -o modes.json
```

atomic definitionはpayload-localな`definition_id`と構造化されたmode identityを持ちます。
identityにはkind、k point、k vector、irrep、gid、direction、site、Wyckoff position、site irrepが含まれます。
`role`はprimaryと誘起secondaryを区別します。

人間向けmode tableが必要な場合:

```bash
apo modes structure.cif \
  --k R --irrep R4- --displacive O --opd P1 \
  --format text
```

## 保存JSON Case

caseは順序付きの科学selectionを保存します。
parent sourceは`structure`かrepository localな`cif` content IDのどちらか1つです。
`cif`はdevelopment checkout内の既存`Assets/cif` content IDを解決する形式で、standalone利用では通常`structure`を使います。

```json
{
  "structure": "structure.cif",
  "sites": {
    "O": {"displacive": true, "magnetic": false}
  },
  "strain": false,
  "k": [
    {"label": "R", "ir": "R4-"}
  ],
  "opd": "P1"
}
```

指定段まで実行できます。

```bash
apo run --case case.json --upto kpoints
apo run --case case.json --upto irreps
apo run --case case.json --upto opds
apo run --case case.json --upto modes
```

保存caseや自動処理では、stageごとに別processでtableを読み直さない`run`を使います。
個別commandはcaseへ入れるidentifierの探索・確認に使います。

標準入力には`--case -`を使います。
相対structure pathはcase fileのdirectoryから解決します。
k itemは`label`、optionalなexact `params`、OPD段以降で必要な`ir`を持ちます。
保存caseの形式はJSONのみです。

## JSONとOutput Control

compact resultはunversioned schema名でshapeを示します。
consumerはoptional keyから stageを推測せず、このfieldを使います。

```text
APOSTRUCT.cli.info
APOSTRUCT.cli.kpoints
APOSTRUCT.cli.irreps
APOSTRUCT.cli.opds
APOSTRUCT.cli.modes
```

共通option:

```bash
--indent N       JSON indentation
-o, --output     stdoutの代わりにpathへ保存
--full-state     対応commandで完全な計算stateを出力
```

errorはstderrへ出力され、status 2で終了します。
未知のk、irrep、OPDから別結果へsilent fallbackしません。

## Local Interface

```bash
apo serve --host 127.0.0.1 --port 8300 --open-browser
```

GUIはCIFのworkflowを同じ計算経路で提供します。

すでに開いているGUIへcaseを送れます。

```bash
apo show --case case.json --server http://127.0.0.1:8300
```

case内で最も深いselectionまで表示します。
Kだけならirrep選択、irrepまでならOPD選択、`opd`まであればmode detailsとstructure viewerを開きます。
mode sliderの初期値は payload-localなdefinition順で指定します。

```bash
apo show --case case.json \
  --amplitude magnetic-1=1 --amplitude magnetic-2=19/20
```

相対structure pathは送信側commandが解決し、CIF本文をserverへ渡します。
caseと計算stateは稼働中serverのmemoryにだけ保持されます。
