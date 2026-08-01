# Distortropy CLIガイド

[English](CLI.md) | [README](README_ja.md)

CLIはlocal GUIと同じ計算serviceを公開します。各commandはstatelessです。前段の
JSONからlabelやselectorを選び、それを次のcommandへ渡します。

```text
parent -> k points -> irreps -> OPDs -> modes
                                  `----> invariants
```

全optionは`distortropy COMMAND --help`で確認できます。この文書では各commandに
共通する科学的・machine-readable契約を説明します。

## Command一覧

| Command | 結果 |
|---|---|
| `settings` | 1つのspace groupに対するCOPL互換International setting一覧 |
| `directions` | exactなG:H embeddingと両立するordinary/time-odd OP direction |
| `info` | 正規化されたparent cellと選択可能なcrystallographic site |
| `kpoints` | k-point label、式、star、parameter名 |
| `irreps` | 選択したmode kindで発現するordinary/magnetic irrep |
| `opds` | OPD、subgroup、basis、origin、domain、ferroic metadata |
| `modes` | 1つのOPDに対するatomic/strain mode definition |
| `invariants` | 選択したOPD domain上のinvariant basis |
| `combine-modes` | compact mode definitionを重ねた原子ごとのvector |
| `run` | 保存したJSON caseを指定段まで実行 |
| `serve` | local GUIを起動 |

## Subgroupと両立するDirection

`directions`は、具体構造に依存しない逆向きの問いに答えます。すなわち、exactな
parent-to-subgroup embeddingと両立するordinaryまたはtime-odd
order-parameter directionを返します。

parentまたはsubgroupがdefault settingでない場合、先に受理されるInternational
setting IDを確認します。

```bash
distortropy settings --sg 14
```

```bash
distortropy directions \
  --parent-sg 14 --parent-setting 64 \
  --subgroup-sg 2 --subgroup-setting 2 \
  --basis=0,0,3 --basis=0,1,0 --basis=-1,0,0
```

setting optionを省略するとdefaultを使います。

```bash
distortropy directions \
  --parent-sg 221 --subgroup-sg 140 \
  --basis=-1,1,0 --basis=-1,-1,0 --basis=0,0,2 \
  --origin=0,0,0
```

magnetic childではordinary subgroup番号の代わりにBNS番号を指定します。parentは
`--parent-sg`のparamagnetic gray extensionです。

```bash
distortropy directions \
  --parent-sg 62 --subgroup-msg 62.448 \
  --basis=1,0,0 --basis=0,1,0 --basis=0,0,1
```

basis vectorとoriginはparent/subgroup双方の選択したInternational settingで指定します。
setting IDはCOPLと同じものです。同じsubgroup型でもparent内のorientation、cell、
origin、domainが異なり得るため、space-group番号だけでは入力として不十分です。
CIF、Wyckoff site、distortion modeの選択は使いません。

各rowはk point、irrep、OPD label、structured direction matrix、stabilizer subgroup、group index、
cell indexを返します。direction matrixのrowはfull irrep coordinate、columnは
`parameters`に対応し、parametric rowはpublic/Miller-Love両方のparameter値を保持します。
Sourceに対応するOPD-family labelがないdirectionでは`opd`はnullです。
`primary`はそのdirection単独が指定embeddingと完全に同じsubgroupを持つことを示し、
該当rowが複数の場合もあります。その非nullな`domain`は`modes --from-directions`が
OPD catalogを再列挙せず再利用するSource domain番号です。`secondary`はstabilizerがそのsupergroupだが指定subgroupで
許容されることを示します。ordinary rowはISOTROPY `DISPLAY DIRECTION`に対応します。
magnetic rowは同じlocal irrep matrixとexact BNS operationを用い、bundled static magnetic
OPD catalog全体に対して自己検証します。ISO 9.6.1には対応するmagnetic
`DISPLAY DIRECTION` oracleがありません。superspaceの逆引きは対象外です。

## Parent入力

### CIF

CIF pathを位置引数へ渡します。

```bash
distortropy info structure.cif
distortropy kpoints structure.cif
```

parent setting、site、occupancy、座標はCIFをauthorityとします。

### Space groupとWyckoff site

具体構造を必要としない場合は`--sg`を使います。

```bash
distortropy kpoints --sg 221
distortropy info --sg 221 --wyckoff 1a 1b 3c
```

symbolic routeはdefault International Tables settingを使います。Wyckoff letterと
multiplicity+letterの両方を受け付けます。入力の自由座標は`invariants`まで未指定の
まま進められます。

```bash
distortropy opds \
  --sg 137 --wyckoff 2a 4d \
  --k M --irrep M1 --displacive a d
```

`info`はそのsiteを未実体化として示し、識別用のgeneric placeholder座標を表示する
場合があります。これは物理構造ではなく、k point、irrep、OPD、invariant計算の
geometryとしては使用されません。`modes`はatomic geometryを実体化するため、全ての
自由座標が必要です。

```bash
distortropy modes \
  --sg 137 --wyckoff 2a 4d:z=1/7 \
  --k M --irrep M1 --displacive a d --opd P1
```

ここで与えた値は具体的なgeometryです。generic placeholderとしては扱いません。

## SiteとModeの選択

利用可能なselectorは`info.sites`から取得します。固有の`label`は1つの
crystallographic siteを選び、`type`は同じtypeのsiteを全て選びます。

```bash
distortropy info structure.cif | jq '.sites[] | {label, type, wyckoff}'
distortropy irreps structure.cif --k GM --displacive O --magnetic Fe
```

同じsymbolic Wyckoff orbitを複数指定した場合、`i1`、`i2`のようなlabelで個別に
選べます。type `i`は両方を選びます。`--strain`はhomogeneous strainを追加します。
site-free strain-onlyの`irreps --sg SG --strain`は`GM`を推論します。atomic modeを
含む通常のselectionではsiteとk pointを明示します。

## K pointとParameter

選択可能なlabelとparameter名のauthorityは`kpoints`です。

```bash
distortropy kpoints --sg 225 | \
  jq '.kpoints[] | {label, kvector, parameter_names, miller_love_kvector}'
```

CLI入力はGUIと同じ`kvector`と`parameter_names`の規約を使います。parameterはその
順序で位置指定するか、名前付きで指定します。

```bash
distortropy irreps structure.cif --k DT 1/3 --displacive O
distortropy irreps structure.cif --k DT b=1/3 --displacive O
```

`miller_love_kvector`と解決後の`miller_love_parameters`はSource内部座標のprovenance
です。そのsymbolをCLI入力へ直接代入しないでください。parametric kがspecial
k pointへ一致する値は拒否されます。その場合は対応するfixed-point labelを選びます。

1から4個のk/irrep factorをcoupleできます。繰り返した`--k`と`--irrep`は指定順で
対応します。

```bash
distortropy opds structure.cif \
  --k R --irrep R4- \
  --k M --irrep M3+ \
  --displacive O
```

## OPDとDomain

`opds`が返すlabelが後段で受理されるexact labelです。表示されたdirectionから
labelを再構成しないでください。

```bash
distortropy opds structure.cif \
  --k R --irrep R4- --displacive O \
  -o opds.json

jq '.opds[] | {label, opd, subgroup, basis, origin, index, cell_index}' \
  opds.json
```

coupled selectionのlabelにはfactorごとのdomain番号が含まれます。例えば
`P1(1)P1(2)`です。丸括弧を含むlabelはshell quoteしてください。

```bash
distortropy modes structure.cif \
  --k R --irrep R4- \
  --k M --irrep M3+ \
  --displacive O --opd 'P1(1)P1(2)'
```

`index`はgroup index、`cell_index`はsupercell factorです。分類済みだが空の
`ferroic_properties`と、未分類の結果は`ferroic_classified`で区別できます。

## Mode Details

`modes`は既定でcompact JSONを返します。

```bash
distortropy modes structure.cif \
  --k R --irrep R4- --displacive O --opd P1 \
  -o modes.json
```

報告されたsubgroup embeddingがcatalog OPD代表とは異なるsymmetry-equivalent domainを
選ぶ場合、OPD labelへ戻さず、exactな`directions`結果を再利用できます。

```bash
distortropy directions \
  --parent-sg 225 --subgroup-sg 87 \
  --basis=1/2,-1/2,0 --basis=1/2,1/2,0 --basis=0,0,1 \
  -o directions.json

distortropy modes structure.cif \
  --from-directions directions.json --direction-row 3 \
  --displacive Br -o modes.json
```

`direction-row`は1-basedの`directions[].row`値で、primary rowを選ぶ必要があります。
保存されたexact basis、origin、setting、direction subspace、Source domainを既存mode核へ
そのまま渡します。この経路ではdirectな`--k`、`--irrep`、`--opd`は不要で、併用できません。

atomic definitionはpayload-localな`definition_id`と構造化されたmode identityを
持ちます。identityにはkind、k point、k vector、irrep、gid、direction、site、
Wyckoff position、site irrepが含まれます。`role`はprimaryと誘起secondaryを区別し、
帰属が一意なprimaryはinvariant factorのslotとglobal parameter名も持ちます。

人間向けmode tableが必要な場合:

```bash
distortropy modes structure.cif \
  --k R --irrep R4- --displacive O --opd P1 \
  --format text
```

### Definitionの線形結合

compact modes payloadからdefinition IDを確認します。

```bash
jq '.mode_details.displacive_definitions[]?,
    .mode_details.magnetic_definitions[]? |
    {definition_id, normfactor, role, mode}' modes.json
```

`--weight`は公開されたrowへ係数を直接掛けます。

```bash
distortropy combine-modes modes.json \
  --weight magnetic-1=1 --weight magnetic-2=3/4
```

`--amplitude`はdefinitionごとのpublished mode normfactorを先に適用します。

```bash
distortropy combine-modes modes.json \
  --amplitude magnetic-1=1 --amplitude magnetic-2=1
```

このamplitudeは各definitionのmode normalization規約であり、無関係なmode間に
共通の物理単位を与えるものではありません。返されるvector componentはchild
crystallographic `dxyz` basisです。net magnetic vectorは、返されたconventional
child-cell atomについての和です。

## Invariants

通常の`invariants`は、選択したprimary factorとexact OPD domainを使います。

```bash
distortropy invariants structure.cif \
  --k GM --irrep GM4- \
  --k GM --irrep GM3+ \
  --displacive O --strain \
  --opd 'P1(1)P1(3)' \
  --minimum-degree 1 --maximum-degree 4
```

degree 1から12を計算できます。要求範囲の全degreeが出力され、不変量がないdegreeも
`count: 0`、`invariants: []`、`polynomials: []`として残ります。表示文字列に加え、
machine用のstructured polynomialを返します。

```json
{
  "terms": [
    {"coefficient": "1", "exponents": [3, 0]},
    {"coefficient": "-3", "exponents": [1, 2]}
  ]
}
```

exponentの位置はtop-level `variables`の順です。coefficientはfloatではなくexactな
SymPy expression stringです。

### 誘起secondary factorを含める

`modes`はprimaryとsecondaryの`invariant_factors`を順序付きで保持します。各factor
にはexact domainとglobal parameter offsetが含まれます。このpayloadを再利用すれば、
parent、OPD、mode計算を再実行せずにsecondaryを追加できます。

```bash
jq '.invariant_factors[] |
    {slot, role, label, opd, domain, parameter_offset}' modes.json

distortropy invariants \
  --from-modes modes.json \
  --secondary 4 \
  --minimum-degree 2 --maximum-degree 4
```

`--secondary`は繰り返し可能で、payload-localなfactor slotを取ります。primaryは常に
含まれます。secondary domainの手動overrideは認めず、modes結果が保持したdomainを
authorityとします。

## 保存JSON Case

caseは順序付きの科学selectionを保存します。parent sourceは`structure`、repository
localな`cif` content ID、またはoptional `wyckoff`を伴う`sg`のどれか1つです。
`cif`はdevelopment checkout内の既存`Assets/cif` content IDを解決する形式で、
standalone利用では通常`structure`を使います。

```json
{
  "sg": 221,
  "wyckoff": ["1a", "1b", "3c"],
  "sites": {
    "c": {"displacive": true, "magnetic": false}
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
distortropy run --case case.json --upto kpoints
distortropy run --case case.json --upto irreps
distortropy run --case case.json --upto opds
distortropy run --case case.json --upto modes
distortropy run --case case.json --upto invariants \
  --minimum-degree 2 --maximum-degree 6
```

標準入力には`--case -`を使います。相対structure pathはcase fileのdirectoryから解決
します。k itemは`label`、optionalなexact `params`、OPD段以降で必要な`ir`を持ちます。
保存caseの形式はJSONのみです。

## JSONとOutput Control

compact resultはunversioned schema名でshapeを示します。consumerはoptional keyから
stageを推測せず、このfieldを使います。

```text
distortropy.cli.settings
distortropy.cli.directions
distortropy.cli.info
distortropy.cli.kpoints
distortropy.cli.irreps
distortropy.cli.opds
distortropy.cli.modes
distortropy.cli.invariants
distortropy.cli.mode_combination
```

共通option:

```bash
--indent N       JSON indentation
-o, --output     stdoutの代わりにpathへ保存
--full-state     対応commandで完全な内部pipeline stateを出力
```

`invariants`と`combine-modes`はcompact-onlyです。そのため
`run --upto invariants`は`--full-state`を拒否します。errorはstderrへ出力され、
status 2で終了します。未知のk、irrep、OPDから別結果へsilent fallbackしません。

## Local Interface

```bash
distortropy serve --host 127.0.0.1 --port 8300 --open-browser
```

GUIはCIFとsymbolic parentの両workflowを同じbackend上で提供します。debug表示は
provenanceを追加するだけで、別の科学計算経路を選ぶものではありません。
