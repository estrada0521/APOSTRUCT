# Magnetic group generators

`magnetic_group_generators.json.gz` contains exact fractional-coordinate generators for all 1,651 BNS magnetic space groups, 4,551 coordinate presentations, and an index from the 230 ordinary space-group numbers to their type-I BNS groups.
It is compiled with spglib 2.7.0 from the tracked Hall and BNS tables in `upstream/spglib-2.7.0`, taken byte-for-byte from spglib commit `12355c77fb7c505a55f52cae36341d73b781a065` (tag `v2.7.0`).
The BSD 3-Clause license is retained in `SPGLIB-LICENSE.txt`.

The presentation set is the union of the 527 distinct ordinary Hall operation sets and the 15 current, nonlegacy Source operation sets not represented by a Hall number.
Source's 742 current setting rows collapse to 542 exact operation sets before the difference is taken; aliases are not scanned separately.
Each Source-only presentation is transported from an exact Hall-overlap setting, including the antiunitary signs of magnetic groups.

The upstream magnetic Hall table needs nine time-reversal prime corrections to match the ISO-MAG operation table exactly.
The corrected entries and all input hashes are embedded in the compressed document.

The tracked builder uses exact rational Hall parsing, converts spglib's sixth-rational setting translations back to `Fraction`, reconstructs every closure, selects generators deterministically, and reproduces the committed gzip byte-for-byte:

```bash
python3 Verification/mathematics/build_magnetic_group_data.py
```

For the independent upstream audit, obtain ISO-MAG's computer-readable `magnetic_data.txt` with SHA-256 `2b11217ae10687b0836d8151846db2fad57b0d77211da90bb138ccd123a7b1fc` and run:

```bash
python3 Verification/mathematics/build_magnetic_group_data.py \
  --official-iso-mag /path/to/magnetic_data.txt
```

That audit parses the official table independently and requires exact equality of rotation matrices, translations modulo one, and time-reversal signs for all 1,651 groups and 38,307 operations.
It uses no setting conversion or numerical tolerance.
`--write` is accepted only together with this complete audit.

Twelve BNS groups have a known time-reversal placement disagreement between spglib's setting API and the independently audited ISO-MAG closure.
For those groups the artifact retains only the corrected, audited default setting and marks setting coverage incomplete.
The verifier may satisfy a proof in that setting, but it reports `indeterminate` rather than `refuted` when no retained setting satisfies it.
No unverified alternate setting can create a refutation.

The runtime verifier reconstructs every retained setting closure and checks its operation count and SHA-256 identity.
A displayed atomic field is refuted only when no declared presentation first fits the displayed undistorted structure and then satisfies the rounded mode coordinates and vectors.
The table is proof authority only; it is not imported by APOSTRUCT computation.
