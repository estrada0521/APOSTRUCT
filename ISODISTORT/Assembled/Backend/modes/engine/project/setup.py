"""Setup and count-adjacent helpers for binary ``project_`` entry."""

from __future__ import annotations

from typing import Iterable

from ISODISTORT.Assembled.Backend.modes.engine.records import LittleIrrepRecord, WyckoffRow


class ProjectSetupMixin:
    def project_parent_operation_records_for_little_row(
        self,
        sg: int,
        little: LittleIrrepRecord,
        row: WyckoffRow,
        site_params: Iterable[object] | None = None,
    ) -> tuple[tuple[int, int, int, int, int], ...]:
        """Return full parent operation records in `project_` setting order."""

        _ = little, site_params
        return self.wyc_pg_elements_records(sg, row)

    def map_parent_ops_to_site_project_records(
        self,
        site_pg: int,
        parent_point_ops: Iterable[int],
    ) -> tuple[tuple[int, int, int, int, int], ...]:
        """Map parent-space operations to site-irrep operations.

        This mirrors the `project_` setting loop.  It tries each candidate
        setting from `iwyckoff_pg_elements` until all setting point operations
        are found in the parent operation list.  The returned operation records
        keep the parent operation order but replace each matched parent point
        op by the operation in the canonical first site setting at the same
        position.

        SG205/R/Wyckoff-a example:

        parent `[1,9,5,25,33,29]` matches setting 2 and maps to
        `[49,51,53,61,63,65]`, which is exactly the GDB `local_408` dump.
        """

        parent = [int(op) for op in parent_point_ops]
        settings = self.site_pg_element_settings(site_pg)
        if not settings:
            raise ValueError(f"site point group {site_pg} has no element settings")
        canonical = settings[0]
        for setting in settings:
            if len(setting) != len(canonical):
                continue
            mapped: list[tuple[int, int, int, int, int] | None] = [None] * len(parent)
            ok = True
            for setting_op, canonical_op in zip(setting, canonical):
                try:
                    parent_index = parent.index(setting_op)
                except ValueError:
                    ok = False
                    break
                mapped[parent_index] = (0, 0, 0, 1, canonical_op)
            if ok:
                return tuple(record for record in mapped if record is not None)
        raise ValueError(f"no site setting for PG{site_pg} matches parent ops {parent}")

    def vector_basis_id_for_site_records(
        self,
        site_pg: int,
        site_records: Iterable[tuple[int, int, int, int, int]],
    ) -> int:
        """Return the vector-basis row id used by MAIN__ before `project_vector_`.

        The row id is aligned with the site point-group setting order in
        `iwyckoff_pg_elements`.  MAIN__ later passes this id as the first
        integer of `project_vector_` argument 2; `project_vector_` uses it to
        index `iwyckoff_pg_vector_basis` in 18-int blocks.
        """

        ops = tuple(int(record[4]) for record in site_records)
        pointer = int(self.iso.wyckoff["iwyckoff_pg_elements_pointer"][int(site_pg) - 1])
        for index, setting in enumerate(self.site_pg_element_settings(site_pg)):
            if tuple(setting) == ops:
                return pointer + index
        for index, setting in enumerate(self.site_pg_element_settings(site_pg)):
            if set(setting) == set(ops):
                return pointer + index
        raise ValueError(f"no vector basis setting for PG{site_pg} ops={ops}")

    def project_site_matrix_ops(
        self,
        site_pg: int,
        pg_irrep: int,
        site_project_records: Iterable[tuple[int, int, int, int, int]],
    ) -> tuple[int, ...]:
        """Map local_408 operation records to local matrix-table ops.

        `project_` passes the `local_408` records to `get_irreps_`.  For some
        site point groups, those records are not the same operation labels as
        the finite coset reps stored in `data_irreps` for the corresponding
        site irrep.  The mapping is position-wise between the site setting used
        in `local_408` and the setting whose row equals the matrix op order.
        """

        local_ops = tuple(int(record[4]) for record in site_project_records)
        matrix_ops = self.site_irrep_operation_order(site_pg, pg_irrep)
        settings = self.site_pg_element_settings(site_pg)
        local_setting = next((setting for setting in settings if setting == local_ops), None)
        matrix_setting = next((setting for setting in settings if setting == matrix_ops), None)
        if local_setting is None:
            raise ValueError(f"local_408 ops {local_ops} are not a setting for PG{site_pg}")
        if matrix_setting is None:
            raise ValueError(f"matrix ops {matrix_ops} are not a setting for PG{site_pg}")
        if len(local_setting) != len(matrix_setting):
            raise ValueError(f"setting length mismatch for PG{site_pg}: {local_setting} vs {matrix_setting}")
        return tuple(matrix_setting[local_setting.index(op)] for op in local_ops)

    def estimated_project_counts(self, sg: int, wyckoff: str, label: str) -> list[tuple[int, int, int]]:
        """Return ``(pg_irrep, frequency, estimated_project_count)``.

        This mirrors the count branch observed in ``project_``:

        ``frequency`` is doubled when ``image_type != 1``, then divided by the
        real/complex type divisor.  The offset mapping is from the read_images
        load addresses:

        * ``common_image + 0x20e71 + image_id`` -> ``image_dimension``
        * ``common_image + 0x21587 + image_id`` -> ``image_type``
        """

        little = self.little_record(sg, label)
        divisor = self.TYPE_DIVISOR[little.irrep_type]
        rows = [item for item in self.wyckoff_subductions(sg, label) if item.wyckoff.label == wyckoff]
        if len(rows) != 1:
            raise KeyError(f"expected one Wyckoff row for SG{sg} {wyckoff}, got {rows}")
        out = []
        for pg_irrep, frequency in rows[0].pairs:
            site_old_id = self.site_pg_irrep_old_id(rows[0].wyckoff.site_pg, pg_irrep)
            site_image = self.image_record(site_old_id)
            factor = 1 if site_image.image_type == 1 else 2
            numerator = frequency * factor
            if numerator % divisor:
                raise ValueError(
                    f"nonintegral project count for SG{sg} {wyckoff} {label}: "
                    f"frequency={frequency} factor={factor} divisor={divisor}"
                )
            out.append((pg_irrep, frequency, numerator // divisor))
        return out


