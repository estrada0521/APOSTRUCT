"""International space-group setting catalog."""

from __future__ import annotations

from typing import Any

from distortropy.Backend.source.tables import source_tables


def _public_setting(record: dict[str, Any], *, default_id: int) -> dict[str, Any]:
    return {
        "id": int(record["id"]),
        "symbol": str(record["label_short"]),
        "full_symbol": str(record["label_full"]),
        "unique_axis": str(record["axis"]) or None,
        "basis_choice": str(record["abc"]) or None,
        "cell_choice": int(record["cell"]) or None,
        "origin_choice": int(record["origin"]) or None,
        "default": int(record["id"]) == int(default_id),
    }


def international_settings(space_group: int) -> dict[str, Any]:
    """Return COPL-compatible International settings for one space group."""

    sg = int(space_group)
    tables = source_tables()
    default = tables.default_inter_setting_record(sg)
    rows = tables.international_setting_records_for_space_group(sg)
    return {
        "schema": "distortropy.cli.settings",
        "space_group": {
            "number": sg,
            "symbol": str(default["label_short"]),
            "default_setting": int(default["id"]),
        },
        "settings": [
            _public_setting(record, default_id=int(default["id"]))
            for record in rows
        ],
    }


def public_setting(space_group: int, setting_id: int | None) -> dict[str, Any]:
    """Return one validated setting in the compact directions schema."""

    tables = source_tables()
    default = tables.default_inter_setting_record(int(space_group))
    record = tables.international_setting_record(int(space_group), setting_id)
    return _public_setting(record, default_id=int(default["id"]))


__all__ = ["international_settings", "public_setting"]
