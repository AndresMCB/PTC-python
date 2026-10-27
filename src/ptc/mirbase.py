"""Conversion of miRNA names to miRBase version 21, as PTC.miRv21() in R.

PTC.miRv21() uses miRBaseConverter::miRNAVersionConvert(). The conversion
table in data/miRBase_v21.csv.gz was exported from miRBaseConverter
(scripts/export_from_R.R) and this module applies the same rules, so the
converted names are the same as in R.
"""

import csv
import gzip
import re
from functools import lru_cache
from importlib import resources


@lru_cache(maxsize=1)
def _table():
    """miRNA name -> list of miRBase v21 names, in miRBaseConverter's order."""
    path = resources.files("ptc") / "data" / "miRBase_v21.csv.gz"
    with gzip.open(path, "rt", newline="") as f:
        reader = csv.reader(f)
        next(reader)
        return {name: targets.split("&") for name, targets in reader}


def mirv21(names):
    """Returns the miRBase v21 version of a list of miRNA names.

    As in R, the part of a name after a dot is removed, names that cannot be
    converted are kept, and when a name maps to several miRNAs the first one
    is used. Like miRNAVersionConvert(exact = TRUE), a name that is itself the
    new name of another miRNA in the list is only kept as it is.
    """
    table = _table()
    names = [re.sub(r"\..*", "", str(n)) for n in names]
    keys = [n.replace(" ", "") for n in names]
    # all v21 names that the miRNAs in this list map to
    in_list = {t for k in set(keys) if k in table for t in table[k]}
    converted = []
    for name, key in zip(names, keys):
        targets = [t for t in table.get(key, []) if t == key or key not in in_list]
        converted.append(targets[0] if targets else name)
    return converted
