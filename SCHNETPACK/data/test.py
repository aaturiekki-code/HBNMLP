#!/usr/bin/env python3

from pathlib import Path
from ase.db import connect

path = Path(
    "/scratch/work/riekkia3/Special_Assignment/SCHNETPACK/data/qm9.db"
)

if not path.is_file():
    raise FileNotFoundError(path)

with connect(str(path)) as db:
    print("Number of structures:", db.count())
    print("Metadata keys:", list(db.metadata))