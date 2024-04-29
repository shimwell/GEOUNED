

from pathlib import Path

import pytest

from geouned import CadToCsg
from geouned import reverse

path_to_cad = Path("testing/inputSTEP")
step_files = list(path_to_cad.rglob("*.stp")) + list(path_to_cad.rglob("*.step"))


# def test_conversion_and_reverse():
"""Test that step files can be converted to openmc and mcnp files"""

# sets up an output folder for the results
input_step_file = Path('tests/cuboid.stp')
output_dir = Path("tests_outputs") / input_step_file.with_suffix("")
output_dir.mkdir(parents=True, exist_ok=True)
output_filename_stem = output_dir / input_step_file.stem

# deletes the output openmc and mcnp output files if it already exists
output_filename_stem.with_suffix(".mcnp").unlink(missing_ok=True)

GEO = CadToCsg(
    title = 'Input Test',
    stepFile = f"{input_step_file.resolve()}" ,
    geometryName = f"{output_filename_stem.resolve()}" ,
    outFormat = ('mcnp')
)

GEO.Start()

assert output_filename_stem.with_suffix(".mcnp").exists()


template = (
    "[Setting]\n"
    f"inputFile = {output_filename_stem.with_suffix('.mcnp')}\n"
    "CADFile = modelCAD\n"
    "inFormat = mcnp\n"
    "outBox = -500 500 -500 500 0 1000\n"
    "[Levels]\n"
    "UStart = 0\n"
    "levelMax = all\n"
    "[Cells]\n"
    "rangeType = all\n"
    "[Materials]\n"
    "rangeType = exclude\n"
    "range = 0\n"
    "[Options]\n"
    "splitTolerance = 0.01\n"
)

with open(output_dir / "config.ini", mode="w") as file:
    file.write(template)

reverse(output_dir / "config.ini")