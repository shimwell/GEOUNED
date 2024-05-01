import os
from pathlib import Path

import pytest

from geouned import CadToCsg

path_to_cad = Path("testing/inputSTEP")
step_files = list(path_to_cad.rglob("*.stp")) + list(path_to_cad.rglob("*.step"))
# removing two geometries that are particularly slow to convert from CI testing
# these two geometries remain in the test suite for locally testing
if os.getenv("GITHUB_ACTIONS"):
    step_files.remove(Path("testing/inputSTEP/large/SCDR.stp"))
    step_files.remove(Path("testing/inputSTEP/large/Triangle.stp"))


@pytest.mark.parametrize("input_step_file", step_files)
def test_conversion(input_step_file):
    """Test that step files can be converted to openmc and mcnp files"""

    # sets up an output folder for the results
    output_dir = Path("tests_outputs") / input_step_file.with_suffix("")
    output_dir.mkdir(parents=True, exist_ok=True)
    output_filename_stem = output_dir / input_step_file.stem

    # creates the config file contents
    template = {
        "title": "Input Test",
        "step_file": f"{input_step_file.resolve()}",
        "geometry_name": f"{output_filename_stem.resolve()}",
        "out_format": ("mcnp", "openMC_XML", "openMC_PY", "serpent", "phits"),
        "comp_solids": False,
        "vol_card": False,
        "vol_sdef": True,
        "void_gen": True,
        "dummy_mat": True,
        "min_void_size": 100,
        "cell_summary_file": False,
        "cell_comment_file": False,
        "debug": False,
        "simplify": "no",
        "forceCylinder": False,
        "splitTolerance": 0,
        "newSplitPlane": True,
        "nPlaneReverse": 0,
    }

    # deletes the output MC files if they already exists
    suffixes = (".mcnp", ".xml", ".inp", ".py", ".serp")
    for suffix in suffixes:
        output_filename_stem.with_suffix(suffix).unlink(missing_ok=True)

    GEO = CadToCsg("Input Test")

    # set parameters values stored in template dictionary
    for key, value in template.items():
        GEO.set(key, value)

    GEO.start()

    for suffix in suffixes:
        assert output_filename_stem.with_suffix(suffix).exists()
