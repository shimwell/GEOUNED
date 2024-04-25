from pathlib import Path

import pytest

from geouned import GEOUNED

path_to_cad = Path("testing/inputSTEP")
step_files = list(path_to_cad.rglob("*.stp")) + list(path_to_cad.rglob("*.step"))


@pytest.mark.parametrize("input_step_file", step_files)
def test_conversion(input_step_file):
    """Test that step files can be converted to openmc and mcnp files"""

    # sets up an output folder for the results
    output_dir = Path("tests_outputs") / input_step_file.with_suffix("")
    output_dir.mkdir(parents=True, exist_ok=True)
    output_filename_stem = output_dir / input_step_file.stem

    # creates the config file contents
    template = (
        "[Files]\n"
        "title = 'Input Test'\n"
        f"stepFile = {input_step_file.resolve()}\n"
        f"geometryName = {output_filename_stem.resolve()}\n"
        "outFormat = ('mcnp', 'openMC_XML')\n"
        "[Parameters]\n"
        "compSolids = False\n"
        "volCARD = False\n"
        "volSDEF = True\n"
        "voidGen = True\n"
        "dummyMat = True\n"
        "minVoidSize = 100\n"
        "cellSummaryFile = False\n"
        "cellCommentFile = False\n"
        "debug = False\n"
        "simplify = no\n"
        "[Options]\n"
        "forceCylinder = False\n"
        "splitTolerance = 0\n"
        "newSplitPlane = True\n"
        "nPlaneReverse = 0\n"
    )

    with open(output_dir / "config.ini", mode="w") as file:
        file.write(template)

    # deletes the output openmc and mcnp output files if it already exists
    output_filename_stem.with_suffix(".mcnp").unlink(missing_ok=True)
    output_filename_stem.with_suffix(".xml").unlink(missing_ok=True)

    inifile = f"{output_dir/'config.ini'}"
    geo = GEOUNED(inifile)
    geo.SetOptions()
    geo.outFormat = ("mcnp", "openMC_XML")
    geo.Start()

    assert output_filename_stem.with_suffix(".mcnp").exists()
    assert output_filename_stem.with_suffix(".xml").exists()


def test_tolerances_are_no_persistant():
    template1 = (
        "[Files]\n"
        "title = 'Input Test'\n"
        f"stepFile = {step_files[2]}\n"
    )

    with open("config1.ini", mode="w") as file:
        file.write(template1)

    geo1 = GEOUNED('config1.ini')
    geo1.SetOptions()  # this should set the default tolerances.relativePrecision

    assert geo1.tolerances.relativePrecision == 1e-6

    template2 = (
        "[Files]\n"
        "title = 'Input Test'\n"
        f"stepFile = {step_files[2]}\n"
        "[Tolerances]\n"
        "relativePrecision = 2e-6"  # here we specify a new number for relativePrecision
    )

    with open("config2.ini", mode="w") as file:
        file.write(template2)
    geo2 = GEOUNED('config2.ini')
    geo2.SetOptions()  # this should set 2e-6 as the tolerances.relativePrecision

    assert geo2.tolerances.relativePrecision == 2e-6

    geo3 = GEOUNED('config1.ini')
    geo3.SetOptions()  # this does not change the relativePrecision back to the default
    assert geo3.tolerances.relativePrecision == 1e-6
    assert geo2.tolerances.relativePrecision == 2e-6
    assert geo1.tolerances.relativePrecision == 1e-6
