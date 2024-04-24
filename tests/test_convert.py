from pathlib import Path

import pytest

from geouned import Geouned

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
        "[files]\n"
        "title = 'Input Test'\n"
        f"step_file = {input_step_file.resolve()}\n"
        f"geometry_name = {output_filename_stem.resolve()}\n"
        "out_format = ('mcnp', 'openMC_XML')\n"
        "[parameters]\n"
        "comp_solids = False\n"
        "vol_card = False\n"
        "vol_sdef = True\n"
        "void_gen = True\n"
        "dummy_mat = True\n"
        "min_void_size = 100\n"
        "cell_summary_file = False\n"
        "cell_comment_file = False\n"
        "debug = False\n"
        "simplify = no\n"
        "[options]\n"
        "force_cylinder = False\n"
        "split_tolerance = 0\n"
        "new_split_plane = True\n"
        "n_plane_reverse = 0\n"
    )

    with open(output_dir / "config.ini", mode="w") as file:
        file.write(template)

    # deletes the output openmc and mcnp output files if it already exists
    output_filename_stem.with_suffix(".mcnp").unlink(missing_ok=True)
    output_filename_stem.with_suffix(".xml").unlink(missing_ok=True)

    inifile = f"{output_dir/'config.ini'}"
    GEO = Geouned(inifile)
    GEO.set_options()
    GEO.out_format = ("mcnp", "openMC_XML")
    GEO.start()

    assert output_filename_stem.with_suffix(".mcnp").exists()
    assert output_filename_stem.with_suffix(".xml").exists()
