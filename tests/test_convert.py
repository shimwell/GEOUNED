import os
from pathlib import Path

import pytest

import geouned

path_to_cad = Path("testing/inputSTEP")
step_files = list(path_to_cad.rglob("*.stp")) + list(path_to_cad.rglob("*.step"))
# removing two geometries that are particularly slow to convert from CI testing
# these two geometries remain in the test suite for locally testing
if os.getenv("GITHUB_ACTIONS"):
    step_files.remove(Path("testing/inputSTEP/large/SCDR.stp"))
    step_files.remove(Path("testing/inputSTEP/large/Triangle.stp"))
suffixes = (".mcnp", ".xml", ".inp", ".py", ".serp")


@pytest.mark.parametrize("input_step_file", step_files)
def test_conversion(input_step_file):
    """Test that step files can be converted to openmc and mcnp files"""

    # sets up an output folder for the results
    output_dir = Path("tests_outputs") / input_step_file.with_suffix("")
    output_dir.mkdir(parents=True, exist_ok=True)
    output_filename_stem = output_dir / input_step_file.stem

    # deletes the output MC files if they already exists
    for suffix in suffixes:
        output_filename_stem.with_suffix(suffix).unlink(missing_ok=True)

    my_options = geouned.Options(
        force_cylinder=False,
        new_split_plane=True,
        del_last_number=False,
        enlarge_box=2,
        n_plane_reverse=0,
        split_tolerance=0,
        scale_up=True,
        quadric_py=False,
        facets=False,
        prnt_3p_plane=False,
        forceNoOverlap=False,
    )

    my_tolerances = geouned.Tolerances(
        relative_tol=False,
        relative_precision=0.000001,
        value=0.000001,
        distance=0.0001,
        angle=0.0001,
        pln_distance=0.0001,
        pln_angle=0.0001,
        cyl_distance=0.0001,
        cyl_angle=0.0001,
        sph_distance=0.0001,
        kne_distance=0.0001,
        kne_angle=0.0001,
        tor_distance=0.0001,
        tor_angle=0.0001,
        min_area=0.01,
    )

    my_numeric_format = geouned.NumericFormat(
        p_abc="14.7e",
        p_d="14.7e",
        p_xyz="14.7e",
        s_r="14.7e",
        s_xyz="14.7e",
        c_r="12f",
        c_xyz="12f",
        k_xyz="13.6e",
        k_tan2="12f",
        t_r="14.7e",
        t_xyz="14.7e",
        gq_1_to_6="18.15f",
        gq_7_to_9="18.15f",
        gq_10="18.15f",
    )

    my_settings = geouned.Settings(
        mat_file="",
        void_gen=True,
        debug=False,
        comp_solids=True,
        simplify="no",
        cell_range=[],
        export_solids="",
        min_void_size=200.0,  # units mm
        max_surf=50,
        max_bracket=30,
        void_mat=[],
        void_exclude=[],
        start_cell=1,
        start_surface=1,
        sort_enclosure=False,
    )

    geo = geouned.CadToCsg(
        step_filename=f"{input_step_file.resolve()}",
        options=my_options,
        settings=my_settings,
        tolerances=my_tolerances,
        numeric_format=my_numeric_format,
    )

    geo.start()

    geo.export_csg(
        title="Converted with GEOUNED",
        geometryName=f"{output_filename_stem.resolve()}",
        outFormat=(
            "openmc_xml",
            "openmc_py",
            "serpent",
            "phits",
            "mcnp",
        ),
        volSDEF=True,  # changed from the default
        volCARD=False,  # changed from the default
        UCARD=None,
        dummyMat=True,  # changed from the default
        cellCommentFile=False,
        cellSummaryFile=False,  # changed from the default
    )

    for suffix in suffixes:
        assert output_filename_stem.with_suffix(suffix).exists()


@pytest.mark.parametrize(
    "input_json_file",
    ["tests/config_complete_defaults.json", "tests/config_minimal.json"],
)
def test_cad_to_csg_from_json_with_defaults(input_json_file):

    # deletes the output MC files if they already exists
    for suffix in suffixes:
        Path("csg").with_suffix(suffix).unlink(missing_ok=True)

    my_cad_to_csg = geouned.CadToCsg.from_json(input_json_file)
    assert isinstance(my_cad_to_csg, geouned.CadToCsg)

    assert my_cad_to_csg.step_filename == "testing/inputSTEP/BC.stp"
    assert my_cad_to_csg.options.force_cylinder == False
    assert my_cad_to_csg.tolerances.relative_tol == False
    assert my_cad_to_csg.numeric_format.p_abc == "14.7e"
    assert my_cad_to_csg.settings.mat_file == ""

    for suffix in suffixes:
        assert Path("csg").with_suffix(suffix).exists()

    # deletes the output MC files if they already exists
    for suffix in suffixes:
        Path("csg").with_suffix(suffix).unlink(missing_ok=True)

    my_cad_to_csg.start()
    my_cad_to_csg.export_csg()


def test_cad_to_csg_from_json_with_non_defaults():

    # deletes the output MC files if they already exists
    for suffix in suffixes:
        Path("csg").with_suffix(suffix).unlink(missing_ok=True)

    my_cad_to_csg = geouned.CadToCsg.from_json("tests/config_non_defaults.json")
    assert isinstance(my_cad_to_csg, geouned.CadToCsg)

    assert my_cad_to_csg.step_filename == "testing/inputSTEP/BC.stp"
    assert my_cad_to_csg.options.force_cylinder == True
    assert my_cad_to_csg.tolerances.relative_precision == 2e-6
    assert my_cad_to_csg.numeric_format.p_abc == "15.7e"
    assert my_cad_to_csg.settings.mat_file == "non default"

    for suffix in suffixes:
        assert Path("csg").with_suffix(suffix).exists()

    # deletes the output MC files if they already exists
    for suffix in suffixes:
        Path("csg").with_suffix(suffix).unlink(missing_ok=True)

    my_cad_to_csg.start()
    my_cad_to_csg.export_csg()


def test_writing_to_new_folders():
    """Checks that a folder is created prior to writing output files"""

    geo = geouned.CadToCsg(step_filename="testing/inputSTEP/BC.stp")
    geo.start()

    for outformat in ["mcnp", "phits", "serpent", "openmc_xml", "openmc_py"]:
        geo.export_csg(
            geometryName=f"tests_outputs/new_folder_for_testing_{outformat}/csg",
            cellCommentFile=False,
            cellSummaryFile=False,
            outFormat=[outformat],
        )
        geo.export_csg(
            geometryName=f"tests_outputs/new_folder_for_testing_{outformat}_cell_comment/csg",
            cellCommentFile=True,
            cellSummaryFile=False,
            outFormat=[outformat],
        )
        geo.export_csg(
            geometryName=f"tests_outputs/new_folder_for_testing_{outformat}_cell_summary/csg",
            cellCommentFile=False,
            cellSummaryFile=True,
            outFormat=[outformat],
        )


def test_with_relative_tol_true():

    # test to protect against incorrect attribute usage in FreeCAD
    # more details https://github.com/GEOUNED-org/GEOUNED/issues/154

    geo = geouned.CadToCsg(
        step_filename=f"{step_files[1].resolve()}",
        tolerances=geouned.Tolerances(relative_tol=False),
    )
    geo.start()
    geo = geouned.CadToCsg(
        step_filename=f"{step_files[1].resolve()}",
        tolerances=geouned.Tolerances(relative_tol=True),
    )
    geo.start()
