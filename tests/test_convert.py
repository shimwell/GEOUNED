from pathlib import Path

import pytest

from geouned import CadToCsg, Options

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
    template = {
        "compSolids" : False , # Parameters
        "volCARD" : False , # Parameters
        "volSDEF" : True , # Parameters
        "voidGen" : True , # Parameters
        "dummyMat" : True , # Parameters
        "minVoidSize" : 100 , # Parameters
        "cellSummaryFile" : False , # Parameters
        "cellCommentFile" : False , # Parameters
        "debug" : False , # Parameters
        "simplify" : 'no' , # Parameters
    }
  
    my_options = Options(
      forceCylinder=False , # Options
      splitTolerance=0 , # Options
      newSplitPlane=True , # Options
      nPlaneReverse=0 , # Options
    )

    GEO = CadToCsg(
        title= 'Input Test' ,
        stepFile= f"{input_step_file.resolve()}" ,
        geometryName= f"{output_filename_stem.resolve()}" ,
        outFormat= ('mcnp', 'openMC_XML'),
        options = my_options
    )

    # set parameters values stored in template dictionary
    for key,value in template.items():
      GEO.set(key, value)

    # deletes the output openmc and mcnp output files if it already exists
    output_filename_stem.with_suffix(".mcnp").unlink(missing_ok=True)
    output_filename_stem.with_suffix(".xml").unlink(missing_ok=True)

    GEO.Start()

    assert output_filename_stem.with_suffix(".mcnp").exists()
    assert output_filename_stem.with_suffix(".xml").exists()
