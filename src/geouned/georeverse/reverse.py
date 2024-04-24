import FreeCAD
import Import

from .code_version import *
from .Modules.build_cad import build_cad, make_tree
from .Modules.McnpInput import McnpInput
from .Modules.objects import CADCell
from .Modules.process_inp import set_options
from .Modules.xml_input import XmlInput


def reverse(optFile="configRevese.ini"):

    print_code_version()

    setting = set_options(optFile)

    geomfile = setting["fileIn"]
    outname = setting["fileOut"]
    outBox = setting["outBox"]
    inFormat = setting["inFormat"]

    CADselection = {
        "Ustart": setting["UStart"],
        "levelMax": setting["levelMax"],
        "cell": setting["cell"],
        "mat": setting["mat"],
        "format": setting["inFormat"],
    }

    UnivCell = CADCell()
    UnivCell.shape = UnivCell.makeBox(FreeCAD.BoundBox(*outBox))

    # get geometry definition from mcnp input
    if inFormat == "mcnp":
        geom = McnpInput(geomfile)
    elif inFormat == "openMC_XML":
        geom = XmlInput(geomfile)
    else:
        msg = (
            f"input format type {inFormat} is not supported."
            'Supported options are "mcnp" or "openMC_XML"'
        )
        raise ValueError(msg)

    CADCells, fails = build_cad(UnivCell, geom, CADselection)

    if fails:
        print("failed in conversion", fails)

    CADdoc = FreeCAD.newDocument("WorkingDoc")

    make_tree(CADdoc, CADCells)
    Import.export(CADdoc.Objects[0:1], outname + ".stp")
    CADdoc.saveAs(outname + ".FCStd")


def print_code_version():

    FreeCAD_Version = "{V[0]:}.{V[1]:}.{V[2]:}".format(V=FreeCAD.Version())
    title = """\
#########################################################################
#                                                                       # 
#      GEOReverse version {:<11}{}{:>26}
#      FreeCAD    version {:<11}{:>36}  
#                                                                       # 
#########################################################################""".format(
        GEOReverse_Version, GEOReverse_ReleaseDate, "#", FreeCAD_Version, "#"
    )
    print(title)


if __name__ == "__main__":
    reverse()
