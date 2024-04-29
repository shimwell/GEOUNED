import FreeCAD
import Import

from .CodeVersion import *
from .Modules.buildCAD import build_cad, makeTree
from .Modules.MCNPinput import McnpInput
from .Modules.Objects import CadCell
from .Modules.processInp import setOptions
from .Modules.XMLinput import XmlInput


def cad_to_csg_from_config(filename="config_cad_to_csg_example.ini"):
    setting = setOptions(filename)

    geom_file = setting["inputFile"]
    out_name = setting["fileOut"]
    out_box = setting["outBox"]
    in_format = setting["inFormat"]
    cad_to_csg()

def cad_to_csg(geom_file, ustart, level_max, cell, mat, in_format, out_box, out_name):

    printCodeVersion()

    CADselection = {
        "Ustart": ustart,
        "levelMax": level_max,
        "cell": cell,
        "mat": mat,
        "format": in_format,
    }

    UnivCell = CadCell()
    UnivCell.shape = UnivCell.makeBox(FreeCAD.BoundBox(*out_box))

    # get geometry definition from MCNP input
    if in_format == "mcnp":
        geom = McnpInput(geom_file)
    elif in_format == "openMC_XML":
        geom = XmlInput(geom_file)
    else:
        msg = (
            f"input format type {in_format} is not supported."
            'Supported options are "mcnp" or "openMC_XML"'
        )
        raise ValueError(msg)

    CADCells, fails = build_cad(UnivCell, geom, CADselection)

    if fails:
        print("failed in conversion", fails)

    CADdoc = FreeCAD.newDocument("WorkingDoc")

    makeTree(CADdoc, CADCells)
    Import.export(CADdoc.Objects[0:1], outname + ".stp")
    CADdoc.saveAs(outname + ".FCStd")


def printCodeVersion():

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
