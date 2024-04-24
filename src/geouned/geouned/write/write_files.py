from . import additional_files as OutFiles
from .mcnp_format import McnpInput
from .openmc_format import OpenmcInput
from .phits_format import PhitsInput
from .serpent_format import SerpentInput


def write_geometry(UniverseBox, MetaList, Surfaces, code_setting):

    baseName = code_setting["geometry_name"]

    # write cells comments in file
    if code_setting["cell_comment_file"]:
        OutFiles.comments_write(baseName, MetaList)
    if code_setting["cell_summary_file"]:
        OutFiles.summary_write(baseName, MetaList)

    if "mcnp" in code_setting["out_format"]:
        mcnpFilename = baseName + ".mcnp"
        outBox = (
            UniverseBox.XMin,
            UniverseBox.XMax,
            UniverseBox.YMin,
            UniverseBox.YMax,
            UniverseBox.ZMin,
            UniverseBox.ZMax,
        )
        if code_setting["void_gen"]:
            outSphere = (Surfaces["Sph"][-1].Index, Surfaces["Sph"][-1].Surf.Radius)
        else:
            outSphere = None

        mcnpfile = McnpInput(MetaList, Surfaces, code_setting)
        mcnpfile.set_sdef((outSphere, outBox))
        mcnpfile.write_input(mcnpFilename)

    if (
        "openMC_XML" in code_setting["out_format"]
        or "openMC_PY" in code_setting["out_format"]
    ):
        OMCFile = OpenmcInput(MetaList, Surfaces, code_setting)

    if "openMC_XML" in code_setting["out_format"]:
        omcFilename = baseName + ".xml"
        OMCFile.write_xml(omcFilename)

    if "openMC_PY" in code_setting["out_format"]:
        omcFilename = baseName + ".py"
        OMCFile.write_py(omcFilename)

    if "serpent" in code_setting["out_format"]:
        serpentFilename = baseName + ".serp"
        outBox = (
            UniverseBox.XMin,
            UniverseBox.XMax,
            UniverseBox.YMin,
            UniverseBox.YMax,
            UniverseBox.ZMin,
            UniverseBox.ZMax,
        )
        if code_setting["void_gen"]:
            outSphere = (Surfaces["Sph"][-1].Index, Surfaces["Sph"][-1].Surf.Radius)
        else:
            outSphere = None

        Serpentfile = SerpentInput(MetaList, Surfaces, code_setting)
        # Serpentfile.set_sdef((outSphere,outBox))
        Serpentfile.write_input(serpentFilename)

    if "phits" in code_setting["out_format"]:
        phitsFilename = baseName + ".inp"
        phits_outBox = (
            UniverseBox.XMin,
            UniverseBox.XMax,
            UniverseBox.YMin,
            UniverseBox.YMax,
            UniverseBox.ZMin,
            UniverseBox.ZMax,
        )
        if code_setting["void_gen"]:
            phits_outSphere = (
                Surfaces["Sph"][-1].Index,
                Surfaces["Sph"][-1].Surf.Radius,
            )
        else:
            phits_outSphere = None

        phitsfile = PhitsInput(MetaList, Surfaces, code_setting)
        # phitsfile.set_sdef_phits((phits_outSphere,phits_outBox))
        phitsfile.write_phits(phitsFilename)
