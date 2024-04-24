from . import additional_files as OutFiles
from .mcnp_format import McnpInput
from .openmc_format import OpenmcInput
from .phits_format import PhitsInput
from .serpent_format import SerpentInput


def write_geometry(
    UniverseBox,
    MetaList,
    Surfaces,
    void_gen,
    cell_comment_file,
    cell_summary_file,
    out_format,
    geometry_name,
    step_file,
    title,
    vol_sdef,
    vol_card,
    ucard,
    dummy_mat,
    mat_file,
    void_mat,
    start_cell,
):

    # write cells comments in file
    if cell_comment_file:
        OutFiles.comments_write(geometry_name, MetaList)
    if cell_summary_file:
        OutFiles.summary_write(geometry_name, MetaList)

    if "mcnp" in out_format:
        mcnpFilename = geometry_name + ".mcnp"
        outBox = (
            UniverseBox.XMin,
            UniverseBox.XMax,
            UniverseBox.YMin,
            UniverseBox.YMax,
            UniverseBox.ZMin,
            UniverseBox.ZMax,
        )
        if void_gen:
            outSphere = (Surfaces["Sph"][-1].Index, Surfaces["Sph"][-1].Surf.Radius)
        else:
            outSphere = None

        mcnpfile = McnpInput(
            Meta=MetaList,
            Surfaces=Surfaces,
            step_file=step_file, # todo remove arg, this is just used if the title is none
            title=title,
            vol_sdef=vol_sdef,
            vol_card=vol_card,
            ucard=ucard,
            dummy_mat=dummy_mat,
        )
        mcnpfile.set_sdef((outSphere, outBox))
        mcnpfile.write_input(mcnpFilename)

    if "openmc_xml" in out_format or "openmc_py" in out_format:
        OMCFile = OpenmcInput(Meta=MetaList, Surfaces=Surfaces)

    if "openmc_xml" in out_format:
        omcFilename = geometry_name + ".xml"
        OMCFile.write_xml(omcFilename)

    if "openmc_py" in out_format:
        omcFilename = geometry_name + ".py"
        OMCFile.write_py(omcFilename)

    if "serpent" in out_format:
        serpentFilename = geometry_name + ".serp"
        outBox = (
            UniverseBox.XMin,
            UniverseBox.XMax,
            UniverseBox.YMin,
            UniverseBox.YMax,
            UniverseBox.ZMin,
            UniverseBox.ZMax,
        )
        if void_gen:
            outSphere = (Surfaces["Sph"][-1].Index, Surfaces["Sph"][-1].Surf.Radius)
        else:
            outSphere = None

        Serpentfile = SerpentInput(
            MetaList, Surfaces, step_file, title, vol_sdef, vol_card, ucard, dummy_mat
        )
        # Serpentfile.set_sdef((outSphere,outBox))
        Serpentfile.write_input(serpentFilename)

    if "phits" in out_format:
        phitsFilename = geometry_name + ".inp"
        phits_outBox = (
            UniverseBox.XMin,
            UniverseBox.XMax,
            UniverseBox.YMin,
            UniverseBox.YMax,
            UniverseBox.ZMin,
            UniverseBox.ZMax,
        )
        if void_gen:
            phits_outSphere = (
                Surfaces["Sph"][-1].Index,
                Surfaces["Sph"][-1].Surf.Radius,
            )
        else:
            phits_outSphere = None

        phitsfile = PhitsInput(
            MetaList,
            Surfaces,
            step_file,
            title,
            vol_sdef,
            vol_card,
            ucard,
            dummy_mat,
            mat_file,
            void_mat,
            start_cell,
        )
        # phitsfile.set_sdef_phits((phits_outSphere,phits_outBox))
        phitsfile.write_phits(phitsFilename)
