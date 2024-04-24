# Init file of GEOUNED module
#

# We load the STEP and the materials

# this try except attempts to import freecad (lowercase) which is the conda
# package name for FreeCAD (mixed case) upon import the conda package appends
# the sys path for Conda installed FreeCAD, consequently FreeCAD can then be
# found by subsequent import statements through out the code base
try:
    import freecad
except:
    pass
import configparser
from datetime import datetime
from os import mkdir, path

import FreeCAD
import Part

from .code_version import *
from .conversion import cell_definition as Conv
from .cuboid.translate import translate
from .decompose import decom_one as Decom
from .load_file import load_step as Load
from .utils import functions as UF
from .utils.boolean_solids import build_c_table_from_solids
from .utils.Options.classes import mcnp_numeric_format, Options, Tolerances
from .void import void as void
from .write.functions import write_mcnp_cell_def
from .write.write_files import write_geometry


class Geouned:

    def __init__(self, title="Geouned conversion"):
        """ """
        self.__dict__ = dict()
        self.__dict__["step_file"] = ""
        self.__dict__["geometry_name"] = ""
        self.__dict__["matFile"] = ""
        self.__dict__["out_format"] = ("mcnp",)
        self.__dict__["title"] = title
        self.__dict__["void_gen"] = True
        self.__dict__["debug"] = False
        self.__dict__["comp_solids"] = True
        self.__dict__["vol_sdef"] = False
        self.__dict__["dummy_mat"] = False
        self.__dict__["vol_card"] = True
        self.__dict__["UCARD"] = None
        self.__dict__["simplify"] = "No"
        self.__dict__["cellRange"] = []
        self.__dict__["exportSolids"] = ""
        self.__dict__["min_void_size"] = 200  # units mm
        self.__dict__["maxSurf"] = 50
        self.__dict__["maxBracket"] = 30
        self.__dict__["voidMat"] = []
        self.__dict__["voidExclude"] = []
        self.__dict__["startCell"] = 1
        self.__dict__["startSurf"] = 1
        self.__dict__["cell_comment_file"] = False
        self.__dict__["cell_summary_file"] = True
        self.__dict__["sort_enclosure"] = False

    def set_options(self):
        toleranceKwrd = (
            "relativeTolerance",
            "relativePrecision",
            "singleValue",
            "generalDistance",
            "generalAngle",
            "planeDistance",
            "planeAngle",
            "cylinderDistance",
            "cylinderAngle",
            "sphereDistance",
            "coneDistance",
            "coneAngle",
            "torusDistance",
            "torusAngle",
            "minArea",
        )
        numericKwrd = (
            "P_abc",
            "P_d",
            "P_xyz",
            "S_r",
            "S_xyz",
            "C_r",
            "C_xyz",
            "K_tan2",
            "K_xyz",
            "T_r",
            "T_xyz",
            "GQ_1to6",
            "GQ_7to9",
            "GQ_10",
        )
        tolKwrdEquiv = {
            "relativeTolerance": "relativeTol",
            "relativePrecision": "relativePrecision",
            "singleValue": "value",
            "generalDistance": "distance",
            "generalAngle": "angle",
            "planeDistance": "pln_distance",
            "planeAngle": "pln_angle",
            "cylinderDistance": "cyl_distance",
            "cylinderAngle": "cyl_angle",
            "sphereDistance": "sph_distance",
            "coneDistance": "kne_distance",
            "coneAngle": "kne_angle",
            "torusDistance": "tor_distance",
            "torusAngle": "tor_angle",
            "minArea": "min_area",
        }

        config = configparser.ConfigParser()
        config.optionxform = str
        config.read(self.__dict__["title"])
        for section in config.sections():
            if section == "Files":
                for key in config["Files"].keys():
                    if key in ("geometry_name", "matFile", "title"):
                        self.set(key, config.get("Files", key))

                    elif key == "step_file":
                        value = config.get("Files", key).strip()
                        lst = value.split()
                        if value[0] in ("(", "[") and value[-1] in ("]", ")"):
                            data = value[1:-1].split(",")
                            data = [x.strip() for x in data]
                            self.set(key, data)
                        elif len(lst) > 1:
                            self.set(key, lst)
                        else:
                            self.set(key, value)

                    elif key == "out_format":
                        raw = config.get("Files", key).strip()
                        values = tuple(x.strip() for x in raw.split(","))
                        out_format = []
                        for v in values:
                            if v.lower() == "mcnp":
                                out_format.append("mcnp")
                            elif v.lower() == "openmc_xml":
                                out_format.append("openMC_XML")
                            elif v.lower() == "openmc_py":
                                out_format.append("openMC_PY")
                            elif v.lower() == "serpent":
                                out_format.append("serpent")
                            elif v.lower() == "phits":
                                out_format.append("phits")
                        self.set(key, tuple(out_format))

            elif section == "parameters":
                for key in config["parameters"].keys():
                    if key in (
                        "void_gen",
                        "debug",
                        "comp_solids",
                        "vol_sdef",
                        "vol_card",
                        "dummy_mat",
                        "cell_summary_file",
                        "cell_comment_file",
                        "sort_enclosure",
                    ):
                        self.set(key, config.getboolean("parameters", key))
                    elif key in (
                        "min_void_size",
                        "maxSurf",
                        "maxBracket",
                        "startCell",
                        "startSurf",
                    ):
                        self.set(key, config.getint("parameters", key))
                    elif key in ("exportSolids", "UCARD", "simplify"):
                        self.set(key, config.get("parameters", key))
                    elif key == "voidMat":
                        value = config.get("parameters", key).strip()
                        data = value[1:-1].split(",")
                        self.set(key, (int(data[0]), float(data[1]), data[2]))
                    else:
                        value = config.get("parameters", key).strip()
                        data = value[1:-1].split(",")
                        self.set(key, tuple(map(int, data)))

            elif section == "options":
                for key in config["options"].keys():
                    if key in (
                        "force_cylinder",
                        "new_split_plane",
                        "delLastNumber",
                        "verbose",
                        "scaleUP",
                        "quadricPY",
                        "Facets",
                        "prnt3PPlane",
                        "forceNoOverlap",
                    ):
                        setattr(Options, key, config.getboolean("options", key))
                    elif key in ("enlargeBox", "n_plane_reverse", "split_tolerance"):
                        setattr(Options, key, config.getfloat("options", key))

            elif section == "Tolerances":
                for key in config["Tolerances"].keys():
                    if key == "relativeTolerance":
                        setattr(Tolerances, key, config.getboolean("Tolerances", key))
                    elif key in toleranceKwrd:
                        setattr(
                            Tolerances,
                            tolKwrdEquiv[key],
                            config.getfloat("Tolerances", key),
                        )

            elif section == "mcnp_Numeric_Format":
                PdEntry = False
                for key in config["mcnp_Numeric_Format"].keys():
                    if key in numericKwrd:
                        setattr(
                            mcnp_numeric_format,
                            key,
                            config.get("mcnp_Numeric_Format", key),
                        )
                    if key == "P_d":
                        PdEntry = True

            else:
                print("bad section name : {}".format(section))

        if self.__dict__["geometry_name"] == "":
            self.__dict__["geometry_name"] = self.__dict__["step_file"][:-4]

        if Options.prnt3PPlane and not PdEntry:
            mcnp_numeric_format.P_d = "22.15e"

        print(self.__dict__)

    def set(self, kwrd, value):

        if kwrd not in self.__dict__.keys():
            print("Bad entry : {}".format(kwrd))
            return

        if kwrd == "step_file":
            if isinstance(value, (list, tuple)):
                for v in value:
                    if not isinstance(v, str):
                        print("elemt in {} list should be string".format(kwrd))
                        return
            elif not isinstance(value, str):
                print("{} should be string or tuple of strings".format(kwrd))
                return

        elif kwrd == "UCARD":
            if value == "None":
                value = None
            elif value.isdigit():
                value = int(value)
            else:
                print("{} value should be None or integer".format(kwrd))
                return
        elif kwrd == "out_format":
            if len(value) == 0:
                return
        elif kwrd in ("geometry_name", "matFile", "exportSolids"):
            if not isinstance(value, str):
                print("{} value should be str instance".format(kwrd))
                return
        elif kwrd in ("cellRange", "voidMat", "voidExclude"):
            if not isinstance(value, (list, tuple)):
                print("{} value should be list or tuple".format(kwrd))
                return
        elif kwrd in ("min_void_size", "maxSurf", "maxBracket", "startCell", "startSurf"):
            if not isinstance(value, int):
                print("{} value should be integer".format(kwrd))
                return
        elif kwrd in (
            "void_gen",
            "debug",
            "comp_solids",
            "simplifyCTable",
            "vol_sdef",
            "vol_card",
            "dummy_mat",
            "cell_summary_file",
            "cell_comment_file",
            "sort_enclosure",
        ):
            if not isinstance(value, bool):
                print("{} value should be boolean".format(kwrd))
                return

        self.__dict__[kwrd] = value
        if kwrd == "step_file" and self.__dict__["geometry_name"] == "":
            if isinstance(value, (tuple, list)):
                self.__dict__["geometry_name"] == "joined_step_files"
            else:
                self.__dict__["geometry_name"] == value[:-4]

    def start(self):

        print("start")
        FreeCAD_Version = "{V[0]:}.{V[1]:}.{V[2]:}".format(V=FreeCAD.Version())
        print(
            "GEOUNED version {} {} \nFreeCAD version {}".format(
                GEOUNED_Version, GEOUNED_ReleaseDate, FreeCAD_Version
            )
        )

        code_setting = self.__dict__
        if code_setting is None:
            raise ValueError("Cannot run the code. Input are missing")
        if code_setting["step_file"] == "":
            raise ValueError("Cannot run the code. Step file name is missing")

        step_file = code_setting["step_file"]
        matfile = code_setting["matFile"]

        if isinstance(step_file, (tuple, list)):
            for stp in step_file:
                if not path.isfile(stp):
                    raise FileNotFoundError(f"Step file {stp} not found.\nStop.")
        else:
            if not path.isfile(step_file):
                raise FileNotFoundError(f"Step file {step_file} not found.\nStop.")

        startTime = datetime.now()

        if isinstance(step_file, (list, tuple)):
            MetaChunk = []
            EnclosureChunk = []
            for stp in step_file:
                print("read step file : {}".format(stp))
                Meta, Enclosure = Load.LoadCAD(stp, matfile)
                MetaChunk.append(Meta)
                EnclosureChunk.append(Enclosure)
            MetaList = join_meta_lists(MetaChunk)
            EnclosureList = join_meta_lists(EnclosureChunk)
        else:
            print("read step file : {}".format(step_file))
            MetaList, EnclosureList = Load.LoadCAD(
                step_file, matfile, code_setting["voidMat"], code_setting["comp_solids"]
            )

        print("End of loading phase")
        tempstr1 = str(datetime.now() - startTime)
        print(tempstr1)
        tempTime = datetime.now()

        # Select a specific solid range from original STEP solids
        if code_setting["cellRange"]:
            MetaList = MetaList[
                code_setting["cellRange"][0] : code_setting["cellRange"][1]
            ]

        # export in STEP format solids read from input file
        # terminate excution
        if code_setting["exportSolids"] != "":
            solids = []
            for m in MetaList:
                if m.IsEnclosure:
                    continue
                solids.extend(m.Solids)
            Part.makeCompound(solids).exportStep(code_setting["exportSolids"])
            msg = (
                f'Solids exported in file {code_setting["exportSolids"]}\n'
                "GEOUNED Finish. No solid translation performed."
            )
            raise ValueError(msg)

        # set up Universe
        if EnclosureList:
            UniverseBox = get_universe(MetaList + EnclosureList)
        else:
            UniverseBox = get_universe(MetaList)
        Comsolids = []

        surfOffset = code_setting["startSurf"] - 1
        Surfaces = UF.Surfaces_dict(offset=surfOffset)

        warnSolids = []
        warnEnclosures = []
        coneInfo = dict()
        tempTime0 = datetime.now()
        if not Options.Facets:

            # decompose all solids in elementary solids (convex ones)
            warningSolidList = decompose_solids(
                MetaList, Surfaces, UniverseBox, code_setting, True
            )

            # decompose Enclosure solids
            if code_setting["void_gen"] and EnclosureList:
                warningEnclosureList = decompose_solids(
                    EnclosureList, Surfaces, UniverseBox, code_setting, False
                )

            print("End of decomposition phase")

            # start Building CGS cells phase

            for j, m in enumerate(MetaList):
                if m.IsEnclosure:
                    continue
                print("Building cell: ", j + 1)
                cones = Conv.cell_def(m, Surfaces, UniverseBox)
                if cones:
                    coneInfo[m.__id__] = cones
                if j in warningSolidList:
                    warnSolids.append(m)
                if not m.Solids:
                    print("none", j, m.__id__)
                    print(m.Definition)

            if Options.forceNoOverlap:
                Conv.no_overlapping_cell(MetaList, Surfaces)

        else:
            translate(MetaList, Surfaces, UniverseBox, code_setting)
            # decompose Enclosure solids
            if code_setting["void_gen"] and EnclosureList:
                warningEnclosureList = decompose_solids(
                    EnclosureList, Surfaces, UniverseBox, code_setting, False
                )

        tempstr2 = str(datetime.now() - tempTime)
        print(tempstr2)

        #  building enclosure solids

        if code_setting["void_gen"] and EnclosureList:
            for j, m in enumerate(EnclosureList):
                print("Building Enclosure Cell: ", j + 1)
                cones = Conv.cell_def(m, Surfaces, UniverseBox)
                if cones:
                    coneInfo[m.__id__] = cones
                if j in warningEnclosureList:
                    warnEnclosures.append(m)

        tempTime1 = datetime.now()

        # void generation phase
        MetaVoid = []
        if code_setting["void_gen"]:
            print("Build Void")
            print(code_setting["voidExclude"])
            if not code_setting["voidExclude"]:
                MetaReduced = MetaList
            else:
                MetaReduced = exclude_cells(MetaList, code_setting["voidExclude"])

            if MetaList:
                init = MetaList[-1].__id__ - len(EnclosureList)
            else:
                init = 0
            MetaVoid = void.void_generation(
                MetaReduced, EnclosureList, Surfaces, UniverseBox, code_setting, init
            )

        # if code_setting['simplify'] == 'full' and not Options.forceNoOverlap:
        if code_setting["simplify"] == "full":
            Surfs = {}
            for lst in Surfaces.values():
                for s in lst:
                    Surfs[s.Index] = s

            for c in MetaList:
                if c.Definition.level == 0 or c.IsEnclosure:
                    continue
                print("simplify cell", c.__id__)
                Box = UF.get_box(c)
                CT = build_c_table_from_solids(Box, (c.Surfaces, Surfs), option="full")
                c.Definition.simplify(CT)
                c.Definition.clean()
                if type(c.Definition.elements) is bool:
                    print(
                        "unexpected constant cell {} :{}".format(
                            c.__id__, c.Definition.elements
                        )
                    )

        tempTime2 = datetime.now()
        print("build Time:", tempTime2 - tempTime1)

        print(datetime.now() - startTime)

        cellOffSet = code_setting["startCell"] - 1
        if EnclosureList and code_setting["sort_enclosure"]:
            # sort group solid cell / void cell sequence in each for each enclosure
            # if a solid belong to several enclosure, its definition will be written
            # for the highest enclosure level or if same enclosure level in the first
            # enclosure found
            MetaList = sort_enclosure(MetaList, MetaVoid, cellOffSet)
        else:
            # remove Null Cell and apply cell numbering offset
            deleted = []
            idLabel = {0: 0}
            icount = cellOffSet
            for i, m in enumerate(MetaList):
                if m.NullCell or m.IsEnclosure:
                    deleted.append(i)
                    continue

                icount += 1
                m.label = icount
                idLabel[m.__id__] = m.label

            for i in reversed(deleted):
                del MetaList[i]

            lineComment = """\
##########################################################
             VOID CELLS
##########################################################"""
            mc = UF.GeounedSolid(None)
            mc.Comments = lineComment
            MetaList.append(mc)

            deleted = []
            for i, m in enumerate(MetaVoid):
                if m.NullCell:
                    deleted.append(i)
                    continue
                icount += 1
                m.label = icount
                update_comment(m, idLabel)
            for i in reversed(deleted):
                del MetaVoid[i]

            MetaList.extend(MetaVoid)

        print_warning_solids(warnSolids, warnEnclosures)

        # add plane definition to cone
        process_cones(MetaList, coneInfo, Surfaces, UniverseBox)

        # write outputformat input
        write_geometry(UniverseBox, MetaList, Surfaces, code_setting)

        print("End of mcnp, OpenMC, Serpent and phits translation phase")

        print("Process finished")
        print(datetime.now() - startTime)

        print("Translation time of solid cells", tempTime1 - tempTime0)
        print("Translation time of void cells", tempTime2 - tempTime1)


def decompose_solids(MetaList, Surfaces, UniverseBox, setting, meta):
    totsolid = len(MetaList)
    warningSolids = []
    for i, m in enumerate(MetaList):
        if meta and m.IsEnclosure:
            continue
        print("Decomposing solid: {}/{} ".format(i + 1, totsolid))
        if setting["debug"]:
            print(m.Comments)
            if not path.exists("debug"):
                mkdir("debug")
            if m.IsEnclosure:
                m.Solids[0].exportStep("debug/origEnclosure_{}.stp".format(i))
            else:
                m.Solids[0].exportStep("debug/origSolid_{}.stp".format(i))

        comsolid, err = Decom.split_solid(Part.makeCompound(m.Solids), UniverseBox)

        if err != 0:
            if not path.exists("Suspicious_solids"):
                mkdir("Suspicious_solids")
            if m.IsEnclosure:
                Part.CompSolid(m.Solids).exportStep(
                    "Suspicious_solids/Enclosure_original_{}.stp".format(i)
                )
                comsolid.exportStep(
                    "Suspicious_solids/Enclosure_split_{}.stp".format(i)
                )
            else:
                Part.CompSolid(m.Solids).exportStep(
                    "Suspicious_solids/Solid_original_{}.stp".format(i)
                )
                comsolid.exportStep("Suspicious_solids/Solid_split_{}.stp".format(i))

            warningSolids.append(i)

        if setting["debug"]:
            if m.IsEnclosure:
                comsolid.exportStep("debug/compEnclosure_{}.stp".format(i))
            else:
                comsolid.exportStep("debug/compSolid_{}.stp".format(i))
        Surfaces.extend(
            Decom.extract_surfaces(comsolid, "All", UniverseBox, MakeObj=True)
        )
        m.setCADSolid()
        m.updateSolids(comsolid.Solids)

    return warningSolids


def update_comment(meta, idLabel):
    if meta.__commentInfo__ is None:
        return
    if meta.__commentInfo__[1] is None:
        return
    newLabel = (idLabel[i] for i in meta.__commentInfo__[1])
    meta.setComments(void.void_commentLine((meta.__commentInfo__[0], newLabel)))


def process_cones(MetaList, coneInfo, Surfaces, UniverseBox):
    cellId = tuple(coneInfo.keys())
    for m in MetaList:
        if m.__id__ not in cellId and not m.Void:
            continue

        if m.Void and m.__commentInfo__ is not None:
            if m.__commentInfo__[1] is None:
                continue
            cones = set()
            for Id in m.__commentInfo__[1]:
                if Id in cellId:
                    cones.update(-x for x in coneInfo[Id])
            Conv.add_cone_plane(m.Definition, cones, Surfaces, UniverseBox)
        elif not m.Void:
            Conv.add_cone_plane(m.Definition, coneInfo[m.__id__], Surfaces, UniverseBox)


def get_universe(MetaList):
    d = 10
    Box = MetaList[0].BoundBox
    xmin = Box.XMin
    xmax = Box.XMax
    ymin = Box.YMin
    ymax = Box.YMax
    zmin = Box.ZMin
    zmax = Box.ZMax
    for m in MetaList[1:]:
        # MIO. This was removed since in HELIAS the enclosure cell is the biggest one
        # if m.IsEnclosure: continue
        xmin = min(m.BoundBox.XMin, xmin)
        xmax = max(m.BoundBox.XMax, xmax)
        ymin = min(m.BoundBox.YMin, ymin)
        ymax = max(m.BoundBox.YMax, ymax)
        zmin = min(m.BoundBox.ZMin, zmin)
        zmax = max(m.BoundBox.ZMax, zmax)

    return FreeCAD.BoundBox(
        FreeCAD.Vector(xmin - d, ymin - d, zmin - d),
        FreeCAD.Vector(xmax + d, ymax + d, zmax + d),
    )


def print_warning_solids(warnSolids, warnEnclosures):

    if warnSolids or warnEnclosures:
        fic = open("Warning_Solids_definition.txt", "w")
    else:
        return

    if warnSolids:
        lines = "Solids :\n"
        for sol in warnSolids:
            lines += "\n"
            lines += "{}\n".format(sol.label)
            lines += "{}\n".format(sol.Comments)
            lines += "{}\n".format(write_mcnp_cell_def(sol.Definition))
        fic.write(lines)

    if warnEnclosures:
        lines = "Enclosures :\n"
        for sol in warnEnclosures:
            lines += "\n"
            lines += "{}\n".format(sol.label)
            lines += "{}\n".format(sol.Comments)
            lines += "{}\n".format(write_mcnp_cell_def(sol.Definition))

        fic.write(lines)

    fic.close()


def join_meta_lists(MList):

    newMetaList = MList[0]
    if MList[0]:
        for M in MList[1:]:
            lastID = newMetaList[-1].__id__ + 1
            for i, meta in enumerate(M):
                meta.__id__ = lastID + i
                newMetaList.append(meta)

    return newMetaList


def exclude_cells(MetaList, labelList):
    voidMeta = []
    for m in MetaList:
        if m.IsEnclosure:
            continue
        found = False
        for label in labelList:
            if label in m.Comments:
                found = True
                break
        if not found:
            voidMeta.append(m)

    return voidMeta


def sort_enclosure(MetaList, MetaVoid, offSet=0):

    newList = {}
    for m in MetaVoid:
        if m.EnclosureID in newList.keys():
            newList[m.EnclosureID].append(m)
        else:
            newList[m.EnclosureID] = [m]

    icount = offSet
    idLabel = {0: 0}
    newMeta = []
    for m in MetaList:
        if m.NullCell:
            continue
        if m.IsEnclosure:
            lineComment = """\
##########################################################
             ENCLOSURE {}
##########################################################""".format(
                m.EnclosureID
            )
            mc = UF.GeounedSolid(None)
            mc.Comments = lineComment
            newMeta.append(mc)
            for e in newList[m.EnclosureID]:
                if e.NullCell:
                    continue
                icount += 1
                e.label = icount
                idLabel[e.__id__] = e.label
                newMeta.append(e)
            lineComment = """\
##########################################################
            END  ENCLOSURE {}
##########################################################""".format(
                m.EnclosureID
            )
            mc = UF.GeounedSolid(None)
            mc.Comments = lineComment
            newMeta.append(mc)

        else:
            icount += 1
            m.label = icount
            idLabel[m.__id__] = m.label
            newMeta.append(m)

    lineComment = """\
##########################################################
             VOID CELLS 
##########################################################"""
    mc = UF.GeounedSolid(None)
    mc.Comments = lineComment
    newMeta.append(mc)

    for v in newList[0]:
        if v.NullCell:
            continue
        icount += 1
        v.label = icount
        idLabel[v.__id__] = v.label
        newMeta.append(v)

    for m in newMeta:
        if not m.Void:
            continue
        if m.IsEnclosure:
            continue
        update_comment(m, idLabel)

    return newMeta
