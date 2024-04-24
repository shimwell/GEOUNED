import BOPTools.SplitAPI

from .build_solid_cell import fuse_solid
from .utils.booleanFunction import BoolSequence


def build_cad(UnivCell, data, config):

    UniverseCut = True
    if "Ustart" not in config.keys():
        config["Ustart"] = 0
    if "levelMax" not in config.keys():
        config["levelMax"] = "all"
    UnivCell.name = 0
    UnivCell.Fill = config["Ustart"]

    # read all surfaces definition
    if config["format"] == "mcnp":
        factor = 10
    else:
        factor = 1

    modelSurfaces = data.get_surfaces(
        scale=factor
    )  # scale change cm in mcnp to mm in CAD Obj

    # read Cells and group into universes
    print(config)
    levels, UniverseCells, modelSurfaces = data.get_filtered_cells(
        modelSurfaces, config
    )

    # assign to each cell the surfaces belonging to the cell
    assign_surface_to_cell(UniverseCells, modelSurfaces)

    #    print(UniverseCells[0][120].definition.str)
    #    print(UniverseCells[0][120].surfaces)
    #    CT=build_c_table_from_solids(UnivCell.shape,UniverseCells[0][70],option='full')
    #    print(CT)
    #    simply = BoolSequence(UniverseCells[0][70].definition.str)
    #    print('antesSimply',simply)
    #    simply.simplify(CT)
    #    print('despues',simply)
    # exit()

    # dictionnary of the cells filled with a given Universe U
    # universeContainers = get_universe_containers(levels,UniverseCells)

    UnivCell.level = None
    levelMax = config["levelMax"]
    Ustart = config["Ustart"]
    if levelMax == "all":
        levelMax = len(levels)

    for lev, Univ in levels.items():
        if Ustart in Univ:
            UnivCell.level = lev - 1
            break
    startInfo = (Ustart, levelMax)

    return build_universe(startInfo, UnivCell, UniverseCells, universeCut=UniverseCut)


def interferencia(container, cell, mode="slice"):

    if mode == "common":
        return cell.shape.common(container.shape)

    Base = cell.shape
    Tool = (container.shape,)

    solids = BOPTools.SplitAPI.slice(Base, Tool, "Split", tolerance=0).Solids
    cellParts = []
    for s in solids:
        if container.shape.isInside(s.CenterOfMass, 0.0, False):
            cellParts.append(s)

    if not cellParts:
        return cell.shape
    else:
        return fuse_solid(cellParts)


def assign_surface_to_cell(UniverseCells, modelSurfaces):

    for Uid, uniCells in UniverseCells.items():
        for cId, c in uniCells.items():
            c.set_surfaces(modelSurfaces)


def get_universe_containers(levels, Universes):
    Ucontainer = {}
    for lev in range(1, len(levels)):
        for U, name in levels[lev]:
            UFILL = Universes[U][name].FILL
            if UFILL in Ucontainer.keys():
                Ucontainer[UFILL].append((U, name, lev))
            else:
                Ucontainer[UFILL] = [(U, name, lev)]
    return Ucontainer


def build_universe(
    startInfo, ContainerCell, AllUniverses, universeCut=True, duplicate=False
):
    CADUniverse = []

    Ustart, levelMax = startInfo
    Universe = AllUniverses[Ustart]

    print(
        "Build Universe {} in container cell {}".format(
            ContainerCell.FILL, ContainerCell.name
        )
    )
    fails = []
    for NTcell in Universe.values():
        if duplicate:
            if NTcell.shape:
                build_shape = False
                if ContainerCell.CurrentTR:
                    cell = NTcell.copy()
                    cell.transform_solid(ContainerCell.CurrentTR)
                else:
                    cell = NTcell
            else:
                CTRF = None
                build_shape = True
                if ContainerCell.CurrentTR:
                    CC = ContainerCell.copy()
                    CC.transform_solid(CC.CurrentTR, reverse=True)
                else:
                    CC = ContainerCell
        else:
            CTRF = ContainerCell.CurrentTR
            CC = ContainerCell
            build_shape = True
            NTcell = NTcell.copy()

        if build_shape:
            print("Level :{}  build Cell {} ".format(CC.level + 1, NTcell.name))
            if type(NTcell.definition) is not BoolSequence:
                NTcell.definition = BoolSequence(NTcell.definition.str)

            bBox = ContainerCell.shape.BoundBox

            debug = False
            if debug:
                NTcell.build_shape(bBox, surfTR=CTRF, simplify=False)
            else:
                try:
                    NTcell.build_shape(bBox, surfTR=CTRF, simplify=False)
                except:
                    print(f"fail converting cell {NTcell.name}")
                    fails.append(NTcell.name)

            if NTcell.shape is None:
                continue

            if duplicate:
                if ContainerCell.CurrentTR:
                    cell = NTcell.copy()
                    cell.transform_solid(ContainerCell.CurrentTR)
                else:
                    cell = NTcell
            else:
                cell = NTcell

        if universeCut and ContainerCell.shape:
            cell.shape = interferencia(ContainerCell, cell)

        if not cell.FILL or ContainerCell.level + 1 == levelMax:
            CADUniverse.append(cell)
        else:
            if ContainerCell.CurrentTR:
                cell.CurrentTR = ContainerCell.CurrentTR.multiply(cell.TRFL)
            cell.level = ContainerCell.level + 1
            univ, ff = build_universe(
                (cell.FILL, levelMax), cell, AllUniverses, universeCut=universeCut
            )
            CADUniverse.append(univ)
            fails.extend(ff)

    return ((ContainerCell.name, Ustart), CADUniverse), fails


def make_tree(CADdoc, CADCells):

    label, universeCADCells = CADCells
    groupObj = CADdoc.addObject("App::Part", "Materials")

    groupObj.Label = "Universe_{}_Container_{}".format(label[1], label[0])

    CADObj = {}
    for i, c in enumerate(universeCADCells):
        if type(c) is tuple:
            groupObj.addObject(make_tree(CADdoc, c))
        else:
            featObj = CADdoc.addObject("Part::FeaturePython", "solid{}".format(i))
            featObj.Label = "Cell_{}_{}".format(c.name, c.MAT)
            featObj.Shape = c.shape
            if c.MAT not in CADObj.keys():
                CADObj[c.MAT] = [featObj]
            else:
                CADObj[c.MAT].append(featObj)

    for mat, matGroup in CADObj.items():
        groupMatObj = CADdoc.addObject("App::Part", "Materials")
        groupMatObj.Label = "Material_{}".format(mat)
        groupMatObj.addObjects(matGroup)
        groupObj.addObject(groupMatObj)

    return groupObj
