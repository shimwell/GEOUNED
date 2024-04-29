import BOPTools.SplitAPI

from .buildSolidCell import fuse_solid
from .Utils.booleanFunction import BoolSequence


def build_cad(univ_cell, data, config):

    if "Ustart" not in config.keys():
        config["Ustart"] = 0
    if "levelMax" not in config.keys():
        config["levelMax"] = "all"
    univ_cell.name = 0
    univ_cell.Fill = config["Ustart"]

    # read all surfaces definition
    if config["format"] == "mcnp":
        factor = 10
    else:
        factor = 1

    modelSurfaces = data.GetSurfaces(
        scale=factor
    )  # scale change cm in mcnp to mm in CAD Obj

    # read Cells and group into universes
    print(config)
    levels, universe_cells, modelSurfaces = data.GetFilteredCells(modelSurfaces, config)

    # assign to each cell the surfaces belonging to the cell
    assign_surface_to_cell(universe_cells, modelSurfaces)

    univ_cell.level = None
    level_max = config["levelMax"]
    u_start = config["Ustart"]
    if level_max == "all":
        level_max = len(levels)

    for lev, Univ in levels.items():
        if u_start in Univ:
            univ_cell.level = lev - 1
            break
    start_info = (u_start, level_max)

    return build_universe(start_info, univ_cell, universe_cells, universe_cut=True)


def interferencia(container, cell, mode="slice"):

    if mode == "common":
        return cell.shape.common(container.shape)

    base = cell.shape
    tool = (container.shape,)

    solids = BOPTools.SplitAPI.slice(base, tool, "Split", tolerance=0).Solids
    cell_parts = []
    for s in solids:
        if container.shape.isInside(s.CenterOfMass, 0.0, False):
            cell_parts.append(s)

    if not cell_parts:
        return cell.shape
    else:
        return fuse_solid(cell_parts)


def assign_surface_to_cell(universe_cells, model_surfaces):

    for _, uniCells in universe_cells.items():
        for _, c in uniCells.items():
            c.setSurfaces(model_surfaces)


def get_universe_containers(levels, universes):
    u_container = {}
    for lev in range(1, len(levels)):
        for U, name in levels[lev]:
            u_fill = universes[U][name].FILL
            if u_fill in u_container.keys():
                u_container[u_fill].append((U, name, lev))
            else:
                u_container[u_fill] = [(U, name, lev)]
    return u_container


def build_universe(
    start_info, container_cell, all_universes, universe_cut=True, duplicate=False
):
    cad_universe = []

    u_start, level_max = start_info
    universe = all_universes[u_start]

    print(
        "Build Universe {} in container cell {}".format(
            container_cell.FILL, container_cell.name
        )
    )
    fails = []
    for nt_cell in universe.values():
        if duplicate:
            if nt_cell.shape:
                build_shape = False
                if container_cell.CurrentTR:
                    cell = nt_cell.copy()
                    cell.transformSolid(container_cell.CurrentTR)
                else:
                    cell = nt_cell
            else:
                CTRF = None
                build_shape = True
                if container_cell.CurrentTR:
                    CC = container_cell.copy()
                    CC.transformSolid(CC.CurrentTR, reverse=True)
                else:
                    CC = container_cell
        else:
            CTRF = container_cell.CurrentTR
            CC = container_cell
            build_shape = True
            nt_cell = nt_cell.copy()

        if build_shape:
            print("Level :{}  build Cell {} ".format(CC.level + 1, nt_cell.name))
            if type(nt_cell.definition) is not BoolSequence:
                nt_cell.definition = BoolSequence(nt_cell.definition.str)

            bBox = container_cell.shape.BoundBox

            debug = False
            if debug:
                nt_cell.buildShape(bBox, surfTR=CTRF, simplify=False)
            else:
                try:
                    nt_cell.buildShape(bBox, surfTR=CTRF, simplify=False)
                except:
                    print(f"fail converting cell {nt_cell.name}")
                    fails.append(nt_cell.name)

            if nt_cell.shape is None:
                continue

            if duplicate:
                if container_cell.CurrentTR:
                    cell = nt_cell.copy()
                    cell.transformSolid(container_cell.CurrentTR)
                else:
                    cell = nt_cell
            else:
                cell = nt_cell

        if universe_cut and container_cell.shape:
            cell.shape = interferencia(container_cell, cell)

        if not cell.FILL or container_cell.level + 1 == level_max:
            cad_universe.append(cell)
        else:
            if container_cell.CurrentTR:
                cell.CurrentTR = container_cell.CurrentTR.multiply(cell.TRFL)
            cell.level = container_cell.level + 1
            univ, ff = build_universe(
                (cell.FILL, level_max), cell, all_universes, universe_cut=universe_cut
            )
            cad_universe.append(univ)
            fails.extend(ff)

    return ((container_cell.name, u_start), cad_universe), fails


def makeTree(cad_doc, cad_cells):

    label, universe_cad_cells = cad_cells
    group_obj = cad_doc.addObject("App::Part", "Materials")

    group_obj.Label = "Universe_{}_Container_{}".format(label[1], label[0])

    cad_obj = {}
    for i, cell in enumerate(universe_cad_cells):
        if type(cell) is tuple:
            group_obj.addObject(makeTree(cad_doc, cell))
        else:
            feat_obj = cad_doc.addObject("Part::FeaturePython", "solid{}".format(i))
            feat_obj.Label = "Cell_{}_{}".format(cell.name, cell.MAT)
            feat_obj.Shape = cell.shape
            if cell.MAT not in cad_obj.keys():
                cad_obj[cell.MAT] = [feat_obj]
            else:
                cad_obj[cell.MAT].append(feat_obj)

    for mat, mat_group in cad_obj.items():
        group_mat_obj = cad_doc.addObject("App::Part", "Materials")
        group_mat_obj.Label = "Material_{}".format(mat)
        group_mat_obj.addObjects(mat_group)
        group_obj.addObject(group_mat_obj)

    return group_obj
