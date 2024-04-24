import FreeCAD
import Part

from ..load_file import load_functions as LF
from ..utils.basic_functions_part1 import is_opposite
from ..utils.boolean_function import BoolSequence
from ..utils.functions import GeounedSolid, GEOUNED_Surface
from ..utils.options.classes import Options as opt
from . import void_functions as VF
from .void_box import VoidBox


def void_generation(
    MetaList,
    EnclosureList,
    Surfaces,
    UniverseBox,
    sort_enclosure,
    void_mat,
    max_surf,
    max_bracket,
    min_void_size,
    simplify,
    init,
):
    voidList = []

    if EnclosureList:
        NestedEnclosure = LF.setEnclosureLevels(EnclosureList)
        VF.assign_enclosure(MetaList, NestedEnclosure)

        # add to Metalist Level 1 enclosures, remove from list material cells totally embedded in Level 1 enclosures
        newMetaList = VF.select_solids(MetaList, NestedEnclosure[0], UniverseBox)
    else:
        newMetaList = MetaList[:]
        NestedEnclosure = []

    Box = Part.makeBox(
        UniverseBox.XLength,
        UniverseBox.YLength,
        UniverseBox.ZLength,
        FreeCAD.Vector(UniverseBox.XMin, UniverseBox.YMin, UniverseBox.ZMin),
        FreeCAD.Vector(0, 0, 1),
    )

    EnclosureBox = GeounedSolid(None, Box)
    if void_mat:
        EnclosureBox.setMaterial(void_mat[0], void_mat[1], void_mat[2])

    # get voids in 0 Level Enclosure (original Universe)
    # if exist Level 1 enclosures are considered as material cells
    print("Build Void highest enclosure")

    voids = get_void_def(
        newMetaList,
        Surfaces,
        EnclosureBox,
        max_surf,
        max_bracket,
        min_void_size,
        simplify,
        Lev0=True,
    )
    voidList.append(voids)

    # Perform enclosure void
    # Loop until the lowest enclosure level

    for i, Level in enumerate(NestedEnclosure):

        print("Build Void highest enclosure")
        for j, encl in enumerate(Level):
            if encl.CellType == "envelope":
                continue
            newMetaList = VF.select_solids(MetaList, encl.SonEnclosures, encl)
            print("Build Void enclosure {} in enclosure level {}".format(j, i + 1))
            # select solids overlapping current enclosure "encl", and lower level enclosures
            voids = get_void_def(
                newMetaList,
                Surfaces,
                encl,
                max_surf,
                max_bracket,
                min_void_size,
                simplify,
            )
            voidList.append(voids)

    voidList.append(set_graveyard_cell(Surfaces, UniverseBox))

    return VF.update_void_list(init, voidList, NestedEnclosure, sort_enclosure)


def get_void_def(
    MetaList,
    Surfaces,
    Enclosure,
    max_surf,
    max_bracket,
    min_void_size,
    simplify,
    Lev0=False,
):

    if "full" in simplify.lower():
        simplifyVoid = "full"
    elif "void" in simplify.lower():
        simplifyVoid = "diag"
    else:
        simplifyVoid = "no"

    if Lev0:
        Universe = VoidBox(MetaList, Enclosure.BoundBox)
    else:
        Universe = VoidBox(MetaList, Enclosure.CADSolid)

    Initial = [Universe]
    VoidDef = []
    iloop = 0
    while iloop < 50:
        Temp = []
        iloop += 1
        nvoid = len(Initial)
        print("Loop, Box to Split :", iloop, nvoid)

        for iz, z in enumerate(Initial):
            nsurfaces, nbrackets = z.getNumbers()
            if opt.verbose:
                print(
                    "{} {}/{} {} {}".format(iloop, iz + 1, nvoid, nsurfaces, nbrackets)
                )

            if nsurfaces > max_surf and nbrackets > max_bracket:
                newspace = z.split(min_void_size)
            else:
                newspace = None

            if type(newspace) is tuple:
                Temp.extend(newspace)
            else:
                #           if len(z.Objects) >= 50 : z.refine()
                boxDim = (
                    z.BoundBox.XMin * 0.1,
                    z.BoundBox.XMax * 0.1,
                    z.BoundBox.YMin * 0.1,
                    z.BoundBox.YMax * 0.1,
                    z.BoundBox.ZMin * 0.1,
                    z.BoundBox.ZMax * 0.1,
                )

                print("build complementary {} {}".format(iloop, iz))

                cell, CellIn = z.get_void_complementary(Surfaces, simplify=simplifyVoid)
                if cell is not None:
                    VoidCell = (cell, (boxDim, CellIn))
                    VoidDef.append(VoidCell)

        Initial = Temp
        if len(Temp) == 0:
            break

    voidList = []
    for i, vcell in enumerate(VoidDef):
        mVoid = GeounedSolid(i)
        mVoid.Void = True
        mVoid.CellType = "void"
        mVoid.set_definition(vcell[0], simplify=True)
        mVoid.setMaterial(Enclosure.Material, Enclosure.Rho, Enclosure.MatInfo)
        mVoid.setDilution(Enclosure.Dilution)

        mVoid.__commentInfo__ = vcell[1]

        voidList.append(mVoid)

    return voidList


def set_graveyard_cell(Surfaces, UniverseBox):
    Universe = VoidBox([], UniverseBox)

    externalBox = get_universe_complementary(Universe, Surfaces)
    center = UniverseBox.Center
    radius = 0.51 * UniverseBox.DiagonalLength
    sphere = GEOUNED_Surface(("Sphere", (center, radius)), UniverseBox)
    id, exist = Surfaces.addSphere(sphere)

    sphdef = BoolSequence(str(-id))
    sphdef.operator = "AND"
    sphdef.append(externalBox)

    notsph = BoolSequence(str(id))

    mVoidSphIn = GeounedSolid(0)
    mVoidSphIn.Void = True
    mVoidSphIn.CellType = "void"
    mVoidSphIn.set_definition(sphdef)
    mVoidSphIn.setMaterial(0, 0, "Graveyard_in")
    mVoidSphIn.__commentInfo__ = None

    mVoidSphOut = GeounedSolid(1)
    mVoidSphOut.Void = True
    mVoidSphOut.CellType = "void"
    mVoidSphOut.set_definition(notsph)
    mVoidSphOut.setMaterial(0, 0, "Graveyard")
    mVoidSphOut.__commentInfo__ = None

    return (mVoidSphIn, mVoidSphOut)


def get_universe_complementary(Universe, Surfaces):
    Def = BoolSequence(operator="OR")
    for p in Universe.getBoundPlanes():
        id, exist = Surfaces.addPlane(p)
        if not exist:
            Def.elements.append(-id)
        else:
            s = Surfaces.getSurface(id)
            if is_opposite(p.Surf.Axis, s.Surf.Axis):
                Def.elements.append(id)
            else:
                Def.elements.append(-id)
    return Def


def void_commentLine(CellInfo):
    boxDef, cellIn = CellInfo
    cells = ", ".join(map(str, cellIn))
    box = ", ".join(f"{num:.3f}" for num in boxDef)
    line = f"Automatic Generated Void Cell. Enclosure({box})\n"
    line += f"Enclosed cells : ({cells})\n"
    return line
