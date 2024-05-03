import FreeCAD
import Part

from ..LoadFile import LoadFunctions as LF
from ..Utils.BasicFunctions_part1 import is_opposite
from ..Utils.booleanFunction import BoolSequence
from ..Utils.Functions import GeounedSolid, GeounedSurface
from ..Void import voidFunctions as VF
from .VoidBoxClass import VoidBox


def void_generation(
    MetaList,
    EnclosureList,
    Surfaces,
    UniverseBox,
    voidMat,
    maxsurf,
    maxbracket,
    minVoidSize,
    simplify,
    sort_enclosure,
    init,
    enlargeBox,
    nPlaneReverse,
    splitTolerance,
    scale_up,
    verbose,
    pln_distance,
    pln_angle,
    relativeTol,
    sph_distance,
):
    voidList = []

    if EnclosureList:
        NestedEnclosure = LF.set_enclosure_levels(EnclosureList)
        VF.assignEnclosure(MetaList, NestedEnclosure)

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
    if voidMat:
        EnclosureBox.set_material(voidMat[0], voidMat[1], voidMat[2])

    # get voids in 0 Level Enclosure (original Universe)
    # if exist Level 1 enclosures are considered as material cells
    print("Build Void highest enclosure")

    voids = get_void_def(
        MetaList=newMetaList,
        Surfaces=Surfaces,
        Enclosure=EnclosureBox,
        maxsurf=maxsurf,
        maxbracket=maxbracket,
        minVoidSize=minVoidSize,
        simplify=simplify,
        enlargeBox=enlargeBox,
        verbose=verbose,
        nPlaneReverse=nPlaneReverse,
        splitTolerance=splitTolerance,
        scale_up=scale_up,
        pln_distance=pln_distance,
        pln_angle=pln_angle,
        relativeTol=relativeTol,
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
            print(f"Build Void enclosure {j} in enclosure level {i + 1}")
            # select solids overlapping current enclosure "encl", and lower level enclosures
            voids = get_void_def(
                MetaList=newMetaList,
                Surfaces=Surfaces,
                Enclosure=encl,
                maxsurf=maxsurf,
                maxbracket=maxbracket,
                minVoidSize=minVoidSize,
                simplify=simplify,
                enlargeBox=enlargeBox,
                verbose=verbose,
                nPlaneReverse=nPlaneReverse,
                splitTolerance=splitTolerance,
                scale_up=scale_up,
                pln_distance=pln_distance,
                pln_angle=pln_angle,
                relativeTol=relativeTol,
                Lev0=False,
            )
            voidList.append(voids)

    voidList.append(set_graveyard_cell(Surfaces, UniverseBox, 
                                       pln_distance, pln_angle, 
                                       relativeTol, sph_distance))

    return VF.update_void_list(init, voidList, NestedEnclosure, sort_enclosure)


def get_void_def(
    MetaList,
    Surfaces,
    Enclosure,
    maxsurf,
    maxbracket,
    minVoidSize,
    simplify,
    enlargeBox,
    verbose,
    nPlaneReverse,
    splitTolerance,
    scale_up,
    pln_distance,
    pln_angle,
    relativeTol,
    Lev0=False,
):

    # TODO on another PR this can be removed with an attribute setter on the CadToCsg class
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
            nsurfaces, nbrackets = z.get_numbers()
            if verbose:
                print(f"{iloop} {iz + 1}/{nvoid} {nsurfaces} {nbrackets}")

            if nsurfaces > maxsurf and nbrackets > maxbracket:
                newspace = z.split(minVoidSize)
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

                print(f"build complementary {iloop} {iz}")

                cell, CellIn = z.get_void_complementary(
                    Surfaces=Surfaces,
                    enlargeBox=enlargeBox,
                    verbose=verbose,
                    nPlaneReverse=nPlaneReverse,
                    splitTolerance=splitTolerance,
                    scale_up=scale_up,
                    pln_distance=pln_distance,
                    pln_angle=pln_angle,
                    relativeTol=relativeTol,
                    simplify=simplifyVoid,
                )
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
        mVoid.set_definition(vcell[0])
        mVoid.set_material(Enclosure.Material, Enclosure.Rho, Enclosure.MatInfo)
        mVoid.set_dilution(Enclosure.Dilution)

        mVoid.__commentInfo__ = vcell[1]

        voidList.append(mVoid)

    return voidList


def set_graveyard_cell(Surfaces, UniverseBox, pln_distance,
                       pln_angle, relativeTol, sph_distance):
    Universe = VoidBox([], UniverseBox)

    externalBox = get_universe_complementary(Universe, Surfaces, pln_distance, pln_angle, relativeTol)
    center = UniverseBox.Center
    radius = 0.51 * UniverseBox.DiagonalLength
    sphere = GeounedSurface(("Sphere", (center, radius)), UniverseBox)
    id, _ = Surfaces.add_sphere(sphere, sph_distance, relativeTol)

    sphdef = BoolSequence(str(-id))
    sphdef.operator = "AND"
    sphdef.append(externalBox)

    notsph = BoolSequence(str(id))

    mVoidSphIn = GeounedSolid(0)
    mVoidSphIn.Void = True
    mVoidSphIn.CellType = "void"
    mVoidSphIn.set_definition(sphdef)
    mVoidSphIn.set_material(0, 0, "Graveyard_in")
    mVoidSphIn.__commentInfo__ = None

    mVoidSphOut = GeounedSolid(1)
    mVoidSphOut.Void = True
    mVoidSphOut.CellType = "void"
    mVoidSphOut.set_definition(notsph)
    mVoidSphOut.set_material(0, 0, "Graveyard")
    mVoidSphOut.__commentInfo__ = None

    return (mVoidSphIn, mVoidSphOut)


def get_universe_complementary(Universe, Surfaces, pln_distance, pln_angle, relativeTol):
    Def = BoolSequence(operator="OR")
    for p in Universe.get_bound_planes():
        id, exist = Surfaces.add_plane(
            plane=p,
            pln_distance=pln_distance,
            pln_angle=pln_angle,
            relativeTol=relativeTol,
            fuzzy=False,
        )
        if not exist:
            Def.elements.append(-id)
        else:
            s = Surfaces.get_surface(id)
            if is_opposite(p.Surf.Axis, s.Surf.Axis):
                Def.elements.append(id)
            else:
                Def.elements.append(-id)
    return Def


def void_comment_line(CellInfo):
    boxDef, cellIn = CellInfo
    cells = ", ".join(map(str, cellIn))
    box = ", ".join(f"{num:.3f}" for num in boxDef)
    line = f"Automatic Generated Void Cell. Enclosure({box})\n"
    line += f"Enclosed cells : ({cells})\n"
    return line
