import FreeCAD
import Part

from ..decompose import decom_one as Decom
from ..utils import basic_functions_part2 as BF
from ..utils import geometry_gu as GU
from ..utils.basic_functions_part1 import is_opposite, is_parallel
from ..utils.boolean_function import BoolSequence
from ..utils.Options.classes import Options as opt
from ..utils.Options.classes import Tolerances as tol


def commonEdge(face1, face2):
    for e1 in face1.Edges:
        for e2 in face2.Edges:
            if e1.isSame(e2):
                return e1
    return None


def isConvex(face1, face2, edge):
    de = 0.1
    tol = 1.0e-5
    e = edge.Vertexes[1].Point - edge.Vertexes[0].Point
    e.normalize()
    V = e.cross(face2.Surface.Axis)

    P = edge.CenterOfMass + de * V
    if not face2.__face__.isInside(P, tol, True):
        V = -V

    convex = False
    if face1.Surface.Axis.dot(V) < 0:
        if face1.Orientation == "Forward":
            convex = True
    else:
        if face1.Orientation == "Reversed":
            convex = True
    return convex


def removeElement(Faces, idf):
    for i, f in enumerate(Faces):
        if f[0] == idf:
            del Faces[i]
            break


def is_inverted(solid):
    face = solid.Faces[0]
    Range = face.ParameterRange
    u = (Range[1] + Range[0]) / 2.0
    v = (Range[3] + Range[2]) / 2.0

    point2 = face.CenterOfMass.add(face.normalAt(u, v).multiply(1.0e-6))

    if solid.isInside(point2, 1e-7, False):
        return True
    else:
        return False


def get_id(facein, Surfaces):

    if is_parallel(facein.Axis, FreeCAD.Vector(1, 0, 0), tol.pln_angle):
        P = "PX"
    elif is_parallel(facein.Axis, FreeCAD.Vector(0, 1, 0), tol.pln_angle):
        P = "PY"
    elif is_parallel(facein.Axis, FreeCAD.Vector(0, 0, 1), tol.pln_angle):
        P = "PZ"
    else:
        P = "P"

    for s in Surfaces[P]:
        if BF.is_same_plane(
            facein,
            s.Surf,
            dtol=tol.pln_distance,
            atol=tol.pln_angle,
            relTol=tol.relativeTol,
        ):
            return s.Index

    return 0


def translate(MetaList, Surfaces, UniverseBox, setting):
    totsolid = len(MetaList)
    for i, m in enumerate(MetaList):
        if m.IsEnclosure:
            continue
        print("Decomposing solid: {}/{} ".format(i, totsolid))
        if setting["debug"]:
            print(m.Comments)
            if m.IsEnclosure:
                m.Solids[0].exportStep("origEnclosure_{}.stp".format(i))
            else:
                m.Solids[0].exportStep("origSolid_{}.stp".format(i))

        Surfaces.extend(
            Decom.extract_surfaces(
                Part.makeCompound(m.Solids), "Plane3Pts", UniverseBox, MakeObj=False
            )
        )
        set_definition(m, Surfaces)


def set_definition(metaObj, Surfaces):
    solids = metaObj.Solids
    sDef = BoolSequence(operator="OR")

    for sol in solids:
        subSol = BoolSequence(operator="AND")
        flag_inv = is_inverted(sol)
        SolidGu = GU.SolidGu(sol, plane3Pts=True)

        Faces = []
        for face in SolidGu.Faces:
            if abs(face.Area) < 1e-2:
                continue
            if face.Area < 0:
                if opt.verbose:
                    print("Warning : Negative surface Area")
            if str(face.Surface) != "<Plane object>":
                print("Warning : All surfaces must be plane")
                continue
            if face.Orientation not in ("Forward", "Reversed"):
                continue

            id = get_id(face.Surface, Surfaces)
            s = Surfaces.getSurface(id)
            if is_opposite(face.Surface.Axis, s.Surf.Axis, tol.pln_angle):
                id = -id
            if face.Orientation == "Forward":
                id = -id
            if flag_inv:
                id = -id
            Faces.append((id, face))

        while len(Faces) > 0:
            id1, face1 = Faces[0]
            noConvex = []
            for id2, face2 in Faces[1:]:
                edge = commonEdge(face1, face2)
                if edge is None:
                    continue
                if not isConvex(face1, face2, edge):
                    noConvex.append(id2)

            if noConvex != []:
                noConvex.insert(0, id1)
                for i in noConvex:
                    removeElement(Faces, i)
                orPlanes = BoolSequence(operator="OR")
                orPlanes.append(*noConvex)
                subSol.append(orPlanes)
            else:
                removeElement(Faces, id1)
                subSol.append(id1)

        sDef.append(subSol)

    metaObj.set_definition(sDef)
