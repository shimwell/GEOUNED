############################
# Module for Cell definiton #
#############################
import math

import FreeCAD
import Part

from ..utils import basic_functions_part2 as BF
from ..utils import functions as UF
from ..utils import geometry_gu as GU
from ..utils.basic_functions_part1 import (
    is_in_line,
    is_opposite,
    is_parallel,
    is_same_value,
    sign_plane,
)
from ..utils.boolean_function import BoolSequence, insert_in_sequence
from ..utils.boolean_solids import build_c_table_from_solids, remove_extra_surfaces
from ..utils.functions import GEOUNED_Surface
from ..utils.options.classes import Options as opt
from ..utils.options.classes import Tolerances as tol


def get_id(facein, Surfaces):

    surfin = str(facein)
    if surfin == "<Plane object>":
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

    elif surfin == "<Cylinder object>":
        for s in Surfaces["Cyl"]:
            if BF.is_same_cylinder(
                facein,
                s.Surf,
                dtol=tol.cyl_distance,
                atol=tol.cyl_angle,
                relTol=tol.relativeTol,
            ):
                return s.Index

    elif surfin == "<Cone object>":
        for s in Surfaces["Cone"]:
            if BF.is_same_cone(
                facein,
                s.Surf,
                dtol=tol.kne_distance,
                atol=tol.kne_angle,
                relTol=tol.relativeTol,
            ):
                return s.Index

    elif surfin[0:6] == "Sphere":
        for s in Surfaces["Sph"]:
            if BF.is_same_sphere(
                facein, s.Surf, tol.sph_distance, relTol=tol.relativeTol
            ):
                return s.Index

    elif surfin == "<Toroid object>":
        for s in Surfaces["Tor"]:
            if BF.is_same_torus(
                facein,
                s.Surf,
                dtol=tol.tor_distance,
                atol=tol.tor_angle,
                relTol=tol.relativeTol,
            ):
                return s.Index

    return 0


def is_inverted(solid):

    face = solid.Faces[0]

    # u=(face.Surface.bounds()[0]+face.Surface.bounds()[1])/2.0 # entre 0 y 2pi si es completo
    # v=face.Surface.bounds()[0]+(face.Surface.bounds()[3]-face.Surface.bounds()[2])/3.0 # a lo largo del eje
    Range = face.ParameterRange
    u = (Range[1] + Range[0]) / 2.0
    v = (Range[3] + Range[2]) / 2.0

    if str(face.Surface) == "<Cylinder object>":
        dist1 = face.Surface.value(u, v).distanceToLine(
            face.Surface.Center, face.Surface.Axis
        )
        dist2 = (
            face.Surface.value(u, v)
            .add(face.Surface.normal(u, v).multiply(1.0e-6))
            .distanceToLine(face.Surface.Center, face.Surface.Axis)
        )
        if (dist2 - dist1) < 0.0:
            # The normal of the cylinder is going inside
            return True
    elif str(face.Surface) == "<Cone object>":
        dist1 = face.Surface.value(u, v).distanceToLine(
            face.Surface.Apex, face.Surface.Axis
        )
        dist2 = (
            face.Surface.value(u, v)
            .add(face.Surface.normal(u, v).multiply(1.0e-6))
            .distanceToLine(face.Surface.Apex, face.Surface.Axis)
        )
        if (dist2 - dist1) < 0.0:
            # The normal of the cylinder is going inside
            return True
    # MIO
    elif str(face.Surface)[0:6] == "Sphere":
        # radii = point - center
        radii = face.Surface.value(u, v).add(face.Surface.Center.multiply(-1))
        radiiB = (
            face.Surface.value(u, v)
            .add(face.Surface.normal(u, v).multiply(1.0e-6))
            .add(face.Surface.Center.multiply(-1))
        )
        # radiiB  = radii.add( face.Surface.normal(u,v).multiply(1.0e-6) )
        if (radiiB.Length - radii.Length) < 0.0:
            # An increasing of the radii vector in the normal direction decreases the radii: oposite normal direction
            return True

    elif str(face.Surface) == "<Plane object>":
        dist1 = face.CenterOfMass.distanceToPoint(solid.BoundBox.Center)
        dist2 = face.CenterOfMass.add(
            face.normalAt(u, v).multiply(1.0e-6)
        ).distanceToPoint(solid.BoundBox.Center)
        point2 = face.CenterOfMass.add(face.normalAt(u, v).multiply(1.0e-6))
        if solid.isInside(point2, 1e-7, False):
            return True

    return False


def gen_plane(face, solid):
    """Generate an additional plane when convex surfaces of second order are presented as a face of the solid"""

    surf = face.Surface
    if str(surf) == "<Cylinder object>":
        return gen_plane_cylinder(face, solid)
    if str(surf) == "<Cone object>":
        return gen_plane_cone(face, solid)
    if str(surf) == "Sphere":
        return gen_plane_sphere(face, solid)


def get_closed_ranges(solid, face_index):

    UNodes = []
    for index in face_index:
        URange = solid.Faces[index].ParameterRange
        UNodes.append((URange[0], index))
        UNodes.append((URange[1], index))
    UNodes.sort()

    closedRange = get_intervals(UNodes)

    aMin = closedRange[0][0][0]
    aMax = closedRange[-1][1][0]

    if abs(aMax - aMin - 2.0 * math.pi) < 1e-2:
        if len(closedRange) == 1:
            closedFace = True
        else:
            endPoint = (closedRange[-1][0][0] - 2 * math.pi, closedRange[-1][0][1])
            closedRange[0][0] = endPoint
            closedRange[0][2].update(closedRange[-1][2])
            del closedRange[-1]

            if len(closedRange) == 1:
                if (
                    abs(closedRange[0][1][0] - closedRange[0][0][0] - 2.0 * math.pi)
                    < 1e-2
                ):
                    closedFace = True
                else:
                    closedFace = False
            else:
                closedFace = False
    else:
        closedFace = False
    return closedRange, closedFace


def get_intervals(UNodes):
    closedRange = []
    posMin = dict()
    posMax = dict()
    for i, node in enumerate(UNodes):
        if node[1] not in posMin.keys():
            posMin[node[1]] = i
        else:
            posMax[node[1]] = i

    UMin = UNodes[0]
    iPos = posMax[UMin[1]]

    while True:
        x = UNodes[iPos]
        end = True
        for i in range(iPos + 1, len(UNodes)):
            nxtInt = UNodes[i][1]
            if (
                UNodes[posMin[nxtInt]][0] - x[0]
            ) < 1e-5:  # x pos is > min boundary of the next inteval inside precision 1e-5
                iPos = posMax[nxtInt]
                end = False
                break

        if end:
            UMax = x
            closedRange.append([UMin, UMax])
            iPos += 1
            if iPos < len(UNodes):
                UMin = UNodes[iPos]
                iPos = posMax[UMin[1]]
            else:
                break

    for rnge in closedRange:
        index = set()
        xmin = rnge[0][0]
        xmax = rnge[1][0]
        for interval in UNodes:
            x = interval[0]
            if (xmin - x) < 1.0e-5 and (x - xmax) < 1.0e-5:
                index.add(interval[1])
        rnge.append(index)

    return closedRange


def get_u_value_boundary(solid, face_index, myIndex):

    faceURange, closedFace = get_closed_ranges(solid, face_index)
    if closedFace:
        return None, None

    for rnge in faceURange:
        if myIndex in rnge[2]:
            UMin, UMax = rnge[0:2]
            return UMin, UMax


def gen_plane_sphere(face, solid):
    Same_Faces = []
    Same_Faces.append(face)

    for f in solid.Faces:
        if f.is_equal(face) or str(f.Surface) != "Sphere":
            continue
        if (
            f.Surface.Center == face.Surface.Center
            and f.Surface.Radius == face.Surface.Radius
        ):
            # print 'Warning: coincident sphere faces are the same'
            for f2 in Same_Faces:
                if f.dist_to_shape(f2)[0] < 1e-6:
                    Same_Faces.append(f)
                    break

    # print Same_Faces
    normal = FreeCAD.Vector(0, 0, 0)
    for f in Same_Faces:
        normal += f.Area * (f.CenterOfMass - face.Surface.Center)

    return Part.Plane(face.Surface.Center, normal).toShape()


def gen_plane_cylinder(face, solid):

    Surf = face.Surface
    rad = Surf.Radius

    if str(Surf) != "<Cylinder object>":
        return None

    myIndex = solid.Faces.index(face)
    face_index = [myIndex]

    for i, face2 in enumerate(solid.Faces):
        if face2.Area < tol.min_area:
            if opt.verbose:
                print(
                    f"Warning: {str(Surf)} surface removed from cell definition. Face area < Min area ({face2.Area} < {tol.min_area}) "
                )
            continue
        if str(face2.Surface) == "<Cylinder object>" and not (face2.is_equal(face)):
            if (
                face2.Surface.Axis.is_equal(face.Surface.Axis, 1e-5)
                and face2.Surface.Radius == rad
                and is_in_line(
                    face2.Surface.Center, face.Surface.Axis, face.Surface.Center
                )
            ):
                # print 'Warning: coincident cylinder faces are the same'
                face_index.append(i)

    UMin, UMax = get_u_value_boundary(solid, face_index, myIndex)
    if UMin is None:
        return None

    U1, i1 = UMin
    U2, i2 = UMax

    V1 = solid.Faces[i1].ParameterRange[2]
    V2 = solid.Faces[i2].ParameterRange[2]

    P1 = solid.Faces[i1].valueAt(U1, V1)
    P2 = solid.Faces[i2].valueAt(U2, V2)

    if P1.is_equal(P2, 1e-5):
        if opt.verbose:
            print("Error in the additional place definition")
        return None

    normal = P2.sub(P1).cross(face.Surface.Axis)
    plane = Part.Plane(P1, normal).toShape()

    return plane


def gen_plane_cylinder_old(face, solid):

    Surf = face.Surface
    rad = Surf.Radius

    if str(Surf) != "<Cylinder object>":
        return None

    face_index = [solid.Faces.index(face)]

    for i, face2 in enumerate(solid.Faces):
        if face2.Area < tol.min_area:
            if opt.verbose:
                print(
                    f"Warning: {str(Surf)} surface removed from cell definition. Face area < Min area ({face2.Area} < {tol.min_area}) "
                )
            continue
        if str(face2.Surface) == "<Cylinder object>" and not (face2.is_equal(face)):
            if (
                face2.Surface.Axis.is_equal(face.Surface.Axis, 1e-5)
                and face2.Surface.Radius == rad
                and is_in_line(
                    face2.Surface.Center, face.Surface.Axis, face.Surface.Center
                )
            ):
                # print 'Warning: coincident cylinder faces are the same'
                face_index.append(i)

    AngleRange = 0.0
    Uval = []
    for index in face_index:
        Range = solid.Faces[index].ParameterRange
        AngleRange = AngleRange + abs(Range[1] - Range[0])
        if not (Range[0] in Uval) and not (Range[1] in Uval):
            Uval.append(Range[0])
            Uval.append(Range[1])
    if 2.0 * math.pi - AngleRange < 1e-2:
        return None

    UVNodes = []
    for index in face_index:
        face2 = solid.Faces[index]
        try:
            UVNodes.append(face2.get_uv_nodes())
        except RuntimeError:
            tess = face.tessellate(1.0, True)
            UVNodes.append(face2.get_uv_nodes())

    Uval_str_cl = []
    for i, elem1 in enumerate(Uval):
        num_str1 = "%11.4E" % elem1
        if abs(elem1) < 1.0e-5:
            num_str1 = "%11.4E" % 0.0
        if not (BF.is_duplicate_in_list(num_str1, i, Uval)):
            Uval_str_cl.append(num_str1)

    face_index_2 = [face_index[0], face_index[0]]

    Node_min = UVNodes[0][0]
    Node_max = UVNodes[0][1]

    dif1_0 = abs(float(Uval_str_cl[0]) - Node_min[0])
    dif2_0 = abs(float(Uval_str_cl[1]) - Node_max[0])

    # searching for minimum and maximum angle points
    for j, Nodes in enumerate(UVNodes):
        for elem in Nodes:
            val = "%11.4E" % elem[0]
            dif1 = abs(float(Uval_str_cl[0]) - elem[0])
            dif2 = abs(float(Uval_str_cl[1]) - elem[0])

            if abs(elem[0]) < 1.0e-5:
                val = "%11.4E" % 0.0
            if dif1 < dif1_0:
                Node_min = elem
                face_index_2[0] = face_index[j]
                dif1_0 = dif1
            if dif2 < dif2_0:
                Node_max = elem
                face_index_2[1] = face_index[j]
                dif2_0 = dif2

    V1 = solid.Faces[face_index_2[0]].valueAt(Node_min[0], Node_min[1])
    V2 = solid.Faces[face_index_2[1]].valueAt(Node_max[0], Node_max[1])

    if V1.is_equal(V2, 1e-5):
        if opt.verbose:
            print("Error in the additional place definition")
        return None

    normal = V2.sub(V1).cross(face.Surface.Axis)
    plane = Part.Plane(V1, normal).toShape()

    return plane


def gen_plane_cone(face, solid):

    Surf = face.Surface
    if str(Surf) != "<Cone object>":
        return None

    myIndex = solid.Faces.index(face)
    face_index = [myIndex]

    for i, face2 in enumerate(solid.Faces):
        if face2.Area < tol.min_area:
            if opt.verbose:
                print(
                    f"Warning: {str(Surf)} surface removed from cell definition. Face area < Min area ({face2.Area} < {tol.min_area}) "
                )
            continue
        if str(face2.Surface) == "<Cone object>" and not (face2.is_equal(face)):
            if (
                face2.Surface.Axis.is_equal(face.Surface.Axis, 1e-5)
                and face2.Surface.Apex.is_equal(face.Surface.Apex, 1e-5)
                and (face2.Surface.SemiAngle - face.Surface.SemiAngle) < 1e-6
            ):
                face_index.append(i)

    UMin, UMax = get_u_value_boundary(solid, face_index, myIndex)
    if UMin is None:
        return None

    U1, i1 = UMin
    U2, i2 = UMax

    V1 = solid.Faces[i1].ParameterRange[2]
    V2 = solid.Faces[i2].ParameterRange[2]

    P1 = solid.Faces[i1].valueAt(U1, V1)
    P2 = solid.Faces[i2].valueAt(U2, V2)

    if P1.is_equal(P2, 1e-5):
        if opt.verbose:
            print("Error in the additional place definition")
        return None

    plane = Part.Plane(P1, P2, face.Surface.Apex).toShape()

    return plane


def gen_plane_cone_old(face, solid):

    Surf = face.Surface
    rad = Surf.Radius
    Axis = face.Surface.Axis
    if str(Surf) != "<Cone object>":
        return None

    face_index = [solid.Faces.index(face)]

    for i, face2 in enumerate(solid.Faces):
        if face2.Area < tol.min_area:
            if opt.verbose:
                print(
                    f"Warning: {str(Surf)} surface removed from cell definition. Face area < Min area ({face2.Area} < {tol.min_area}) "
                )
            continue
        if str(face2.Surface) == "<Cone object>" and not (face2.is_equal(face)):
            if (
                face2.Surface.Axis.is_equal(face.Surface.Axis, 1e-5)
                and face2.Surface.Apex.is_equal(face.Surface.Apex, 1e-5)
                and (face2.Surface.SemiAngle - face.Surface.SemiAngle) < 1e-6
            ):
                face_index.append(i)

    AngleRange = 0.0
    Uval = []
    for index in face_index:
        Range = solid.Faces[index].ParameterRange
        AngleRange = AngleRange + abs(Range[1] - Range[0])
        Uval.append(Range[0])
        Uval.append(Range[1])

    if 2.0 * math.pi - AngleRange < 1e-2:
        return None

    UVNodes = []
    for index in face_index:
        face2 = solid.Faces[index]
        try:
            UVNodes.append(face2.get_uv_nodes())
        except RuntimeError:
            face.tessellate(1.0, True)
            UVNodes.append(face2.get_uv_nodes())

    Uval_str_cl = []

    for i, elem1 in enumerate(Uval):
        num_str1 = "%11.4E" % elem1
        if abs(elem1) < 1.0e-5:
            num_str1 = "%11.4E" % 0.0
        if not (BF.is_duplicate_in_list(num_str1, i, Uval)):
            Uval_str_cl.append(num_str1)

    face_index_2 = [face_index[0], face_index[0]]

    Node_min = UVNodes[0][0]
    Node_max = UVNodes[0][1]
    dif1_0 = abs(float(Uval_str_cl[0]) - Node_min[0])
    dif2_0 = abs(float(Uval_str_cl[1]) - Node_max[0])

    # searching for minimum and maximum angle points
    for j, Nodes in enumerate(UVNodes):
        for elem in Nodes:
            val = "%11.4E" % elem[0]
            dif1 = abs(float(Uval_str_cl[0]) - elem[0])
            dif2 = abs(float(Uval_str_cl[1]) - elem[0])
            if abs(elem[0]) < 1.0e-5:
                val = "%11.4E" % 0.0
            if dif1 < dif1_0:
                Node_min = elem
                face_index_2[0] = face_index[j]
                dif1_0 = dif1
            if dif2 < dif2_0:
                Node_max = elem
                face_index_2[1] = face_index[j]
                dif2_0 = dif2

    V1 = solid.Faces[face_index_2[0]].valueAt(Node_min[0], Node_min[1])
    V2 = solid.Faces[face_index_2[1]].valueAt(Node_max[0], Node_max[1])

    if V1.is_equal(V2, 1e-5):
        if opt.verbose:
            print("Error in the additional place definition")
        return None

    plane = Part.Plane(V1, V2, face.Surface.Apex).toShape()

    return plane


def gen_torus_annex_u_planes(face, Uparams):

    if is_parallel(face.Surface.Axis, FreeCAD.Vector(1, 0, 0), tol.tor_angle):
        axis = FreeCAD.Vector(1, 0, 0)
    elif is_parallel(face.Surface.Axis, FreeCAD.Vector(0, 1, 0), tol.tor_angle):
        axis = FreeCAD.Vector(0, 1, 0)
    elif is_parallel(face.Surface.Axis, FreeCAD.Vector(0, 0, 1), tol.tor_angle):
        axis = FreeCAD.Vector(0, 0, 1)

    center = face.Surface.Center
    p1 = face.valueAt(Uparams[0], 0.0)
    p2 = face.valueAt(Uparams[1], 0.0)
    pmid = face.valueAt(0.5 * (Uparams[0] + Uparams[1]), 0.0)

    if is_same_value(abs(Uparams[1] - Uparams[0]), math.pi, tol.value):
        d = axis.cross(p2 - p1)
        d.normalize()
        if d.dot(pmid - center) < 0:
            d = -d
        return (
            (center, d, face.Surface.MajorRadius, face.Surface.MajorRadius),
            None,
        ), False

    elif Uparams[1] - Uparams[0] < math.pi:
        d = axis.cross(p2 - p1)
        d.normalize()
        if d.dot(pmid - center) < 0:
            d = -d
        return (
            (center, d, face.Surface.MajorRadius, face.Surface.MajorRadius),
            None,
        ), False

    else:
        d1 = axis.cross(p1)
        d1.normalize()
        if d1.dot(pmid - center) < 0:
            d1 = -d1

        d2 = axis.cross(p2)
        d2.normalize()
        if d2.dot(pmid - center) < 0:
            d2 = -d2

        return (
            (center, d1, face.Surface.MajorRadius, face.Surface.MajorRadius),
            (center, d2, face.Surface.MajorRadius, face.Surface.MajorRadius),
        ), True  # (d1 : d2)


def gen_torus_annex_u_planes_org(face, Uparams):

    if is_parallel(face.Surface.Axis, FreeCAD.Vector(1, 0, 0), tol.tor_angle):
        axis = FreeCAD.Vector(1, 0, 0)
    elif is_parallel(face.Surface.Axis, FreeCAD.Vector(0, 1, 0), tol.tor_angle):
        axis = FreeCAD.Vector(0, 1, 0)
    elif is_parallel(face.Surface.Axis, FreeCAD.Vector(0, 0, 1), tol.tor_angle):
        axis = FreeCAD.Vector(0, 0, 1)

    center = face.Surface.Center
    p1 = face.valueAt(Uparams[0], 0.0)
    p2 = face.valueAt(Uparams[1], 0.0)
    pmid = face.valueAt(0.5 * (Uparams[0] + Uparams[1]), 0.0)

    if is_same_value(abs(Uparams[1] - Uparams[0]), math.pi, tol.value):
        d = axis.cross(p2 - p1)
        d.normalize()
        if pmid.dot(d) < 0:
            d = -d
        return ((center, d, face.Surface.MajorRadius), None), False

    else:
        d1 = axis.cross(p1)
        d1.normalize()
        if pmid.dot(d1) < 0:
            d1 = -d1

        d2 = axis.cross(p2)
        d2.normalize()
        if pmid.dot(d2) < 0:
            d2 = -d2

        if Uparams[1] - Uparams[0] < math.pi:
            return (
                (center, d1, face.Surface.MajorRadius, face.Surface.MajorRadius),
                (center, d2, face.Surface.MajorRadius, face.Surface.MajorRadius),
            ), False  # ( d1 d2 )
        else:
            return (
                (center, d1, face.Surface.MajorRadius, face.Surface.MajorRadius),
                (center, d2, face.Surface.MajorRadius, face.Surface.MajorRadius),
            ), True  # (d1 : d2)


def gen_torus_annex_v_surface(face, Vparams, force_cylinder=False):
    if is_parallel(face.Surface.Axis, FreeCAD.Vector(1, 0, 0), tol.tor_angle):
        axis = FreeCAD.Vector(1, 0, 0)
    elif is_parallel(face.Surface.Axis, FreeCAD.Vector(0, 1, 0), tol.tor_angle):
        axis = FreeCAD.Vector(0, 1, 0)
    elif is_parallel(face.Surface.Axis, FreeCAD.Vector(0, 0, 1), tol.tor_angle):
        axis = FreeCAD.Vector(0, 0, 1)

    p1 = face.valueAt(0.0, Vparams[0]) - face.Surface.Center
    z1 = p1.dot(axis)
    d1 = p1.cross(axis).Length

    p2 = face.valueAt(0.0, Vparams[1]) - face.Surface.Center
    z2 = p2.dot(axis)
    d2 = p2.cross(axis).Length

    if is_same_value(z1, z2, tol.distance):
        surfType = "Plane"
        center = face.Surface.Center + z1 * axis
        Vmid = (Vparams[0] + Vparams[1]) * 0.5
        pMid = face.valueAt(0, Vmid) - face.Surface.Center
        if pMid.dot(axis) < z1:
            inSurf = True
        else:
            inSurf = False
        return (
            (center, axis, face.Surface.MajorRadius, face.Surface.MajorRadius),
            surfType,
            inSurf,
        )

    elif is_same_value(d1, d2, tol.distance) or force_cylinder:
        surfType = "Cylinder"
        radius = min(d1, d2)
        center = face.Surface.Center
        if is_same_value(d1, face.Surface.MajorRadius, tol.distance):
            Vmid = (Vparams[0] + Vparams[1]) * 0.5
            pMid = face.valueAt(0, Vmid) - center
            if pMid.cross(axis).Length < face.Surface.MajorRadius:
                inSurf = True
            else:
                inSurf = False
        else:
            if d1 < face.Surface.MajorRadius:
                inSurf = True
            else:
                inSurf = False
        return (center, axis, radius, face.Surface.MinorRadius), surfType, inSurf

    else:
        surfType = "Cone"
        za = (z2 * d1 - z1 * d2) / (d1 - d2)
        Apex = face.Surface.Center + za * axis
        semiAngle = abs(math.atan(d1 / (z1 - za)))

        ConeAxis = axis if za < 0 else -axis

        Vmid = (Vparams[0] + Vparams[1]) * 0.5
        pMid = face.valueAt(0, Vmid) - face.Surface.Center
        zMid = pMid.dot(axis)
        dMid = pMid.cross(axis).Length

        dCone = d1 * (zMid - za) / (z1 - za)
        inSurf = True if dMid < dCone else False

        return (
            (
                Apex,
                ConeAxis,
                semiAngle,
                face.Surface.MinorRadius,
                face.Surface.MajorRadius,
            ),
            surfType,
            inSurf,
        )


def cell_def(metaObj, Surfaces, UniverseBox):

    solids = metaObj.Solids
    delList = []

    PieceDef = BoolSequence(operator="OR")
    PieceObj = []
    cones = set()
    for isol, solid in enumerate(solids):
        SurfPiece = []
        SurfObj = []
        extraPlaneReverse = dict()

        flag_inv = is_inverted(solid)
        SolidGu = GU.SolidGu(solid)
        lastTorus = -1
        for iface, face in enumerate(SolidGu.Faces):
            surfaceType = str(face.Surface)
            if abs(face.Area) < tol.min_area:
                if opt.verbose:
                    print(
                        f"Warning: {surfaceType} surface removed from cell definition. Face area < Min area ({face.Area} < {tol.min_area}) "
                    )
                continue
            if face.Area < 0:
                if opt.verbose:
                    print("Warning : Negative surface Area")
            if face.Orientation not in ("Forward", "Reversed"):
                continue
            if flag_inv:
                orient_temp = face.Orientation
                if orient_temp == "Forward":
                    orient = "Reversed"
                elif orient_temp == "Reversed":
                    orient = "Forward"
            else:
                orient = face.Orientation

            if "Sphere" in surfaceType:
                surfaceType = "Sphere"

            # cone additional plane is added afterward
            if (
                surfaceType in ("<Cylinder object>", "<Cone object>", "Sphere")
                and orient == "Reversed"
            ):
                # cone additional plane is added afterward
                idFace = get_id(face.Surface, Surfaces)
                if surfaceType == "<Cone object>":
                    cones.add(idFace)
                if str(idFace) not in SurfPiece:
                    SurfPiece.append(str(idFace))
                    SurfObj.append(face)

                try:
                    plane = gen_plane(face, SolidGu)
                    if plane is not None:
                        plane = GU.PlaneGu(plane)
                except:
                    plane = None
                    if opt.verbose:
                        print("Warning: generation of additional plane has failed")

                if plane is not None:
                    p = GEOUNED_Surface(
                        ("Plane", (plane.Position, plane.Axis, plane.dim1, plane.dim2)),
                        UniverseBox,
                        Face="Build",
                    )

                    id, exist = Surfaces.addPlane(p)
                    sign = sign_plane(face.CenterOfMass, p)
                    if exist:
                        pp = Surfaces.getSurface(id)
                        if is_opposite(p.Surf.Axis, pp.Surf.Axis, tol.angle):
                            id = -id
                    id *= sign

                    if idFace not in extraPlaneReverse.keys():
                        extraPlaneReverse[idFace] = [str(id)]
                        SurfObj.append(p.shape)
                    else:
                        if str(id) not in extraPlaneReverse[idFace]:
                            extraPlaneReverse[idFace].append(str(id))
                            SurfObj.append(p.shape)

            elif surfaceType == "<Toroid object>":

                if (
                    is_parallel(face.Surface.Axis, FreeCAD.Vector(1, 0, 0), tol.angle)
                    or is_parallel(
                        face.Surface.Axis, FreeCAD.Vector(0, 1, 0), tol.angle
                    )
                    or is_parallel(
                        face.Surface.Axis, FreeCAD.Vector(0, 0, 1), tol.angle
                    )
                ):

                    idT = get_id(face.Surface, Surfaces)

                    index, Uparams = SolidGu.TorusUParams[iface]
                    if index == lastTorus:
                        continue
                    lastTorus = index

                    # add if necesary additional planes following U variable
                    UClosed, UminMax = Uparams
                    # UClosed = True
                    if not UClosed:
                        planes, ORop = gen_torus_annex_u_planes(face, UminMax)
                        plane1, plane2 = planes

                        plane = GEOUNED_Surface(
                            ("Plane", plane1), UniverseBox, Face="Build"
                        )
                        id1, exist = Surfaces.addPlane(plane)
                        if exist:
                            p = Surfaces.getSurface(id1)
                            if is_opposite(plane.Surf.Axis, p.Surf.Axis, tol.pln_angle):
                                id1 = -id1

                        if plane2 is None:
                            UVar = "%i" % id1
                        else:
                            plane = GEOUNED_Surface(
                                ("Plane", plane2), UniverseBox, Face="Build"
                            )
                            id2, exist = Surfaces.addPlane(plane)
                            if exist:
                                p = Surfaces.getSurface(id2)
                                if is_opposite(
                                    plane.Surf.Axis, p.Surf.Axis, tol.pln_angle
                                ):
                                    id2 = -id2

                            UVar = (
                                "(%i : %i)" % (id1, id2)
                                if ORop
                                else "%i %i" % (id1, id2)
                            )

                    else:
                        UVar = ""

                    # add if necesary additional surface following V variable
                    if orient == "Forward":
                        VVar = "-%i" % idT

                    else:
                        index, Vparams = SolidGu.TorusVParams[iface]
                        VClosed, VminMax = Vparams
                        if VClosed:
                            VVar = "%i" % idT
                        else:
                            surfParams, surfType, inSurf = gen_torus_annex_v_surface(
                                face, VminMax, opt.force_cylinder
                            )

                            if surfType == "Cone":
                                cone = GEOUNED_Surface(
                                    ("Cone", surfParams), UniverseBox, Face="Build"
                                )
                                id2, exist = Surfaces.addCone(cone)

                            elif surfType == "Cylinder":
                                cyl = GEOUNED_Surface(
                                    ("Cylinder", surfParams), UniverseBox, Face="Build"
                                )
                                id2, exist = Surfaces.addCylinder(cyl)

                            elif surfType == "Plane":
                                plane = GEOUNED_Surface(
                                    ("Plane", surfParams), UniverseBox, Face="Build"
                                )
                                id2, exist = Surfaces.addPlane(plane)
                                if exist:
                                    p = Surfaces.getSurface(id2)
                                    if is_opposite(
                                        plane.Surf.Axis, p.Surf.Axis, tol.pln_angle
                                    ):
                                        id2 = -id2

                            VVar = "%i %i" % (idT, -id2 if inSurf else id2)

                    var = VVar if UClosed else " ".join((VVar, UVar))
                    if var not in SurfPiece:
                        SurfPiece.append(var)
                        SurfObj.append(face)
                else:
                    if opt.verbose:
                        print(
                            "Only Torus with axis along X, Y , Z axis can be reproduced"
                        )
            else:
                id = get_id(face.Surface, Surfaces)
                if surfaceType == "<Cone object>":
                    cones.add(-id)

                surf = face
                if id == 0:
                    if opt.verbose:
                        print("Warning: ", surfaceType, " not found in surface list")
                    if surfaceType == "<Plane object>":
                        dim1 = face.ParameterRange[1] - face.ParameterRange[0]
                        dim2 = face.ParameterRange[3] - face.ParameterRange[2]
                        plane = GEOUNED_Surface(
                            (
                                "Plane",
                                (face.Surface.Position, face.Surface.Axis, dim1, dim2),
                            ),
                            UniverseBox,
                            Face="Build",
                        )
                        id, exist = Surfaces.addPlane(plane)
                        surf = plane.shape
                    elif surfaceType == "<Cylinder object>":
                        dimL = face.ParameterRange[3] - face.ParameterRange[2]
                        cylinder = GEOUNED_Surface(
                            (
                                "Cylinder",
                                (
                                    face.Surface.Center,
                                    face.Surface.Axis,
                                    face.Surface.Radius,
                                    dimL,
                                ),
                            ),
                            UniverseBox,
                            Face="Build",
                        )
                        id, exist = Surfaces.addCylinder(cylinder)
                        surf = cylinder.shape

                if orient == "Reversed":
                    var = id
                elif orient == "Forward":
                    var = -id

                if surfaceType == "<Plane object>":
                    s = Surfaces.getSurface(id)
                    if is_opposite(face.Surface.Axis, s.Surf.Axis, tol.pln_angle):
                        var = -var

                if str(var) in SurfPiece:
                    continue

                SurfPiece.append(str(var))
                SurfObj.append(surf)

        if extraPlaneReverse:
            for extra in extraPlaneReverse.values():
                if len(extra) == 1:
                    if extra[0] not in SurfPiece:
                        SurfPiece.append(extra[0])
                else:
                    SurfPiece.append("({})".format(":".join(extra)))

        SurfPieceBool = BoolSequence(" ".join(SurfPiece))
        # possible expresion for e
        #  i1
        #  i1 i2SurfPiece
        #  i1 i2 i3
        #  (i1:i2)
        #  i1 (i2:i3)

        #         semi = e.find(':')
        #         blk  = e.find(' ')
        #         print (e)
        #         if semi != -1 :
        #           orTerm = expOR(int(e[1:semi]),int(e[semi+1:-1]))
        #           SurfPieceBool.add(orTerm)
        #         elif blk != -1 :
        #           SurfPieceBool.add(int(e[1:blk]),int(e[blk+1:-1]))
        #         else :
        #          SurfPieceBool.add(int(e))

        if SurfPieceBool.elements:
            PieceDef.append(SurfPieceBool)
            PieceObj.append(SurfObj)
        else:
            delList.append(isol)

    for isol in reversed(delList):
        del metaObj.Solids[isol]
    metaObj.set_definition(PieceDef)
    metaObj.setFaces(PieceObj)
    return tuple(cones)


def get_surf_value(Definition, reverse=False):

    if Definition.level == 0:
        if reverse:
            surf = {-i for i in Definition.elements}
        else:
            surf = set(Definition.elements)
    else:
        surf = set()
        for e in Definition.elements:
            if e.operator == "AND":
                if reverse:
                    surf = {-i for i in e.elements}
                else:
                    surf = set(e.elements)
                break
    return surf


def append_comp(newCell, cellDef, cellCAD, metaComplementary):

    surfCell = get_surf_value(cellDef, True)
    if metaComplementary.Definition.operator == "AND":
        if not cellCAD.BoundBox.intersect(metaComplementary.CADSolid.BoundBox):
            return False
        Seq = metaComplementary.Definition
        surfComp = get_surf_value(Seq, False)
        if len(surfComp & surfCell) > 0:
            return False
        newCell.append(Seq.get_complementary())
        return True

    else:
        append = False
        for i, compPart in enumerate(metaComplementary.Solids):
            if not cellCAD.BoundBox.intersect(compPart.BoundBox):
                continue
            Seq = metaComplementary.Definition.elements[i]
            surfComp = get_surf_value(Seq, False)
            if len(surfComp & surfCell) > 0:
                continue
            append = True
            newCell.append(Seq.get_complementary())
        return append


def no_overlapping_cell(metaList, Surfaces):

    Surfs = {}
    for lst in Surfaces.values():
        for s in lst:
            Surfs[s.Index] = s

    newDefinitionList = []
    metaList[0].setCADSolid()

    for i, m in enumerate(metaList[1:]):
        m.setCADSolid()

        if m.Definition.operator == "AND":
            newDef = BoolSequence(operator="AND")
            newDef.append(m.Definition.copy())
            simplify = False
            for mm in metaList[: i + 1]:
                simp = append_comp(newDef, m.Definition, m.CADSolid, mm)
                if simp:
                    simplify = True
            simpTerm = [simplify]

        else:
            newDef = BoolSequence(operator="OR")
            simpTerm = []
            for j, partSolid in enumerate(m.Solids):
                subDef = BoolSequence(operator="AND")
                subDef.append(m.Definition.elements[j].copy())
                simplify = False
                for mm in metaList[: i + 1]:
                    simp = append_comp(subDef, m.Definition.elements[j], partSolid, mm)
                    if simp:
                        simplify = True
                simpTerm.append(simplify)
                newDef.append(subDef)
        newDefinitionList.append((newDef, simpTerm))

    for m, tDef in zip(metaList[1:], newDefinitionList):
        Def, simplify = tDef
        if True in simplify:
            print(f"reduce cell {m.__id__}")
            Box = UF.get_box(m)

            # evaluate only diagonal elements of the Constraint Table (fastest) and remove surface not
            # crossing in the solid boundBox
            CT = build_c_table_from_solids(
                Box, (tuple(Def.get_surfaces_numbers()), Surfs), option="diag"
            )

            newDef = remove_extra_surfaces(Def, CT)

            # evaluate full constraint Table with less surfaces involved
            CT = build_c_table_from_solids(
                Box, (tuple(newDef.get_surfaces_numbers()), Surfs), option="full"
            )

            if newDef.operator == "AND":
                newDef.simplify(CT)
                newDef.clean()
            else:
                for i, s in enumerate(simplify):
                    if not s:
                        continue
                    comp = newDef.elements[i]
                    comp.simplify(CT)
                    comp.clean()

            m.set_definition(newDef)
            m.Definition.join_operators()
            m.Definition.level_update()


def extra_plane_cyl_face(face, Box, Surfaces):
    wire = face.OuterWire
    planesId = []
    for e in wire.OrderedEdges:
        curve = str(e.Curve)
        if curve[0:6] == "Circle" or curve == "<Ellipse object>":
            dir = e.Curve.Axis
            center = e.Curve.Center
            if curve == "<Ellipse object>":
                dim1 = e.Curve.MinorRadius
                dim2 = e.Curve.MajorRadius
            else:
                dim1 = e.Curve.Radius
                dim2 = e.Curve.Radius
            plane = GEOUNED_Surface(
                ("Plane", (center, dir, dim1, dim2)), Box, Face="Build"
            )
            sign = sign_plane(face.CenterOfMass, plane)
            id, exist = Surfaces.addPlane(plane)
            if exist:
                pp = Surfaces.getSurface(id)
                if is_opposite(plane.Surf.Axis, pp.Surf.Axis, tol.pln_angle):
                    id = -id
            planesId.append(id)
    return planesId


def add_cone_plane(Definition, conesList, Surfaces, UniverseBox):
    x_axis = FreeCAD.Vector(1, 0, 0)
    y_axis = FreeCAD.Vector(0, 1, 0)
    z_axis = FreeCAD.Vector(0, 0, 1)

    for cid in conesList:
        cone = Surfaces.getSurface(abs(cid))
        if (
            is_parallel(cone.Surf.Axis, x_axis, tol.angle)
            or is_parallel(cone.Surf.Axis, y_axis, tol.angle)
            or is_parallel(cone.Surf.Axis, z_axis, tol.angle)
        ):
            continue

        plane = GEOUNED_Surface(
            ("Plane", (cone.Surf.Apex, cone.Surf.Axis, 1, 1)), UniverseBox, Face="Build"
        )
        pid, exist = Surfaces.addPlane(plane)

        if exist:
            p = Surfaces.getSurface(pid)
            if is_opposite(plane.Surf.Axis, p.Surf.Axis, tol.pln_angle):
                pid = -pid

        if cid > 0:
            insert_in_sequence(Definition, cid, -pid, "OR")
        else:
            insert_in_sequence(Definition, cid, pid, "AND")
