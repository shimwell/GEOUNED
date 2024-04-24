import math

import FreeCAD
import numpy as np
import Part

from .build_solid_cell import build_solid
from .remh import Cline
from .utils.booleanFunction import BoolSequence, outer_terms


class CADCell:
    def __init__(self, stringCell=None):

        if not stringCell:
            self.surfaces = {}
            self.surfaceList = []
            self.shape = None
            # self.likeCell = None
            self.definition = None
            self.name = 0
            # self.TRCL     = None  # cell transformacion "like-but" cells
            self.TRFL = None  # Universe transformation in fill Universe
            self.U = -1  # Cell Universe number
            self.FILL = 0  # Fill Universe number
            self.MAT = 0  # material number
            self.CurrentTR = None
            self.level = None
            self.__defTerms__ = None
            self.__operator__ = None
        else:
            self.surfaces = None
            self.shape = None
            self.name = stringCell.name
            self.TRFL = stringCell.TR  # Universe transformation in fill Universe
            self.U = stringCell.U  # Cell Universe number
            self.FILL = stringCell.FILL  # Fill Universe number
            self.MAT = stringCell.MAT  # material number
            self.CurrentTR = self.TRFL
            self.level = None

            self.__defTerms__ = None
            self.__operator__ = None
            self.__set_definition__(stringCell)

    def copy(self):
        cpCell = CADCell()
        cpCell.surfaceList = self.surfaceList[:]
        cpCell.surfaces = {}
        for name, s in self.surfaces.items():
            cpCell.surfaces[name] = s.copy()

        if type(self.definition) is Cline:
            cpCell.definition = Cline(self.definition.str)

        elif type(self.definition) is BoolSequence:
            cpCell.definition = self.definition.copy()

        cpCell.name = self.name
        cpCell.TRFL = self.TRFL
        cpCell.U = self.U
        cpCell.FILL = self.FILL
        cpCell.MAT = self.MAT
        cpCell.level = self.level

        if self.CurrentTR is not None:
            cpCell.CurrentTR = self.CurrentTR.submatrix(4)

        if self.shape is not None:
            cpCell.shape = self.shape.copy()

        return cpCell

    def get_sub_cell(self, seq):

        subCell = self.copy()
        subCell.definition = seq.copy()
        subCell.shape = None
        subCell.surfaceList = subCell.definition.get_surfaces_numbers()
        for s in tuple(subCell.surfaces.keys()):
            if s not in subCell.surfaceList:
                del subCell.surfaces[s]

        return subCell

    # not used
    #
    #    def split(self,nparts=2):
    #
    #        if nparts == 1:
    #            return (self,None)
    #        terms,operador = self.get_outer_terms()
    #        nelemts = int(len(terms)/nparts)
    #        subDefList = []
    #
    #        if operador == 'AND':
    #          for i in range(nparts-1):
    #             newdef = ') ('.join(terms[i*nelemts:(i+1)*nelemts])
    #             newdef = '({})'.format(newdef)
    #             subDefList.append(newdef)
    #          newdef = ') ('.join(terms[(nparts-1)*nelemts:])
    #          newdef = '({})'.format(newdef)
    #          subDefList.append(newdef)
    #
    #        else:
    #          for i in range(nparts-1):
    #             newdef = '):('.join(terms[i*nelemts:(i+1)*nelemts])
    #             newdef = '({})'.format(newdef)
    #             subDefList.append(newdef)
    #          newdef = '):('.join(terms[(nparts-1)*nelemts:])
    #          newdef = '({})'.format(newdef)
    #          subDefList.append(newdef)
    #
    #
    #        subCellList=[]
    #        for df in subDefList:
    #           subCell = self.copy()
    #           subCell.definition= Cline(df)
    #           subCell.shape = None
    #           subCell.surfaceList  = subCell.definition.get_surfaces_numbers()
    #           for s in tuple(subCell.surfaces.keys()) :
    #               if s not in subCell.surfaceList: del(subCell.surfaces[s])
    #
    #           subCellList.append(subCell)
    #
    #        return subCellList,operador

    def get_outer_terms(self):
        if not self.__defTerms__:
            self.__defTerms__, self.__operator__ = outer_terms(self.definition.str)
        return self.__defTerms__, self.__operator__

    def make_box(self, boundBox):
        box_origin = FreeCAD.Vector(boundBox.XMin, boundBox.YMin, boundBox.ZMin)
        return Part.makeBox(
            boundBox.XLength, boundBox.YLength, boundBox.ZLength, box_origin
        )

    def build_shape(
        self, boundBox, force=False, surfTR=None, simplify=False, fuse=False
    ):

        if self.shape is not None and not force:
            return
        if surfTR:
            self.transform_surfaces(surfTR)

        cutShape = build_solid(self, boundBox, simplify=simplify)

        if fuse or True:
            self.shape = fuse_solid(cutShape)
        else:
            self.shape = Part.makeCompound(cutShape)

    def build_surface_shape(self, boundBox):
        for s in self.surfaces.values():
            s.build_shape(boundBox)

    def transform_solid(self, matrix, reverse=False):
        if not self.shape:
            return
        if reverse:
            self.shape = self.shape.transformGeometry(matrix.inverse())
        else:
            self.shape = self.shape.transformGeometry(matrix)

    def transform_surfaces(self, matrix):
        for s in self.surfaces.values():
            s.transform(matrix)

    def set_surfaces(self, Surfaces):
        if self.surfaces is not None:
            return
        self.surfaces = {}
        for s in self.surfaceList:
            self.surfaces[s] = Surfaces[s]

    def clean_undefined(self):
        undefined = []
        for s in self.definition.get_surfaces_numbers():
            if self.surfaces[s].params is None:
                undefined.append(s)
        if undefined:
            self.definition.remove_surface(undefined)

        for s in undefined:
            del self.surfaces[s]

    def __set_definition__(self, stringCell):

        self.definition = stringCell.geom
        self.definition.remove_comments(full=True)
        self.definition.remove_cr()
        self.definition.remove_multispace()
        self.definition.remove_redundant()
        self.surfaceList = self.definition.get_surfaces_numbers()


class Plane:
    def __init__(self, Id, params, tr=None):
        self.type = "plane"
        self.id = Id
        self.shape = None
        self.params = params
        if tr:
            self.transform(tr)

    def __str__(self):
        return "plane : {}\nparameters : {}".format(self.id, self.params)

    def copy(self):
        return Plane(self.id, self.params)

    def transform(self, matrix):
        v, d = self.params
        p = d * v  # vector p is d*plane normal
        v = matrix.submatrix(3).multVec(v)
        v.normalize()
        d = matrix.multVec(p) * v
        self.params = (v, d)

    def build_shape(self, boundBox):
        normal, p0 = self.params
        Box = FreeCAD.BoundBox(boundBox)
        Box.enlarge(10)

        pointEdge = []
        for i in range(12):
            edge = Box.getEdge(i)
            p1 = normal.dot(edge[0])
            p2 = normal.dot(edge[1])
            d0 = p0 - p1
            d1 = p2 - p1
            if d1 != 0:
                a = d0 / d1
                if a >= 0 and a <= 1:
                    pointEdge.append(edge[0] + a * (edge[1] - edge[0]))

        if len(pointEdge) == 0:
            return
        s = FreeCAD.Vector((0, 0, 0))
        for v in pointEdge:
            s = s + v
        s = s / len(pointEdge)

        vtxvec = []
        for v in pointEdge:
            vtxvec.append(v - s)

        X0 = vtxvec[0]
        Y0 = normal.cross(X0)

        orden = []
        for i, v in enumerate(vtxvec):
            phi = np.arctan2(v.dot(Y0), v.dot(X0))
            orden.append((phi, i))
        orden.sort()

        self.shape = Part.Face(Part.makePolygon([pointEdge[p[1]] for p in orden], True))


class Sphere:
    def __init__(self, Id, params, tr=None):
        self.type = "sphere"
        self.id = Id
        self.shape = None
        self.params = params
        if tr:
            self.transform(tr)

    def copy(self):
        return Sphere(self.id, self.params)

    def transform(self, matrix):
        p, R = self.params
        p = matrix.multVec(p)
        self.params = (p, R)

    def build_shape(self, boundBox):
        origin, R = self.params
        self.shape = Part.makeSphere(R, origin)


class Cylinder:
    def __init__(self, Id, params, tr=None, truncated=False):
        self.type = "cylinder"
        self.id = Id
        self.shape = None
        self.params = params
        self.truncated = truncated
        if tr:
            self.transform(tr)

    def copy(self):
        return Cylinder(self.id, self.params, truncated=self.truncated)

    def transform(self, matrix):
        p, v, R = self.params
        v = matrix.submatrix(3).multVec(v)
        p = matrix.multVec(p)
        self.params = (p, v, R)

    def build_shape(self, boundBox):

        p, vec, r = self.params

        if not self.truncated:
            dmin = vec.dot(boundBox.getPoint(0) - p)
            dmax = dmin
            for i in range(1, 8):
                d = vec.dot(boundBox.getPoint(i) - p)
                dmin = min(d, dmin)
                dmax = max(d, dmax)

            height = dmax - dmin
            dmin -= 0.1 * height
            dmax += 0.1 * height
            height = dmax - dmin

            point = p + dmin * vec
            self.shape = Part.makeCylinder(r, height, point, vec, 360)
            # self.shape = makeCylinder2( r,height,point,vec)
        else:
            self.shape = Part.makeCylinder(r, vec.Length, p, vec, 360)
            # self.shape = Part.makeCylinder2( r,vec.Length,p,vec)

        return


class Cone:
    def __init__(self, Id, params, tr=None, truncated=False):
        self.type = "cone"
        self.id = Id
        self.shape = None
        self.params = params
        self.truncated = truncated
        if tr:
            self.transform(tr)

    def copy(self):
        return Cone(self.id, self.params, truncated=self.truncated)

    def transform(self, matrix):
        if not self.truncated:
            p, v, t, dbl = self.params
            v = matrix.submatrix(3).multVec(v)
            p = matrix.multVec(p)
            self.params = (p, v, t, dbl)
        else:
            p, v, r1, r2 = self.params
            v = matrix.submatrix(3).multVec(v)
            p = matrix.multVec(p)
            self.params = (p, v, r1, r2)

    def build_shape(self, boundBox):
        if not self.truncated:
            apex, axis, t, dblsht = self.params

            dmin = axis.dot(boundBox.getPoint(0) - apex)
            dmax = dmin
            for i in range(1, 8):
                d = axis.dot(boundBox.getPoint(i) - apex)
                dmin = min(d, dmin)
                dmax = max(d, dmax)

            length = max(abs(dmin), abs(dmax))
            R = length * t
            OneSheetCone = Part.makeCone(0, R, length, apex, axis, 360)
            if not dblsht:
                self.shape = OneSheetCone
            else:
                OtherSheet = Part.makeCone(0, R, length, apex, -axis, 360)
                DoubleSheetCone = OneSheetCone.fuse([OtherSheet])
                DoubleSheetCone.removeSplitter()
                self.shape = DoubleSheetCone
        else:
            center, axis, r1, r2 = self.params
            self.shape = Part.makeCone(r1, r2, axis.Length, center, axis, 360)


class EllipticCone:
    def __init__(self, Id, params, tr=None):
        self.type = "cone_elliptic"
        self.id = Id
        self.shape = None
        self.params = params
        if tr:
            self.transform(tr)

    def copy(self):
        return EllipticCone(self.id, self.params)

    def transform(self, matrix):
        p, v, ra, radii, raxes, dbl = self.params
        v = matrix.submatrix(3).multVec(v)
        raxes[0] = matrix.submatrix(3).multVec(raxes[0])
        raxes[1] = matrix.submatrix(3).multVec(raxes[1])
        p = matrix.multVec(p)
        self.params = (p, v, ra, radii, raxes, dbl)

    def build_shape(self, boundBox):
        apex, axis, ra, radii, raxes, dblsht = self.params

        dmin = axis.dot(boundBox.getPoint(0) - apex)
        dmax = dmin
        for i in range(1, 8):
            d = axis.dot(boundBox.getPoint(i) - apex)
            dmin = min(d, dmin)
            dmax = max(d, dmax)

        length = max(abs(dmin), abs(dmax))
        OneSheetCone = make_elliptic_cone(apex, axis, ra, radii, raxes, length)
        if not dblsht:
            self.shape = OneSheetCone
        else:
            OtherSheet = make_elliptic_cone(apex, -axis, ra, radii, raxes, length)
            DoubleSheetCone = OneSheetCone.fuse([OtherSheet])
            DoubleSheetCone.removeSplitter()
            self.shape = DoubleSheetCone


class Hyperboloid:
    def __init__(self, Id, params, tr=None):
        self.type = "hyperboloid"
        self.id = Id
        self.shape = None
        self.params = params
        if tr:
            self.transform(tr)

    def copy(self):
        return Hyperboloid(self.id, self.params)

    def transform(self, matrix):
        p, v, radii, raxes, onesht = self.params
        v = matrix.submatrix(3).multVec(v)
        raxes[0] = matrix.submatrix(3).multVec(raxes[0])
        raxes[1] = matrix.submatrix(3).multVec(raxes[1])
        p = matrix.multVec(p)
        self.params = (p, v, radii, raxes, onesht)

    def build_shape(self, boundBox):
        center, axis, radii, rAxes, onesht = self.params

        dmin = axis.dot(boundBox.getPoint(0) - center)
        dmax = dmin
        for i in range(1, 8):
            d = axis.dot(boundBox.getPoint(i) - center)
            dmin = min(d, dmin)
            dmax = max(d, dmax)

        length = max(abs(dmin), abs(dmax))
        self.shape = make_hyperboloid(center, radii, rAxes, axis, onesht, length)


class Ellipsoid:
    def __init__(self, Id, params, tr=None):
        self.type = "ellipsoid"
        self.id = Id
        self.shape = None
        self.params = params
        if tr:
            self.transform(tr)

    def copy(self):
        return Ellipsoid(self.id, self.params)

    def transform(self, matrix):
        p, v, radii, raxes = self.params
        v = matrix.submatrix(3).multVec(v)
        raxes[0] = matrix.submatrix(3).multVec(raxes[0])
        raxes[1] = matrix.submatrix(3).multVec(raxes[1])
        p = matrix.multVec(p)
        self.params = (p, v, radii, raxes)

    def build_shape(self, boundBox):
        center, axis, radii, rAxes = self.params
        self.shape = make_ellipsoid(center, radii, rAxes, axis)


class EllipticCylinder:
    def __init__(self, Id, params, tr=None, truncated=False):
        self.type = "cylinder_elliptic"
        self.id = Id
        self.shape = None
        self.params = params
        self.truncated = truncated
        if tr:
            self.transform(tr)

    def copy(self):
        return EllipticCylinder(self.id, self.params, truncated=self.truncated)

    def transform(self, matrix):
        p, v, radii, raxes = self.params
        v = matrix.submatrix(3).multVec(v)
        raxes[0] = matrix.submatrix(3).multVec(raxes[0])
        raxes[1] = matrix.submatrix(3).multVec(raxes[1])
        p = matrix.multVec(p)
        self.params = (p, v, radii, raxes)

    def build_shape(self, boundBox):
        center, axis, radii, rAxes = self.params
        if not self.truncated:
            dmin = axis.dot(boundBox.getPoint(0) - center)
            dmax = dmin
            for i in range(1, 8):
                d = axis.dot(boundBox.getPoint(i) - center)
                dmin = min(d, dmin)
                dmax = max(d, dmax)

            height = dmax - dmin
            dmin -= 0.1 * height
            dmax += 0.1 * height
            height = dmax - dmin
            point = center + dmin * axis

            self.shape = make_elliptic_cylinder(point, radii, rAxes, axis, height)
        else:
            height = axis.Length
            self.shape = make_elliptic_cylinder(
                center, radii, rAxes, axis / height, height
            )


class HyperbolicCylinder:
    def __init__(self, Id, params, tr=None):
        self.type = "cylinder_hyperbolic"
        self.id = Id
        self.shape = None
        self.params = params
        if tr:
            self.transform(tr)

    def copy(self):
        return HyperbolicCylinder(self.id, self.params)

    def transform(self, matrix):
        p, v, radii, raxes = self.params
        v = matrix.submatrix(3).multVec(v)
        raxes[0] = matrix.submatrix(3).multVec(raxes[0])
        raxes[1] = matrix.submatrix(3).multVec(raxes[1])
        p = matrix.multVec(p)
        self.params = (p, v, radii, raxes)

    def build_shape(self, boundBox):
        center, axis, radii, rAxes = self.params

        dmin = axis.dot(boundBox.getPoint(0) - center)
        dmax = dmin
        for i in range(1, 8):
            d = axis.dot(boundBox.getPoint(i) - center)
            dmin = min(d, dmin)
            dmax = max(d, dmax)

        height = dmax - dmin
        dmin -= 0.1 * height
        dmax += 0.1 * height
        height = dmax - dmin
        point = center + dmin * axis

        self.shape = make_hyperbolic_cylinder(point, radii, rAxes, axis, height)


class Paraboloid:
    def __init__(self, Id, params, tr=None):
        self.type = "paraboloid"
        self.id = Id
        self.shape = None
        self.params = params
        if tr:
            self.transform(tr)

    def copy(self):
        return Paraboloid(self.id, self.params)

    def transform(self, matrix):
        p, v, focal = self.params
        v = matrix.submatrix(3).multVec(v)
        p = matrix.multVec(p)
        self.params = (p, v, focal)

    def build_shape(self, boundBox):
        center, axis, focal = self.params

        dmin = axis.dot(boundBox.getPoint(0) - center)
        dmax = dmin
        for i in range(1, 8):
            d = axis.dot(boundBox.getPoint(i) - center)
            dmin = min(d, dmin)
            dmax = max(d, dmax)

        length = max(abs(dmin), abs(dmax))
        self.shape = make_paraboloid(center, axis, focal, length)


class Torus:
    def __init__(self, Id, params, tr=None):
        self.type = "torus"
        self.id = Id
        self.shape = None
        self.params = params
        if tr:
            self.transform(tr)

    def copy(self):
        return Torus(self.id, self.params)

    def transform(self, matrix):
        p, v, Ra, Rb, Rc = self.params
        v = matrix.submatrix(3).multVec(v)
        p = matrix.multVec(p)
        self.params = (p, v, Ra, Rb, Rc)

    def build_shape(self, boundBox):
        center, axis, Ra, Rb, Rc = (
            self.params
        )  # Ra distance from torus axis; R radius of toroidal-cylinder
        if (abs(Rb - Rc) < 1e-5) and Ra > 0:
            self.shape = Part.makeTorus(Ra, Rb, center, axis)  # FreeCAD circular Torus
        else:
            self.shape = make_elliptic_torus(
                Ra, Rb, Rc, center, axis
            )  # Home made elliptic Torus


class Box:
    def __init__(self, Id, params, tr=None):
        self.type = "box"
        self.id = Id
        self.shape = None
        self.params = params
        if tr:
            self.transform(tr)

    def copy(self):
        return Box(self.id, self.params)

    def transform(self, matrix):
        p, v1, v2, v3 = self.params
        p = matrix.multVec(p)
        v1 = matrix.multVec(v1)
        v2 = matrix.multVec(v2)
        v3 = matrix.multVec(v3)
        self.params = (p, v1, v2, v3)

    def build_shape(self, boundBox):
        p, v1, v2, v3 = self.params
        a1 = FreeCAD.Vector(v1)
        a2 = FreeCAD.Vector(v2)
        a3 = FreeCAD.Vector(v3)
        a1.normalize()
        a2.normalize()
        a3.normalize()

        m = FreeCAD.Matrix(
            a1.x,
            a2.x,
            a3.x,
            p.x,
            a1.y,
            a2.y,
            a3.y,
            p.y,
            a1.z,
            a2.z,
            a3.z,
            p.z,
            0,
            0,
            0,
            1,
        )
        box = Part.makeBox(v1.Length, v2.Length, v3.Length)
        self.shape = box.transformGeometry(m)


class Undefined:
    def __init__(self, Id):
        self.type = "Undefined"
        self.id = Id
        self.shape = None
        self.params = None

    def copy(self):
        return Undefined(self.id)

    def build_shape(self, boundBox):
        return

    def transform(self, matrix):
        return


def fuse_solid(parts):
    if (len(parts)) <= 1:
        if parts:
            solid = parts[0]
        else:
            return None
    else:
        try:
            fused = parts[0].fuse(parts[1:])
        except:
            fused = None

        if fused is not None:
            try:
                refinedfused = fused.removeSplitter()
            except:
                refinedfused = fused

            if refinedfused.isValid():
                solid = refinedfused
            else:
                if fused.isValid():
                    solid = fused
                else:
                    solid = Part.makeCompound(parts)
        else:
            solid = Part.makeCompound(parts)

    if solid.Volume < 0:
        solid.reverse()
    return solid


def make_hyperboloid(center, radii, rAxes, axis, onesht, length):

    S1 = center + rAxes[1] * radii[1]  # major axis
    S2 = center + rAxes[0] * radii[0]  # minor axis
    hyperbola = Part.Hyperbola(S1, S2, center)

    Y = length
    X = radii[1] * math.sqrt((Y / radii[0]) ** 2 + 1)
    point = (
        center + X * radii[1] + Y * radii[0]
    )  # point in taken as length is always counted on minor axis
    parameter = abs(hyperbola.parameter(point))

    if onesht:
        shape = hyperbola.toBSpline(-parameter, parameter).toShape(
            -parameter, parameter
        )
        hyperFace = shape.revolve(center, axis, 360)

        StartPoint = hyperFace.Surface.BasisCurve.StartPoint - center
        EndPoint = hyperFace.Surface.BasisCurve.EndPoint - center

        rad1 = StartPoint.dot(hyperbola.XAxis)
        hgt1 = StartPoint.dot(hyperbola.YAxis)
        cc1 = center + hyperbola.YAxis * hgt1
        circle1 = Part.Circle(cc1, -hyperbola.YAxis, rad1).toShape()
        cFace1 = Part.makeFace(circle1, "Part::FaceMakerSimple")

        rad2 = EndPoint.dot(hyperbola.XAxis)
        hgt2 = EndPoint.dot(hyperbola.YAxis)
        cc2 = center + hyperbola.YAxis * hgt2
        circle2 = Part.Circle(cc2, hyperbola.YAxis, rad2).toShape()
        cFace2 = Part.makeFace(circle2, "Part::FaceMakerSimple")

        shell = Part.makeShell((cFace1, hyperFace, cFace2))
        hyperboloid = Part.makeSolid(shell)
    else:
        shape = hyperbola.toBSpline(0, parameter).toShape(0, parameter)
        hyperFace = shape.revolve(center, axis, 360)

        rad = EndPoint.dot(hyperbola.YAxis)
        hgt = EndPoint.dot(hyperbola.XAxis)
        cc = center + hyperbola.XAxis * hgt
        circle = Part.Circle(cc, -hyperbola.XAxis, rad).toShape()
        cFace = Part.makeFace(circle, "Part::FaceMakerSimple")

        shell = Part.makeShell((cFace, hyperFace))
        hyper1 = Part.makeSolid(shell)
        hyper2 = hyper1.rotated(center, hyperbola.YAxis, 180)
        hyperboloid = Part.makeCompound((hyper1, hyper2))

    return hyperboloid


def make_hyperbolic_cylinder(center, radii, rAxes, axis, length):

    S11 = center + rAxes[1] * radii[1]  # major axis
    S12 = center + rAxes[0] * radii[0]  # minor axis
    S21 = center - rAxes[1] * radii[1]  # major axis
    S22 = center - rAxes[0] * radii[0]  # minor axis

    hyperbola1 = Part.Hyperbola(S11, S12, center)
    hyperbola2 = Part.Hyperbola(S21, S22, center)
    d = axis * length

    Y = length
    X = radii[1] * math.sqrt((Y / radii[0]) ** 2 + 1)
    point = (
        center + X * rAxes[1] + Y * rAxes[0]
    )  # point in taken as length is always counted on minor axis
    parameter = abs(hyperbola1.parameter(point))

    shape1 = hyperbola1.toBSpline(-parameter, parameter).toShape(-parameter, parameter)
    shape2 = hyperbola2.toBSpline(-parameter, parameter).toShape(-parameter, parameter)
    surf1 = shape1.extrude(d)
    surf2 = shape2.extrude(d)

    return Part.makeCompound((surf1, surf2))


def make_elliptic_cylinder(center, radii, rAxes, axis, length):

    S1 = center + rAxes[1] * radii[1]  # major axis
    S2 = center + rAxes[0] * radii[0]  # minor axis
    d = axis * length

    ellipse = Part.Ellipse(S1, S2, center)
    ellipse2 = Part.Ellipse(S1 + d, S2 + d, center + d)

    shape = ellipse.toBSpline().toShape()
    shape2 = ellipse2.toBSpline().toShape()
    shell = Part.makeLoft((shape, shape2), True)
    return Part.makeSolid(shell)


def make_ellipsoid(center, radii, rAxes, axis):

    S1 = center + rAxes[1] * radii[1]  # major axis
    S2 = center + rAxes[0] * radii[0]  # minor axis

    ellipse = Part.Ellipse(S1, S2, center)

    if axis.add(-rAxes[0]).Length < 1e-5:
        shape = ellipse.toBSpline().toShape()
        shell = shape.revolve(center, axis, 180)
    else:
        shape = ellipse.toBSpline(0, math.pi).toShape(0, math.pi)
        shell = shape.revolve(center, axis, 360)
    return Part.makeSolid(shell)


def make_elliptic_torus(R, RZ, RX, center, ZAxis):

    rMaj = RZ
    rMin = RX
    XAxis = orto_vect(ZAxis)

    majorAxis = ZAxis
    minorAxis = XAxis
    if rMaj < rMin:
        rMaj, rMin = rMin, rMaj
        majorAxis, minorAxis = minorAxis, majorAxis

    eCenter = center + R * XAxis
    S1 = eCenter + majorAxis * rMaj  # major axis
    S2 = eCenter + minorAxis * rMin  # minor axis

    ellipse = Part.Ellipse(S1, S2, eCenter)
    if abs(R) < RX:  # degenerated Torus
        pz = RZ * math.sqrt(1 - (R / RX) ** 2)
        pz1 = center - pz * ZAxis
        pz2 = center + pz * ZAxis

        p1 = ellipse.parameter(pz1)
        p2 = ellipse.parameter(pz2)
        if p2 < p1:
            p2 += 2 * math.pi
        shape = ellipse.toBSpline(p1, p2).toShape(
            p1, p2
        )  # revolution around Major axis
        rev = shape.revolve(center, ZAxis, 360)
    else:
        shape = ellipse.toBSpline().toShape()  # revolution around Minor axis
        rev = shape.revolve(center, ZAxis, 360)
    shell = Part.makeShell((rev,))
    return Part.makeSolid(shell)


def make_paraboloid(center, axis, focal, length):

    parabola = Part.Parabola()
    parabola.Center = center
    parabola.Axis = orto_vect(axis)
    parabola.XAxis = axis
    parabola.Focal = focal

    R = math.sqrt(4 * focal * length)
    point = center + length * parabola.XAxis + R * parabola.YAxis
    parameter = abs(parabola.parameter(point))

    shape = parabola.toBSpline(0, parameter).toShape(0, parameter)
    paraFace = shape.revolve(center, axis, 360)

    cc = center + length * parabola.XAxis
    circle = Part.Circle(cc, -axis, R).toShape()
    cFace = Part.makeFace(circle, "Part::FaceMakerSimple")

    shell = Part.makeShell((cFace, paraFace))
    return Part.makeSolid(shell)


def make_elliptic_cone(apex, axis, Ra, radii, rAxes, length):

    S1 = apex + rAxes[1] * radii[1] / Ra * length  # major axis
    S2 = apex + rAxes[0] * radii[0] / Ra * length  # minor axis
    d = axis * length

    point = Part.Point(apex).toShape()
    ellipse = Part.Ellipse(S1 + d, S2 + d, apex + d)

    shape = ellipse.toBSpline().toShape()
    shell = Part.makeLoft((point, shape), True)
    return Part.makeSolid(shell)


def orto_vect(v):
    vmax = 0
    vOrto = None
    if abs(v.x) > vmax:
        vOrto = (0, 1, 0)
        vmax = abs(v.x)
    if abs(v.y) > vmax:
        vOrto = (0, 0, 1)
        vmax = abs(v.y)
    if abs(v.z) > vmax:
        vOrto = (1, 0, 0)
        vmax = abs(v.z)

    if vOrto is None:
        return None

    vOrto = v.cross(FreeCAD.Vector(vOrto))
    vOrto.normalize()
    return vOrto
