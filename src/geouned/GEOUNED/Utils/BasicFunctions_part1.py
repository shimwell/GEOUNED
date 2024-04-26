#
# Set of useful functions used in different parts of the code
#
import math

import FreeCAD


def isSameValue(v1, v2, tolerance=1e-6):
    return abs(v1 - v2) < tolerance


def isOposite(Vector1, Vector2, tolerance=1e-6):
    return abs(Vector1.getAngle(-Vector2)) < tolerance


def isParallel(Vector1, Vector2, tolerance=1e-6):
    angle = abs(Vector1.getAngle(Vector2))
    return angle < tolerance or isSameValue(angle, math.pi, tolerance)


def isInLine(point, Dir, Pnt_line, tolerance=1e-6):
    r12 = point - Pnt_line
    return isParallel(Dir, r12) or (r12.Length < tolerance)


def isInPoints(point, points, tolerance=1e-5):
    if len(points) > 0:
        for P in points:
            if point.isEqual(P, tolerance):
                return True
    return False


def isInEdge(edge1, edge2, tolerance=1e-8):
    ver_1 = edge1.Vertexes
    ver_2 = edge2.Vertexes
    con_1 = ver_1[0].Point.isEqual(ver_2[0].Point, tolerance) or ver_1[0].Point.isEqual(
        ver_2[1].Point, tolerance
    )
    con_2 = ver_1[1].Point.isEqual(ver_2[0].Point, tolerance) or ver_1[1].Point.isEqual(
        ver_2[1].Point, tolerance
    )
    return con_1 and con_2


def isInPlane(point, plane, dtolerance=1e-7):
    return abs(point.distanceToPlane(plane.Surf.Position, plane.Surf.Axis)) < dtolerance


def isInTolerance(val, tol, fuzzyLow, fuzzyHigh):
    if abs(val) < fuzzyLow:
        return True, False  # 1) isintolerance 2) fuzzy
    elif abs(val) < tol:
        return True, True
    elif abs(val) > fuzzyHigh:
        return False, False
    else:
        return False, True


def signPlane(point, plane):
    value = plane.Surf.Axis.dot(point - plane.Surf.Position)
    if value >= 0.0:
        sign = 1
    else:
        sign = -1
    return sign


def pointsToCoeffs(points):
    p1, p2, p3 = points[0:3]
    scf = (p1.x, p1.y, p1.z, p2.x, p2.y, p2.z, p3.x, p3.y, p3.z)

    # mcnp implementation to convert 3 point plane to
    # plane parameters

    tpp = [0] * 4
    for i in range(1, 4):
        j = i % 3 + 1
        k = 6 - i - j
        k -= 1
        j -= 1
        tpp[i - 1] = (
            scf[j] * (scf[k + 3] - scf[k + 6])
            + scf[j + 3] * (scf[k + 6] - scf[k])
            + scf[j + 6] * (scf[k] - scf[k + 3])
        )
        tpp[3] += scf[i - 1] * (scf[j + 3] * scf[k + 6] - scf[j + 6] * scf[k + 3])

    xm = 0
    coeff = [0] * 4
    for i in range(1, 5):
        if xm == 0 and tpp[4 - i] != 0:
            xm = 1 / tpp[4 - i]
        coeff[4 - i] = tpp[4 - i] * xm

    axis = FreeCAD.Vector(coeff[0:3])
    distance = coeff[3] / axis.Length
    axis.normalize()

    return axis, distance


class Plane3PtsParams:
    def __init__(self, params, real=True):
        self.position = params[0]
        self.axis = params[1]
        self.dim_l1 = params[2]
        self.dim_l2 = params[3]
        self.points = params[4]
        self.real = real
        self.point_def = True

    def __str__(self):
        #      outstr = '''Plane :
        #    point 1  : {P1[0]}  {P1[1]}  {P1[2]}
        #    point 2  : {P2[0]}  {P2[1]}  {P2[2]}
        #    point 3  : {P3[0]}  {P3[1]}  {P3[2]} '''.format(P1=self.points[0], P2=self.points[1], P3=self.points[2])
        pos = self.axis.dot(self.position)
        outstr = """Plane :
    axis     : {}  {}  {} 
    position : {}  """.format(
            self.axis.x, self.axis.y, self.axis.z, pos
        )
        return outstr


class PlaneParams:
    def __init__(self, params, real=True):
        self.position = params[0]
        self.axis = params[1]
        self.dim_l1 = params[2]
        self.dim_l2 = params[3]
        self.real = real
        self.point_def = False

    def __str__(self):
        pos = self.axis.dot(self.position)
        outstr = """Plane :
    axis     : {}  {}  {} 
    position : {}  """.format(
            self.axis.x, self.axis.y, self.axis.z, pos
        )
        return outstr


class CylinderParams:
    def __init__(self, params, real=True):
        self.center = params[0]
        self.axis = params[1]
        self.Radius = params[2]
        self.dim_l = params[3]
        self.real = real

    def __str__(self):
        outstr = """Cylinder :
    axis     : {}  {}  {} 
    center   : {}  {}  {}
    Radius   : {}  """.format(
            self.axis.x,
            self.axis.y,
            self.axis.z,
            self.center.x,
            self.center.y,
            self.center.z,
            self.Radius,
        )
        return outstr


class ConeParams:
    def __init__(self, params, real=True):
        self.apex = params[0]
        self.axis = params[1]
        self.semi_angle = params[2]
        self.dim_l = params[3]
        self.dimR = params[4]
        self.real = real

    def __str__(self):
        outstr = """Cone :
    axis     : {}  {}  {} 
    center   : {}  {}  {}
    semi_angle: {}  """.format(
            self.axis.x,
            self.axis.y,
            self.axis.z,
            self.apex.x,
            self.apex.y,
            self.apex.z,
            self.semi_angle,
        )
        return outstr


class SphereParams:
    def __init__(self, params):
        self.center = params[0]
        self.radius = params[1]

    def __str__(self):
        outstr = """Sphere :
    center   : {}  {}  {}
    radius   : {}  """.format(
            self.center.x, self.center.y, self.center.z, self.radius
        )
        return outstr


class TorusParams:
    def __init__(self, params):
        self.center = params[0]
        self.Axis = params[1]
        self.major_radius = params[2]
        self.minor_radius = params[3]

    def __str__(self):
        outstr = """Torus :
    Axis     : {}  {}  {} 
    center   : {}  {}  {}
    major_radius: {}
    minor_radius: {} """.format(
            self.Axis.x,
            self.Axis.y,
            self.Axis.z,
            self.center.x,
            self.center.y,
            self.center.z,
            self.major_radius,
            self.minor_radius,
        )
        return outstr
