from numbers import Real


class Options:
    """A class for containing conversion options

    Args:
        force_cylinder (bool, optional): Use cylinder (instead of cones) as
            ancillary surface where unclosed torus surfaces are involved in
            the solid definition. Defaults to False.
        new_split_plane (bool, optional): New method to consider plane as
            cutting surface during the decomposition process. Former method
            split first planes perpendicular to X,Y,Z axis and then the
            other planes involved in the solid definition. New method group
            all parallel planes independently whether their normal are
            along X,Y,Z axes, and start the decomposition process cutting
            first with the group having the highest number of parallel
            planes. Defaults to True.
        del_last_number (bool, optional): Deleting the last word in the
            comment if it is a number. Defaults to False.
        enlarge_box (Real, optional): Enlarge box boundary when evaluating
            the constraint table during the simplification of the void cell
            definition. (unit is millimeter). Defaults to 2.
        n_plane_reverse (int, optional): Threshold value to determine whether
            cut with parallel planes should be carried out first. Defaults
            to 0.
        split_tolerance (Real, optional): Fuzzy tolerance value used in the
            FreeCAD function “BOPTools.SplitAPI.slice”. This function is
            used during the solid decomposition process. Defaults to 0.
        scale_up (bool, optional): Scale up Fuzzy tolerance once get below
            1e-12. Defaults to True.
        quadric_py (bool, optional): In openMC python script format, the
            cones or cylinders no aligned with the X,Y, or Z axis can be
            defined using the openmc.Cone or open.Cylinder methods but can
            also be defined with their quadric parameter. If “quadricPY” is
            11 True then all cones and cylinders will be defined in the
            openMC python script format under their quadric form. Defaults
            to False.
        facets (bool, optional): use alternative conversion module when
            geometry is defined by cells compound by only triangular plane
            faces. Defaults to False.
        prnt_3p_plane (bool, optional): print 3 point plane definition in
            output as 3 points coordinates. Defaults to False.
        forceNoOverlap (bool, optional): force no overlaping cell
            definition. Adjacent cell definition are rested from current
            cell definition. Defaults to False.
    """

    def __init__(
        self,
        force_cylinder: bool = False,
        new_split_plane: bool = True,
        del_last_number: bool = False,
        enlarge_box: Real = 2.0,
        n_plane_reverse: int = 0,
        split_tolerance: Real = 0.0,
        scale_up: bool = True,
        quadric_py: bool = False,
        facets: bool = False,
        prnt_3p_plane: bool = False,
        forceNoOverlap: bool = False,
    ):

        self.force_cylinder = force_cylinder
        self.new_split_plane = new_split_plane
        self.del_last_number = del_last_number
        self.enlarge_box = enlarge_box
        self.n_plane_reverse = n_plane_reverse
        self.split_tolerance = split_tolerance
        self.scale_up = scale_up
        self.quadric_py = quadric_py
        self.facets = facets
        self.prnt_3p_plane = prnt_3p_plane
        self.forceNoOverlap = forceNoOverlap

    @property
    def force_cylinder(self):
        return self._forceCylinder

    @force_cylinder.setter
    def force_cylinder(self, value: bool):
        if not isinstance(value, bool):
            raise TypeError(f"geouned.Options.force_cylinder should be a bool, not a {type(value)}")
        self._forceCylinder = value

    @property
    def new_split_plane(self):
        return self._newSplitPlane

    @new_split_plane.setter
    def new_split_plane(self, value: bool):
        if not isinstance(value, bool):
            raise TypeError(f"geouned.Options.new_split_plane should be a bool, not a {type(value)}")
        self._newSplitPlane = value

    @property
    def del_last_number(self):
        return self._delLastNumber

    @del_last_number.setter
    def del_last_number(self, value: bool):
        if not isinstance(value, bool):
            raise TypeError(f"geouned.Options.del_last_number should be a bool, not a {type(value)}")
        self._delLastNumber = value

    @property
    def enlarge_box(self):
        return self._enlargeBox

    @enlarge_box.setter
    def enlarge_box(self, value: Real):
        if not isinstance(value, Real):
            raise TypeError(f"geouned.Options.enlarge_box should be a Real, not a {type(value)}")
        if value < 0:
            raise ValueError(f"geouned.Options.enlarge_box should be above 0, not {value}")
        self._enlargeBox = value

    @property
    def n_plane_reverse(self):
        return self._nPlaneReverse

    @n_plane_reverse.setter
    def n_plane_reverse(self, value: int):
        if not isinstance(value, int):
            raise TypeError(f"geouned.Options.n_plane_reverse should be a int, not a {type(value)}")
        self._nPlaneReverse = value

    @property
    def split_tolerance(self):
        return self._splitTolerance

    @split_tolerance.setter
    def split_tolerance(self, value: Real):
        if not isinstance(value, Real):
            raise TypeError(f"geouned.Options.split_tolerance should be a Real, not a {type(value)}")
        if value < 0:
            raise ValueError(f"geouned.Options.split_tolerance should be above 0, not {value}")
        self._splitTolerance = value

    @property
    def scale_up(self):
        return self._scaleUp

    @scale_up.setter
    def scale_up(self, value: bool):
        if not isinstance(value, bool):
            raise TypeError(f"geouned.Options.scale_up should be a bool, not a {type(value)}")
        self._scaleUp = value

    @property
    def quadric_py(self):
        return self._quadricPY

    @quadric_py.setter
    def quadric_py(self, value: bool):
        if not isinstance(value, bool):
            raise TypeError(f"geouned.Options.quadric_py should be a bool, not a {type(value)}")
        self._quadricPY = value

    @property
    def facets(self):
        return self._Facets

    @facets.setter
    def facets(self, value: bool):
        if not isinstance(value, bool):
            raise TypeError(f"geouned.Options.facets should be a bool, not a {type(value)}")
        self._Facets = value

    @property
    def prnt_3p_plane(self):
        return self._prnt3PPlane

    @prnt_3p_plane.setter
    def prnt_3p_plane(self, value: bool):
        if not isinstance(value, bool):
            raise TypeError(f"geouned.Options.prnt_3p_plane should be a bool, not a {type(value)}")
        self._prnt3PPlane = value

    @property
    def forceNoOverlap(self):
        return self._forceNoOverlap

    @forceNoOverlap.setter
    def forceNoOverlap(self, value: bool):
        if not isinstance(value, bool):
            raise TypeError(f"geouned.Options.forceNoOverlap should be a bool, not a {type(value)}")
        self._forceNoOverlap = value


class Tolerances:
    """A class for containing tolerances values

    Args:
        relative_tol (bool, optional): _description_. Defaults to False.
        relative_precision (float, optional): relative precision. Defaults to 1.0e-6.
        value (float, optional): Tolerance in single value comparison. Defaults to 1.0e-6.
        distance (float, optional): General Distance Tolerance. Defaults to 1.0e-4.
        angle (float, optional): General Angle Tolerance. Defaults to 1.0e-4.
        pln_distance (float, optional): distance between planes equal planes if distance between parallel planes < 1e-4 cm. Defaults to 1.0e-4.
        pln_angle (float, optional): angle between axis. 1e-4 : planes separate each other 0.1mm each 1m. Defaults to 1.0e-4.
        cyl_distance (float, optional): distance between radius/center. Defaults to 1.0e-4.
        cyl_angle (float, optional): angle between axis. Defaults to 1.0e-4.
        sph_distance (float, optional): distance between radius/center. Defaults to 1.0e-4.
        kne_distance (float, optional): distance between apex. Defaults to 1.0e-4.
        kne_angle (float, optional): angle between semiangles/axis. Defaults to 1.0e-4.
        tor_distance (float, optional): distance between Major/Minor radii/center. Defaults to 1.0e-4.
        tor_angle (float, optional): angle between axis. Defaults to 1.0e-4.
        min_area (float, optional): minimum face area to consider in cell definition. Defaults to 1.0e-2.
    """

    def __init__(
        self,
        relative_tol: bool = False,
        relative_precision: float = 1.0e-6,
        value: float = 1.0e-6,
        distance: float = 1.0e-4,
        angle: float = 1.0e-4,
        pln_distance: float = 1.0e-4,
        pln_angle: float = 1.0e-4,
        cyl_distance: float = 1.0e-4,
        cyl_angle: float = 1.0e-4,
        sph_distance: float = 1.0e-4,
        kne_distance: float = 1.0e-4,
        kne_angle: float = 1.0e-4,
        tor_distance: float = 1.0e-4,
        tor_angle: float = 1.0e-4,
        min_area: float = 1.0e-2,
    ):

        self.relative_tol = relative_tol
        self.relative_precision = relative_precision
        self.value = value
        self.distance = distance
        self.angle = angle
        self.pln_distance = pln_distance
        self.pln_angle = pln_angle
        self.cyl_distance = cyl_distance
        self.cyl_angle = cyl_angle
        self.sph_distance = sph_distance
        self.kne_distance = kne_distance
        self.kne_angle = kne_angle
        self.tor_distance = tor_distance
        self.tor_angle = tor_angle
        self.min_area = min_area

    @property
    def relative_tol(self):
        return self._relativeTol

    @relative_tol.setter
    def relative_tol(self, value: bool):
        if not isinstance(value, bool):
            raise TypeError(f"geouned.Tolerances.relative_tol should be a bool, not a {type(value)}")
        self._relativeTol = value

    @property
    def relative_precision(self):
        return self._relativePrecision

    @relative_precision.setter
    def relative_precision(self, value: float):
        if not isinstance(value, float):
            raise TypeError(f"geouned.Tolerances.relative_precision should be a float, not a {type(value)}")
        self._relativePrecision = value

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, value: float):
        if not isinstance(value, float):
            raise TypeError(f"geouned.Tolerances.value should be a float, not a {type(value)}")
        self._value = value

    @property
    def distance(self):
        return self._distance

    @distance.setter
    def distance(self, distance: float):
        if not isinstance(distance, float):
            raise TypeError(f"geouned.Tolerances.distance should be a float, not a {type(distance)}")
        self._distance = distance

    @property
    def angle(self):
        return self._angle

    @angle.setter
    def angle(self, angle: float):
        if not isinstance(angle, float):
            raise TypeError(f"geouned.Tolerances.angle should be a float, not a {type(angle)}")
        self._angle = angle

    @property
    def pln_distance(self):
        return self._pln_distance

    @pln_distance.setter
    def pln_distance(self, pln_distance: float):
        if not isinstance(pln_distance, float):
            raise TypeError(f"geouned.Tolerances.pln_distance should be a float, not a {type(pln_distance)}")
        self._pln_distance = pln_distance

    @property
    def cyl_distance(self):
        return self._cyl_distance

    @cyl_distance.setter
    def cyl_distance(self, cyl_distance: float):
        if not isinstance(cyl_distance, float):
            raise TypeError(f"geouned.Tolerances.cyl_distance should be a float, not a {type(cyl_distance)}")
        self._cyl_distance = cyl_distance

    @property
    def cyl_angle(self):
        return self._cyl_angle

    @cyl_angle.setter
    def cyl_angle(self, cyl_angle: float):
        if not isinstance(cyl_angle, float):
            raise TypeError(f"geouned.Tolerances.cyl_angle should be a float, not a {type(cyl_angle)}")
        self._cyl_angle = cyl_angle

    @property
    def sph_distance(self):
        return self._sph_distance

    @sph_distance.setter
    def sph_distance(self, sph_distance: float):
        if not isinstance(sph_distance, float):
            raise TypeError(f"geouned.Tolerances.sph_distance should be a float, not a {type(sph_distance)}")
        self._sph_distance = sph_distance

    @property
    def pln_angle(self):
        return self._pln_angle

    @pln_angle.setter
    def pln_angle(self, pln_angle: float):
        if not isinstance(pln_angle, float):
            raise TypeError(f"geouned.Tolerances.pln_angle should be a float, not a {type(pln_angle)}")
        self._pln_angle = pln_angle

    @property
    def kne_distance(self):
        return self._kne_distance

    @kne_distance.setter
    def kne_distance(self, kne_distance: float):
        if not isinstance(kne_distance, float):
            raise TypeError(f"geouned.Tolerances.kne_distance should be a float, not a {type(kne_distance)}")
        self._kne_distance = kne_distance

    @property
    def kne_angle(self):
        return self._kne_angle

    @kne_angle.setter
    def kne_angle(self, kne_angle: float):
        if not isinstance(kne_angle, float):
            raise TypeError(f"geouned.Tolerances.kne_angle should be a float, not a {type(kne_angle)}")
        self._kne_angle = kne_angle

    @property
    def tor_distance(self):
        return self._tor_distance

    @tor_distance.setter
    def tor_distance(self, tor_distance: float):
        if not isinstance(tor_distance, float):
            raise TypeError(f"geouned.Tolerances.tor_distance should be a float, not a {type(tor_distance)}")
        self._tor_distance = tor_distance

    @property
    def tor_angle(self):
        return self._tor_angle

    @tor_angle.setter
    def tor_angle(self, tor_angle: float):
        if not isinstance(tor_angle, float):
            raise TypeError(f"geouned.Tolerances.tor_angle should be a float, not a {type(tor_angle)}")
        self._tor_angle = tor_angle

    @property
    def min_area(self):
        return self._min_area

    @min_area.setter
    def min_area(self, min_area: float):
        if not isinstance(min_area, float):
            raise TypeError(f"geouned.Tolerances.min_area should be a float, not a {type(min_area)}")
        self._min_area = min_area


class NumericFormat:
    """Numerical format options for each of the surface types.

    Args:
        p_abc (str, optional): Plane general a,b,c params. Defaults to "14.7e".
        p_d (str, optional): Plane general d params. Defaults to "14.7e".
        p_xyz (str, optional): PX/PY/PZ params. Defaults to "14.7e".
        s_r (str, optional): SO/SX/SY/SZ/S radius. Defaults to "14.7e".
        s_xyz (str, optional): SO/SX/SY/SZ/S center. Defaults to "14.7e".
        c_r (str, optional): Cylinder radius. Defaults to "12f".
        c_xyz (str, optional): Cylinder center. Defaults to "12f".
        k_xyz (str, optional): Cone apex. Defaults to "13.6e".
        k_tan2 (str, optional): Cone tan^2 value. Defaults to "12f".
        t_r (str, optional): Torus radii. Defaults to "14.7e".
        t_xyz (str, optional): Torus center. Defaults to "14.7e".
        gq_1_to_6 (str, optional): GQ 1 to 6 coefficients (order 2 x2,y2,z2,xy,...). Defaults to "18.15f".
        gq_7_to_9 (str, optional): GQ 7 to 9 coefficients (order 1 x,y,z). Defaults to "18.15f".
        gq_10 (str, optional): GQ 10 coefficient. Defaults to "18.15f".
    """

    def __init__(
        self,
        p_abc: str = "14.7e",
        p_d: str = "14.7e",
        p_xyz: str = "14.7e",
        s_r: str = "14.7e",
        s_xyz: str = "14.7e",
        c_r: str = "12f",
        c_xyz: str = "12f",
        k_xyz: str = "13.6e",
        k_tan2: str = "12f",
        t_r: str = "14.7e",
        t_xyz: str = "14.7e",
        gq_1_to_6: str = "18.15f",
        gq_7_to_9: str = "18.15f",
        gq_10: str = "18.15f",
    ):

        self.p_abc = p_abc
        self.p_d = p_d
        self.p_xyz = p_xyz
        self.s_r = s_r
        self.s_xyz = s_xyz
        self.c_r = c_r
        self.c_xyz = c_xyz
        self.k_xyz = k_xyz
        self.k_tan2 = k_tan2
        self.t_r = t_r
        self.t_xyz = t_xyz
        self.gq_1_to_6 = gq_1_to_6
        self.gq_7_to_9 = gq_7_to_9
        self.gq_10 = gq_10

    @property
    def p_abc(self):
        return self._p_abc

    @p_abc.setter
    def p_abc(self, p_abc: str):
        if not isinstance(p_abc, str):
            raise TypeError(f"geouned.Tolerances.p_abc should be a str, not a {type(p_abc)}")
        self._p_abc = p_abc

    @property
    def p_d(self):
        return self._p_d

    @p_d.setter
    def p_d(self, p_d: str):
        if not isinstance(p_d, str):
            raise TypeError(f"geouned.Tolerances.p_d should be a str, not a {type(p_d)}")
        self._p_d = p_d

    @property
    def p_xyz(self):
        return self._p_xyz

    @p_xyz.setter
    def p_xyz(self, p_xyz: str):
        if not isinstance(p_xyz, str):
            raise TypeError(f"geouned.Tolerances.p_xyz should be a str, not a {type(p_xyz)}")
        self._p_xyz = p_xyz

    @property
    def s_r(self):
        return self._s_r

    @s_r.setter
    def s_r(self, s_r: str):
        if not isinstance(s_r, str):
            raise TypeError(f"geouned.Tolerances.s_r should be a str, not a {type(s_r)}")
        self._s_r = s_r

    @property
    def s_xyz(self):
        return self._s_xyz

    @s_xyz.setter
    def s_xyz(self, s_xyz: str):
        if not isinstance(s_xyz, str):
            raise TypeError(f"geouned.Tolerances.s_xyz should be a str, not a {type(s_xyz)}")
        self._s_xyz = s_xyz

    @property
    def c_r(self):
        return self._c_r

    @c_r.setter
    def c_r(self, c_r: str):
        if not isinstance(c_r, str):
            raise TypeError(f"geouned.Tolerances.c_r should be a str, not a {type(c_r)}")
        self._c_r = c_r

    @property
    def c_xyz(self):
        return self._c_xyz

    @c_xyz.setter
    def c_xyz(self, c_xyz: str):
        if not isinstance(c_xyz, str):
            raise TypeError(f"geouned.Tolerances.c_xyz should be a str, not a {type(c_xyz)}")
        self._c_xyz = c_xyz

    @property
    def k_xyz(self):
        return self._k_xyz

    @k_xyz.setter
    def k_xyz(self, k_xyz: str):
        if not isinstance(k_xyz, str):
            raise TypeError(f"geouned.Tolerances.k_xyz should be a str, not a {type(k_xyz)}")
        self._k_xyz = k_xyz

    @property
    def k_tan2(self):
        return self._k_tan2

    @k_tan2.setter
    def k_tan2(self, k_tan2: str):
        if not isinstance(k_tan2, str):
            raise TypeError(f"geouned.Tolerances.k_tan2 should be a str, not a {type(k_tan2)}")
        self._k_tan2 = k_tan2

    @property
    def t_r(self):
        return self._t_r

    @t_r.setter
    def t_r(self, t_r: str):
        if not isinstance(t_r, str):
            raise TypeError(f"geouned.Tolerances.t_r should be a str, not a {type(t_r)}")
        self._t_r = t_r

    @property
    def t_xyz(self):
        return self._t_xyz

    @t_xyz.setter
    def t_xyz(self, t_xyz: str):
        if not isinstance(t_xyz, str):
            raise TypeError(f"geouned.Tolerances.t_xyz should be a str, not a {type(t_xyz)}")
        self._t_xyz = t_xyz

    @property
    def gq_1_to_6(self):
        return self._gq_1_to_6

    @gq_1_to_6.setter
    def gq_1_to_6(self, gq_1_to_6: str):
        if not isinstance(gq_1_to_6, str):
            raise TypeError(f"geouned.Tolerances.gq_1_to_6 should be a str, not a {type(gq_1_to_6)}")
        self._gq_1_to_6 = gq_1_to_6

    @property
    def gq_7_to_9(self):
        return self._gq_7_to_9

    @gq_7_to_9.setter
    def gq_7_to_9(self, gq_7_to_9: str):
        if not isinstance(gq_7_to_9, str):
            raise TypeError(f"geouned.Tolerances.gq_7_to_9 should be a str, not a {type(gq_7_to_9)}")
        self._gq_7_to_9 = gq_7_to_9

    @property
    def gq_10(self):
        return self._gq_10

    @gq_10.setter
    def gq_10(self, gq_10: str):
        if not isinstance(gq_10, str):
            raise TypeError(f"geouned.Tolerances.gq_10 should be a str, not a {type(gq_10)}")
        self._gq_10 = gq_10


class Settings:
    """Settings for changing the way the CAD to CSG conversion is done

    Args:
        step_filename (str, optional): Name of the CAD file (in STEP format) to
            be converted. Defaults to "".
        mat_file (str, optional): _description_. Defaults to "".
        void_gen (bool, optional): Generate voids of the geometry. Defaults
            to True.
        debug (bool, optional): Write step files of original and decomposed
            solids, for each solid in the STEP file. Defaults to False.
        comp_solids (bool, optional): Join subsolids of STEP file as a single
            compound solid. Step files generated with SpaceClaim have not
            exactly the same level of solids as FreeCAD. It may a happened
            that solids defined has separated solids are read by FreeCAD
            as a single compound solid (and will produce only one MCNP
            cell). In this case comp_solids should be set to False. Defaults
            to True.
        simplify (str, optional): Simplify the cell definition considering
            relative surfaces position and using Boolean logics. Available
            options are: "no" no optimization, "void" only void cells are
            simplified. Algorithm is faster but the simplification is not
            optimal. "voidfull" : only void cells are simplified with the
            most optimal algorithm. The time of the conversion can be
            multiplied by 5 or more. "full" : all the cells (solids and
            voids) are simplified. Defaults to "No".
        cell_range (list, optional): Range of cell to be converted (only one
            range is allowed, e.g [100,220]). Default all solids are
            converted. Defaults to [].
        export_solids (str, optional): Export CAD solid after reading.
            The execution is stopped after export, the translation is not
            carried out. Defaults to "".
        min_void_size (float, optional): Minimum size of the edges of the
            void cell. Units are in mm. Defaults to 200.0.
        max_surf (int, optional): #TODO
        max_bracket (int, optional): Maximum number of brackets (solid
            complementary) allowed in void cell definition. Defaults to 30.
        void_mat (list, optional): Assign a material defined by the user
            instead of void for cells without material definition and the
            cells generated in the automatic void generation. The format
            is a 3 valued tuple (mat_label, mat_density, mat_description).
            Example (100,1e-3,'Air assigned to Void'). Defaults to [].
        void_exclude (list, optional): #TODO see issue 87. Defaults to [].
        start_cell (int, optional): Starting cell numbering label. Defaults to 1.
        start_surface (int, optional): Starting surface numbering label. Defaults to 1.
        sort_enclosure (bool, optional): If enclosures are defined in the
            CAD models, the voids cells of the enclosure will be located in
            the output file in the same location where the enclosure solid
            is located in the CAD solid tree.. Defaults to False.
    """

    def __init__(
        self,
        mat_file: str = "",
        void_gen: bool = True,
        debug: bool = False,
        comp_solids: bool = True,
        simplify: str = "no",
        cell_range: list = [],
        export_solids: str = "",
        min_void_size: float = 200.0,  # units mm
        max_surf: int = 50,
        max_bracket: int = 30,
        void_mat: list = [],
        void_exclude: list = [],
        start_cell: int = 1,
        start_surface: int = 1,
        sort_enclosure: bool = False,
    ):

        self.mat_file = mat_file
        self.void_gen = void_gen
        self.debug = debug
        self.comp_solids = comp_solids
        self.simplify = simplify
        self.cell_range = cell_range
        self.export_solids = export_solids
        self.min_void_size = min_void_size
        self.max_surf = max_surf
        self.max_bracket = max_bracket
        self.void_mat = void_mat
        self.void_exclude = void_exclude
        self.start_cell = start_cell
        self.start_surface = start_surface
        self.sort_enclosure = sort_enclosure

    @property
    def mat_file(self):
        return self._mat_file

    @mat_file.setter
    def mat_file(self, mat_file: str):
        if not isinstance(mat_file, str):
            raise TypeError(f"geouned.Tolerances.mat_file should be a str, not a {type(mat_file)}")
        self._mat_file = mat_file

    @property
    def void_gen(self):
        return self._void_gen

    @void_gen.setter
    def void_gen(self, void_gen: bool):
        if not isinstance(void_gen, bool):
            raise TypeError(f"geouned.Tolerances.void_gen should be a bool, not a {type(void_gen)}")
        self._void_gen = void_gen

    @property
    def debug(self):
        return self._debug

    @debug.setter
    def debug(self, debug: bool):
        if not isinstance(debug, bool):
            raise TypeError(f"geouned.Tolerances.debug should be a bool, not a {type(debug)}")
        self._debug = debug

    @property
    def comp_solids(self):
        return self._comp_solids

    @comp_solids.setter
    def comp_solids(self, comp_solids: bool):
        if not isinstance(comp_solids, bool):
            raise TypeError(f"geouned.Tolerances.comp_solids should be a bool, not a {type(comp_solids)}")
        self._comp_solids = comp_solids

    @property
    def simplify(self):
        return self._simplify

    @simplify.setter
    def simplify(self, simplify: str):
        if not isinstance(simplify, str):
            raise TypeError(f"geouned.Tolerances.simplify should be a str, not a {type(simplify)}")
        self._simplify = simplify

    @property
    def cell_range(self):
        return self._cell_range

    @cell_range.setter
    def cell_range(self, cell_range: list):
        if not isinstance(cell_range, list):
            raise TypeError(f"geouned.Tolerances.cell_range should be a list, not a {type(cell_range)}")
        for entry in cell_range:
            if not isinstance(entry, int):
                raise TypeError(f"geouned.Tolerances.cell_range should be a list of ints, not a {type(entry)}")
        self._cell_range = cell_range

    @property
    def export_solids(self):
        return self._export_solids

    @export_solids.setter
    def export_solids(self, export_solids: str):
        if not isinstance(export_solids, str):
            raise TypeError(f"geouned.Tolerances.export_solids should be a str, not a {type(export_solids)}")
        self._export_solids = export_solids

    @property
    def min_void_size(self):
        return self._min_void_size

    @min_void_size.setter
    def min_void_size(self, min_void_size: float):
        if not isinstance(min_void_size, float):
            raise TypeError(f"geouned.Tolerances.min_void_size should be a float, not a {type(min_void_size)}")
        self._min_void_size = min_void_size

    @property
    def max_surf(self):
        return self._max_surf

    @max_surf.setter
    def max_surf(self, max_surf: int):
        if not isinstance(max_surf, int):
            raise TypeError(f"geouned.Tolerances.max_surf should be a int, not a {type(max_surf)}")
        self._max_surf = max_surf

    @property
    def max_bracket(self):
        return self._max_bracket

    @max_bracket.setter
    def max_bracket(self, max_bracket: int):
        if not isinstance(max_bracket, int):
            raise TypeError(f"geouned.Tolerances.max_bracket should be a int, not a {type(max_bracket)}")
        self._max_bracket = max_bracket

    @property
    def void_mat(self):
        return self._void_mat

    @void_mat.setter
    def void_mat(self, void_mat: list):
        if not isinstance(void_mat, list):
            raise TypeError(f"geouned.Tolerances.void_mat should be a list, not a {type(void_mat)}")
        for entry in void_mat:
            if not isinstance(entry, int):
                raise TypeError(f"geouned.Tolerances.void_mat should be a list of ints, not a {type(entry)}")
        self._void_mat = void_mat

    @property
    def void_exclude(self):
        return self._void_exclude

    @void_exclude.setter
    def void_exclude(self, void_exclude: list):
        if not isinstance(void_exclude, list):
            raise TypeError(f"geouned.Tolerances.void_exclude should be a list, not a {type(void_exclude)}")
        for entry in void_exclude:
            if not isinstance(entry, int):
                raise TypeError(f"geouned.Tolerances.void_exclude should be a list of ints, not a {type(entry)}")
        self._void_exclude = void_exclude

    @property
    def start_cell(self):
        return self._start_cell

    @start_cell.setter
    def start_cell(self, start_cell: int):
        if not isinstance(start_cell, int):
            raise TypeError(f"geouned.Tolerances.start_cell should be a int, not a {type(start_cell)}")
        self._start_cell = start_cell

    @property
    def start_surface(self):
        return self._start_surface

    @start_surface.setter
    def start_surface(self, start_surface: int):
        if not isinstance(start_surface, int):
            raise TypeError(f"geouned.Tolerances.start_surface should be a int, not a {type(start_surface)}")
        self._start_surface = start_surface

    @property
    def sort_enclosure(self):
        return self._sort_enclosure

    @sort_enclosure.setter
    def sort_enclosure(self, sort_enclosure: bool):
        if not isinstance(sort_enclosure, bool):
            raise TypeError(f"geouned.Tolerances.sort_enclosure should be a bool, not a {type(sort_enclosure)}")
        self._sort_enclosure = sort_enclosure
