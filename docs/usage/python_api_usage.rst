Python Package Usage
====================

The Python API has two main classes.
The first main class is ``CadToCsg()`` which converts CAD geometry to Constructive Solid Geometry (CSG).
There are many arguments that can be passed into the ``CadToCsg()`` class which are documented in the `Python API reference section <../python_api.html>`_ of the documentation.


If you have install GEOUNED and FreeCAD into your system Python then you can simply run a .py script with Python.
The most minimal use case below shows GEOUNED being imported and the CadToCsg being used to convert a STEP CAD file called 'cuboid.stp' into a vanity of CSG format. 
The example makes use of default  attributes.

.. code-block:: python

    import geouned
    geo = geouned.CadToCsg(step_filename='cuboid.stp')
    geo.start()
    geo.export_csg()

Users can change :meth:`geouned.Options`, :meth:`geouned.Settings`, :meth:`geouned.Tolerances` and :meth:`geouned.NumericFormat` to suit the conversion desired.
The following example shows a usage with every attributes specified.

.. code-block:: python

    import geouned

    my_options = geouned.Options(
        force_cylinder=False,
        new_split_plane=True,
        del_last_number=False,
        enlarge_box=2,
        n_plane_reverse=0,
        split_tolerance=0,
        scale_up=True,
        quadric_py=False,
        facets=False,
        prnt_3p_plane=False,
        forceNoOverlap=False,
    )

    my_settings = geouned.Settings(
        mat_file="",
        void_gen=True,
        debug=False,
        comp_solids=True,
        simplify="no",
        cell_range=[],
        export_solids="",
        min_void_size=200.0,
        max_surf=50,
        max_bracket=30,
        void_mat=[],
        void_exclude=[],
        start_cell=1,
        start_surface=1,
        sort_enclosure=False,
    )

    my_tolerances = geouned.Tolerances(
        relative_tol=False,
        relative_precision=0.000001,
        value=0.000001,
        distance=0.0001,
        angle=0.0001,
        pln_distance=0.0001,
        pln_angle=0.0001,
        cyl_distance=0.0001,
        cyl_angle=0.0001,
        sph_distance=0.0001,
        kne_distance=0.0001,
        kne_angle=0.0001,
        tor_distance=0.0001,
        tor_angle=0.0001,
        min_area=0.01,
    )
    my_numeric_format = geouned.NumericFormat(
        p_abc="14.7e",
        p_d="14.7e",
        p_xyz="14.7e",
        s_r="14.7e",
        s_xyz="14.7e",
        c_r="12f",
        c_xyz="12f",
        k_xyz="13.6e",
        k_tan2="12f",
        t_r="14.7e",
        t_xyz="14.7e",
        gq_1_to_6="18.15f",
        gq_7_to_9="18.15f",
        gq_10="18.15f",
    )

    geo = geouned.CadToCsg(
        step_filename="cuboid.stp,
        options=my_options,
        settings=my_settings,
        tolerances=my_tolerances,
        numeric_format=my_numeric_format,
    )

    geo.start()

    geo.export_csg(
        title="Converted with GEOUNED",
        geometryName="csg",
        outFormat=(
            "openmc_xml",
            "openmc_py",
            "serpent",
            "phits",
            "mcnp",
        ),
        volSDEF=True,
        volCARD=False,
        UCARD=None,
        dummyMat=True,
        cellCommentFile=False,
        cellSummaryFile=False,
    )
