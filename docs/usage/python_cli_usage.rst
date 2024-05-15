Command Line Tool Usage
=======================

GEOUNED can be used in the command line.

These examples assumes you have a CAD STEP file in the current working directory of the terminal called "cuboid.stp"

The most minimal use case below shows a minimal config.json file being used.

First create a JSON file called "config.json" containing the following.

.. code-block:: json

    {
        "step_filename": "cuboid.stp"
    }

Then execute the command line interface tool to convert your STEP file to CSG files with the default configuration.

.. code-block:: bash

    geouned_cadtocsg -i config.json

The following example shows a usage with every attributes specified in the config.json file.

The contents of the JSON file closely matches the Class arguments and method arguments when using the Python package.

For a full description of each keyword see the `Python API reference section <../python_api.html>`_ of the documentation.

Here is a complete JSON file specification

.. code-block:: json

    {
        "step_filename": "cuboid.stp",
        "Options": {
            "force_cylinder": false,
            "new_split_plane": true,
            "del_last_number": false,
            "enlarge_box": 2.0,
            "n_plane_reverse": 0,
            "split_tolerance": 0.0,
            "scale_up": true,
            "quadric_py": false,
            "facets": false,
            "prnt_3p_plane": false,
            "forceNoOverlap": false
        },
        "Tolerances": {
            "relative_tol": false,
            "relative_precision": 1e-06,
            "value": 1e-06,
            "distance": 0.0001,
            "angle": 0.0001,
            "pln_distance": 0.0001,
            "pln_angle": 0.0001,
            "cyl_distance": 0.0001,
            "cyl_angle": 0.0001,
            "sph_distance": 0.0001,
            "kne_distance": 0.0001,
            "kne_angle": 0.0001,
            "tor_distance": 0.0001,
            "tor_angle": 0.0001,
            "min_area": 0.01
        },
        "NumericFormat": {
            "p_abc": "14.7e",
            "p_d": "14.7e",
            "p_xyz": "14.7e",
            "s_r": "14.7e",
            "s_xyz": "14.7e",
            "c_r": "12f",
            "c_xyz": "12f",
            "k_xyz": "13.6e",
            "k_tan2": "12f",
            "t_r": "14.7e",
            "t_xyz": "14.7e",
            "gq_1_to_6": "18.15f",
            "gq_7_to_9": "18.15f",
            "gq_10": "18.15f"
        },
        "Settings": {
            "mat_file": "",
            "void_gen": true,
            "debug": false,
            "comp_solids": true,
            "simplify": "no",
            "cell_range": [],
            "export_solids": "",
            "min_void_size": 200.0,
            "max_surf": 50,
            "max_bracket": 30,
            "void_mat": [],
            "void_exclude": [],
            "start_cell": 1,
            "start_surface": 1,
            "sort_enclosure": false
        },
        "export_csg":{
            "title": "Converted with GEOUNED",
            "geometryName": "csg",
            "outFormat": ["openmc_xml", "openmc_py", "serpent", "phits", "mcnp"],
            "volSDEF": false,
            "volCARD": true,
            "dummyMat": false,
            "cellCommentFile": false,
            "cellSummaryFile": true
        }
    }

This is converted in the same way as the minimal JSON config file

.. code-block:: bash

    geouned_cadtocsg -i config.json
