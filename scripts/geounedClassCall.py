#!/usr/bin/python

# Path to GEOUNED Package


# only if modules are not in the PYTHONPATH or directly installed in the python distribution used
import sys

geo_path = "C:\\Users\\Patrick\\Documents\\GitHub\\GEOUNED\\src"
sys.path.append(geo_path)
sys.path.append("C:\\Program Files\\FreeCAD 0.19\\bin...")

from geouned import CadToCsg

stepFileName = "placa.stp"

GEO = CadToCsg("Conversion Example")

GEO.set("step_filename", stepFileName)
GEO.set("geometryName", "Placa")
GEO.set("outFormat", ("mcnp", "openmc_xml"))
GEO.set("planeDistance", 0.05)
GEO.set("quadric_py", True)
GEO.set("p_abc", "12f")
GEO.set("p_d", "12f")
GEO.Start()
