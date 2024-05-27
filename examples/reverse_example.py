
# import geouned
# geo=geouned.CsgToCad()
# geo.export_cad(
#     input_filename="/home/j/GEOUNED/tests_outputs/testing/inputSTEP/placa/placa.xml",
#     output_filename="example_reversed_cad.step",
#     csg_format="openmc_xml",
# )

import geouned
geo=geouned.CsgToCad()
geo.export_cad(
    input_filename="/home/j/GEOUNED/tests_outputs/testing/inputSTEP/placa/placa.mcnp",
    output_filename="example.step",
    csg_format="mcnp",
)

