##############################
# Module to write MCNP input #
##############################

import logging

import FreeCAD

from ..CodeVersion import *
from ..Utils.Functions import SurfacesDict
from .Functions import change_surf_sign, open_mc_surface, write_openmc_region

logger = logging.getLogger("general_logger")


class OpenmcInput:
    def __init__(self, Meta, Surfaces, options, tolerances, numeric_format):

        self.Cells = Meta
        self.options = options
        self.tolerances = tolerances
        self.numeric_format = numeric_format

        self.__get_surface_table__()
        self.__simplify_planes__(Surfaces)

        self.Surfaces = self.__sorted_surfaces__(Surfaces)
        self.Materials = set()

    def write_xml(self, filename):
        logger.info(f"write OpenMC xml file {filename}")
        with open(file=filename, mode="w", encoding="utf-8") as self.inpfile:
            self.__write_xml_header__()
            self.inpfile.write("<geometry>\n")
            self.__write_xml_cell_block__()
            self.inpfile.write(" \n")
            self.__write_xml_surface_block__()
            self.inpfile.write("</geometry>\n")
        return

    def __write_xml_header__(self):
        Header = "<?xml version='1.0' encoding='UTF-8'?>\n"
        self.inpfile.write(Header)
        return

    def __write_xml_cell_block__(self):
        for _, cell in enumerate(self.Cells):
            if cell.MatInfo == "Graveyard":
                continue
            self.__write_xml_cells__(cell)
        return

    def __write_xml_surface_block__(self):
        for surf in self.Surfaces[:-1]:
            self.__write_xml_surfaces__(surf)
        self.__write_xml_surfaces__(self.Surfaces[-1], boundary=True)

    def __write_xml_cells__(self, cell):
        """Write the cell in xml OpenMC format"""
        index = cell.label
        cellName = ". ".join(cell.Comments.splitlines())
        if cell.__id__ is None:
            return

        if cell.Material == 0:
            matName = "void"
        else:
            matName = f"{cell.Material}"

        OMCcell = '  <cell id="{}" material="{}" name="{}" region="{}" universe="1"/>\n'.format(
            index,
            matName,
            cellName,
            write_openmc_region(cell.Definition, self.options, "XML"),
        )
        self.inpfile.write(OMCcell)
        return

    def __write_xml_surfaces__(self, surface, boundary=False):
        """Write the surfaces in xml OpenMC format"""

        surfType, coeffs = open_mc_surface(
            surface.Type, surface.Surf, self.tolerances, self.numeric_format
        )

        if not boundary:
            OMCsurf = '  <surface id="{}" type="{}" coeffs="{}" />\n'.format(
                surface.Index, surfType, coeffs
            )
        else:
            OMCsurf = '  <surface id="{}" type="{}" coeffs="{}" boundary="vacuum" />\n'.format(
                surface.Index, surfType, coeffs
            )

        self.inpfile.write(OMCsurf)
        return

    def write_py(self, filename):
        logger.info(f"write OpenMC python script {filename}")

        # get all the materials present in the model
        for cell in self.Cells:
            if cell.Material != 0:
                self.Materials.add(cell.Material)

        with open(file=filename, mode="w", encoding="utf-8") as self.inpfile:
            self.__write_py_header__()

            if len(self.Materials) > 0:
                self.inpfile.write("# Materials setup\n")
                self.__write_py_materials__()
            self.inpfile.write("\n")

            self.inpfile.write("# Surface setup\n")
            self.__write_py_surface_block__()
            self.inpfile.write("\n")

            self.inpfile.write("# Cell definition \n")
            self.__write_py_cell_block__()
        return

    def __write_py_header__(self):

        Header = """\
# openMC geometry script generated by GEOUNED
import openmc

###############################################################################
# Define problem geometry
###############################################################################

"""
        self.inpfile.write(Header)
        return

    def __write_py_materials__(self):
        matList = tuple(sorted(self.Materials))
        strMat = []
        for m in matList:
            material = f"M{m} = openmc.Material(name='M{m}')\n"
            self.inpfile.write(material)
            strMat.append(f"M{m}")

        collect = f"materials = openmc.Materials([{', '.join(strMat)}])\n"
        self.inpfile.write(collect)
        self.inpfile.write("materials.export_to_xml()\n")

    def __write_py_surface_block__(self):

        for surf in self.Surfaces[:-1]:
            self.__write_py_surfaces__(surf)
        self.__write_py_surfaces__(self.Surfaces[-1], boundary=True)

    def __write_py_surfaces__(self, surface, boundary=False):
        """Write the surfaces in python OpenMC format"""

        surfType, coeffs = open_mc_surface(
            surface.Type,
            surface.Surf,
            self.tolerances,
            self.numeric_format,
            out_xml=False,
            quadricForm=self.options.quadricPY,
        )

        if not boundary:
            OMCsurf = f"S{surface.Index} = openmc.{surfType}({coeffs})\n"
        else:
            OMCsurf = 'S{} = openmc.{}({}, boundary_type="vacuum")\n'.format(
                surface.Index, surfType, coeffs
            )

        self.inpfile.write(OMCsurf)
        return

    def __write_py_cell_block__(self):

        cellNames = []
        for i, cell in enumerate(self.Cells):
            if cell.MatInfo == "Graveyard":
                continue
            self.__write_py_cells__(cell)
            if cell.__id__ is None:
                continue
            cellNames.append(f"C{cell.label}")

        geometry = (
            "\ngeometry = openmc.Geometry([{}])\ngeometry.export_to_xml()\n".format(
                ", ".join(cellNames)
            )
        )

        self.inpfile.write(geometry)
        return

    def __write_py_cells__(self, cell):
        """Write the cell in python OpenMC format"""
        index = cell.label
        cellName = ". ".join(cell.Comments.splitlines())
        if cell.__id__ is None:
            return

        if cell.Material == 0:
            OMCcell = 'C{} = openmc.Cell(name="{}", region={})\n'.format(
                index,
                cellName,
                write_openmc_region(cell.Definition, self.options, "PY"),
            )
        else:
            matName = f"M{cell.Material}"
            OMCcell = 'C{} = openmc.Cell(name="{}", fill={}, region={})\n'.format(
                index,
                cellName,
                matName,
                write_openmc_region(cell.Definition, self.options, "PY"),
            )
        self.inpfile.write(OMCcell)
        return

    def __get_surface_table__(self):
        self.surfaceTable = {}
        self.__solidCells__ = 0
        self.__cells__ = 0
        self.__materials__ = set()

        for i, CellObj in enumerate(self.Cells):
            if CellObj.__id__ is None:
                continue
            self.__cells__ += 1
            if CellObj.Material != 0:
                self.__materials__.add(CellObj.Material)

            surf = CellObj.Definition.get_surfaces_numbers()
            if not CellObj.Void:
                self.__solidCells__ += 1
            for index in surf:
                if index in self.surfaceTable.keys():
                    self.surfaceTable[index].add(i)
                else:
                    self.surfaceTable[index] = {i}
        return

    def __simplify_planes__(self, Surfaces):

        for p in Surfaces["PX"]:
            if p.Surf.Axis[0] < 0:
                p.Surf.Axis = FreeCAD.Vector(1, 0, 0)
                self.__change_surf_sign__(p)

        for p in Surfaces["PY"]:
            if p.Surf.Axis[1] < 0:
                p.Surf.Axis = FreeCAD.Vector(0, 1, 0)
                self.__change_surf_sign__(p)

        for p in Surfaces["PZ"]:
            if p.Surf.Axis[2] < 0:
                p.Surf.Axis = FreeCAD.Vector(0, 0, 1)
                self.__change_surf_sign__(p)
        return

    def __sorted_surfaces__(self, Surfaces):
        temp = SurfacesDict(Surfaces)
        surfList = []
        for ind in range(
            Surfaces.IndexOffset, Surfaces.surfaceNumber + Surfaces.IndexOffset
        ):
            s = temp.get_surface(ind + 1)
            if s is not None:
                surfList.append(s)
                temp.del_surface(ind + 1)
        return surfList

    def __change_surf_sign__(self, p):

        if p.Index not in self.surfaceTable.keys():
            logger.info(
                f"{p.Type} Surface {p.Index} not used in cell definition {p.Surf.Axis} {p.Surf.Position}"
            )
            return

        for ic in self.surfaceTable[p.Index]:
            surf = self.Cells[ic].Definition.get_surfaces_numbers()
            for s in surf:
                if s == p.Index:
                    change_surf_sign(s, self.Cells[ic].Definition)
