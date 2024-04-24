#################################
# Module to write Serpent input #
#################################
from datetime import datetime

import FreeCAD

from ..code_version import *
from ..utils.basic_functions_part1 import is_opposite, points_to_coeffs
from ..utils.functions import Surfaces_dict
from ..utils.options.classes import Options as opt
from .functions import SerpentSurface, change_surf_sign, write_serpent_cell_def


class SerpentInput:
    def __init__(
        self, Meta, Surfaces, step_file, title, vol_sdef, vol_card, ucard, dummy_mat
    ):
        self.title = title
        self.vol_sdef = vol_sdef
        self.vol_card = vol_card
        self.ucard = ucard
        self.dummy_mat = dummy_mat
        self.Cells = Meta
        self.Options = {
            "Volume": self.vol_card,
            "Particle": ("n", "p"),
            "Universe": self.ucard,
        }
        self.part = "p"
        self.step_file = step_file

        if isinstance(self.step_file, (tuple, list)):
            self.step_file = "; ".join(self.step_file)

        if self.title == "":
            self.title = self.step_file

        self.__getSurfaceTable__()
        self.__simplifyPlanes__(Surfaces)

        self.Surfaces = self.__sortedSurfaces__(Surfaces)
        self.Materials = set()

        return

    # def setSRC(self,data):

    #     if data[0] is not None:
    #       sdef     = f'src point {self.part}\n'
    #       self.src_sphere = (sdef)
    #     else:
    #       self.src_sphere = None
    #     xmin,xmax,ymin,ymax,zmin,zmax = data[1]

    #     sdef = 'SDEF PAR={} X=D1 Y=D2 Z=D3 \n'.format(self.part)
    #     SI1  = 'SI1 {:13.7e} {:13.7e} \n'.format(xmin*0.1,xmax*0.1)
    #     SI2  = 'SI2 {:13.7e} {:13.7e} \n'.format(ymin*0.1,ymax*0.1)
    #     SI3  = 'SI3 {:13.7e} {:13.7e} \n'.format(zmin*0.1,zmax*0.1)
    #     SP1  = 'SP1 0  1 \n'
    #     SP2  = 'SP2 0  1 \n'
    #     SP3  = 'SP3 0  1 \n'
    #     self.SDEF_box    = (sdef,SI1,SI2,SI3,SP1,SP2,SP3)

    def write_input(self, filename):
        print(f"write Serpent file {filename}")
        self.inpfile = open(filename, "w", encoding="utf-8")
        self.__write_header__(filename)
        cellblockHeader = """\
% --- CELL DEFINITIONS 
"""
        self.inpfile.write(cellblockHeader)
        self.__write_cell_block__()
        self.inpfile.write(" \n")

        surfaceHeader = """\
% --- SURFACE DEFINITIONS 
"""
        self.inpfile.write(surfaceHeader)
        self.__write_surface_block__()
        self.inpfile.write(" \n")

        self.__write_source_block__()

        self.inpfile.close()
        return

    def __write_header__(self, fileout):

        version = GEOUNED_Version
        releaseDate = GEOUNED_ReleaseDate
        freeCAD_Version = "{V[0]:}.{V[1]:}.{V[2]:}".format(V=FreeCAD.Version())

        Header = """{}
%   ______ _______  _____      _     _ __   _ _______ ______  
%  |  ____ |______ |     | ___ |     | | \  | |______ |     \ 
%  |_____| |______ |_____|     |_____| |  \_| |______ |_____/
% Version : {}     {}
% FreeCAD Version : {} \n""".format(
            self.title, version, releaseDate, freeCAD_Version
        )

        Information = """%
% *************************************************************
% Original Step file : {}
%
% Creation Date : {}
% Solid Cells   : {}
% Total Cells   : {}
% Surfaces      : {}
% Materials     : {}
%
% **************************************************************\n""".format(
            self.step_file,
            datetime.now(),
            self.__solidCells__,
            self.__cells__,
            len(self.Surfaces),
            len(self.__materials__),
        )
        self.inpfile.write(Header)
        self.inpfile.write(Information)
        return

    def __write_surface_block__(self):

        for surf in self.Surfaces:
            self.__write_surfaces__(surf)

    def __write_cell_block__(self):

        for i, cell in enumerate(self.Cells):
            self.__write_cells__(cell)
        return

    # to_mod
    def __write_cells__(self, cell):
        index = cell.label

        # If index is None, the object does not contain cell definition
        # but a comment to insert between cells
        if cell.__id__ is None:
            comment = self.comment_line(cell.Comments)
            self.inpfile.write(comment)
            return

        if self.Options["Universe"] is not None:
            if cell.Material == 0:
                cellHeader = (
                    f'cell {index:<5d} {self.Options["Universe"]} {"void":<5d}  '
                )
            else:
                self.Materials.add(cell.Material)
                cellHeader = (
                    f'cell {index:<5d} {self.Options["Universe"]} {cell.Material:<5d} '
                )

            serpent_cell = (
                f"{cellHeader}{self.__cellFormat__(cell.Definition, offset=len(cellHeader))}"
                f"{self.comment_format(cell.Comments, cell.MatInfo)}"
            )
            self.inpfile.write(serpent_cell)
        else:
            if cell.Material == 0 and cell.MatInfo != "Graveyard":
                cellHeader = f"cell {index:<5d} 0 void "
            elif cell.Material == 0 and cell.MatInfo == "Graveyard":
                cellHeader = f"cell {index:<5d} 0 outside "
            else:
                self.Materials.add(cell.Material)
                cellHeader = f"cell {index:<5d}  0  {cell.Material:<5d} "

        serpent_cell = (
            f"{cellHeader}{self.__cellFormat__(cell.Definition, offset=len(cellHeader))}"
            f"{self.comment_format(cell.Comments, cell.MatInfo)}"
        )
        self.inpfile.write(serpent_cell)

        return

    def __write_surfaces__(self, surface):
        """Write the surfaces in Serpent format"""

        Serpent_def = SerpentSurface(surface.Index, surface.Type, surface.Surf)
        if Serpent_def:
            Serpent_def += "\n"
            self.inpfile.write(Serpent_def)
        else:
            print("Surface {} cannot be written in Serpent input".format(surface.Type))
        return

    # No void all option in Serpent. For now remove addition of source.

    def __write_source_block__(self):

        #       if self.SDEF_sphere is None:  return
        MODE = f"\nset nps 1e6\nset bc 1"
        if self.dummy_mat:
            mat = list(self.Materials)
            mat.sort()
            MATCARD = ""
            for m in mat:
                MATCARD += "mat {:<6d} {:11.4e} \n1001 1 \n".format(
                    m, self.cell.Density
                )
            Block = MATCARD + "% \n" + MODE
        else:
            Block = MODE

        # Not included for Serpent
        #    if self.vol_sdef:
        #       Block += 'PRDMP 2J -1\n'
        #       for line in self.SDEF_box :
        #          Block += 'C ' + line
        #       for line in self.SDEF_sphere:
        #          Block += line

        #       celList,volList = self.__get_solidCellVolume__()

        #       F4Tally = CardLine('F4:{} '.format(self.part))
        #       F4Tally.extend(celList)
        #       SD4     = CardLine('SD4  ',fmt='13.7e')
        #       SD4.extend(volList)

        #       Block += F4Tally.str+'\n'
        #       if not self.vol_card :
        #          Block += SD4.str
        #       else :
        #          Block += 'C Cell volume normalization is set in cell cards VOL\n'
        #    else :
        #       for line in self.SDEF_sphere :
        #          Block += 'C ' + line
        #       for line in self.SDEF_box :
        #          Block +=  line

        self.inpfile.write(Block)

    def __cellFormat__(self, Definition, offset=11):
        return write_serpent_cell_def(Definition, tabspace=11, offset=offset)

    # Function not relevant for Serpent : No importance setting, universes assigned elsewhere.
    # Volumes only defined on tally cards.
    # def __optionFormat__(self,cell):

    #    option = ''
    #    if self.Options['Volume']:
    #       if not cell.Void :
    #         option ='{:11s}Vol={:e}\n'.format('',cell.Volume*1e-3)
    #       else:
    #         option = '{:11s}Vol=1.0\n'.format('')

    #    option +=  '{:11s}'.format('')
    #    for p in self.Options['Particle']:
    #      if cell.MatInfo == 'Graveyard' :
    #        option +=  'imp:{}=0     '.format(p)
    #      else:
    #        option +=  'imp:{}=1.0   '.format(p)

    #    if self.Options['Universe'] is not None:
    #         option += 'U={}  '.format(self.Options['Universe'])
    #    option += '\n'

    #    return option

    def comment_format(self, cComment, mComment=None):
        comment = ""
        if mComment:
            mComment = mComment.split("\n")
            for c in mComment:
                if c:
                    comment += "{:11s}%{}\n".format("", c)

        if cComment.strip() != "":
            cComment = cComment.strip().split("\n")
            for c in cComment:
                if c:
                    comment += "{:11s}%{}\n".format("", c)
        return comment

    def comment_line(self, lineComment):
        lineComment = lineComment.strip().split("\n")
        comment = ""
        if lineComment:
            comment = "% \n"
            for c in lineComment:
                if c:
                    comment += "% {}\n".format(c)
            comment += "% \n"
        return comment

    def __getSurfaceTable__(self):
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

    def __simplifyPlanes__(self, Surfaces):
        offset = len(self.Cells)
        keys = self.surfaceTable.keys()

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

        if opt.prnt3PPlane:
            for p in Surfaces["P"]:
                if p.Surf.pointDef:
                    axis, d = points_to_coeffs(p.Surf.Points)
                    if is_opposite(axis, p.Surf.Axis):
                        self.__change_surf_sign__(p)

        return

    def __sortedSurfaces__(self, Surfaces):
        temp = Surfaces_dict(Surfaces)
        surfList = []
        for ind in range(
            Surfaces.IndexOffset, Surfaces.surfaceNumber + Surfaces.IndexOffset
        ):
            s = temp.getSurface(ind + 1)
            if s is not None:
                surfList.append(s)
                temp.delSurface(ind + 1)
        return surfList

    def __change_surf_sign__(self, p):

        if p.Index not in self.surfaceTable.keys():
            print(
                "{} Surface {} not used in cell definition)".format(p.Type, p.Index),
                p.Surf.Axis,
                p.Surf.Position,
            )
            return

        for ic in self.surfaceTable[p.Index]:
            surf = self.Cells[ic].Definition.get_surfaces_numbers()
            for s in surf:
                if s == p.Index:
                    change_surf_sign(s, self.Cells[ic].Definition)

    def __get_solidCellVolume__(self):

        solidList = []
        volumeList = []
        for m in self.Cells:
            if m.CellType == "solid" and m.__id__ is not None:
                solidList.append(m.label)
                volumeList.append(m.Volume * 1e-3)

        return solidList, volumeList
