#
# Module to load a STEP file
#
import os
import re

import FreeCAD
import Part
from FreeCAD import Import

from ..Utils import Functions as UF
from . import LoadFunctions as LF


# Paco mod
def extractMaterials(filename):
    rho_real = []
    m_dict = {}  # _ Material dictionary
    with open(filename, "rt") as file:
        for line in file:
            vals = line.split()
            if vals[0].startswith("#"):
                continue
            mat_label = int(vals[0])
            rho_real = -float(vals[1])
            matname = " ".join(vals[2:])
            m_dict[mat_label] = (rho_real, matname)
    return m_dict


def LoadCAD(filename, mat_filename, default_mat=[], comp_solids=True):

    # Set document solid tree options when opening CAD differing from version 0.18
    if int(FreeCAD.Version()[1]) > 18:
        LF.set_docOptions()

    cad_simplificado_doc = FreeCAD.newDocument("CAD_simplificado")
    Import.insert(filename, "CAD_simplificado")

    if mat_filename != "":
        if os.path.exists(mat_filename):
            m_dict = extractMaterials(mat_filename)
        else:
            print("Material definition file {} does not exist.".format(mat_filename))
            m_dict = {}
    else:
        m_dict = {}

    s = Part.Shape()
    s.read(filename)
    Solids = s.Solids
    meta_list = []
    for i, s in enumerate(Solids):
        meta_list.append(UF.GEOUNED_Solid(i + 1, s))

    i_solid = 0
    missing_mat = set()

    doc_objects = cad_simplificado_doc.Objects

    for elem in doc_objects:
        if elem.TypeId == "Part::Feature":
            comment = LF.getCommentTree(elem)
            if not elem.Shape.Solids:
                print(
                    "Warning: Element {:} has no associated solid".format(
                        comment + "/" + elem.Label
                    )
                )
                continue
            else:
                tempre_mat = None
                tempre_dil = None

                # MIO: lightly modification of label if required
                label = LF.GetLabel(elem.Label)
                comment = comment + "/" + label
                if elem.InList:
                    # MIO: lightly modification of label if required
                    label_in_list = LF.GetLabel(elem.InList[0].Label)
                    encl_label = re.search(
                        "enclosure(?P<encl>[0-9]+)_(?P<parent>[0-9]+)_", label_in_list
                    )
                    if not encl_label:
                        encl_label = re.search(
                            "enclosure(?P<encl>[0-9]+)_(?P<parent>[0-9]+)_", label
                        )

                    envel_label = re.search(
                        "envelope(?P<env>[0-9]+)_(?P<parent>[0-9]+)_", label_in_list
                    )
                    if not envel_label:
                        envel_label = re.search(
                            "envelope(?P<env>[0-9]+)_(?P<parent>[0-9]+)_", label
                        )

                    # tempre_mat = re.search("(m(?P<mat>\d+)_)",elem.Label)
                    # if not tempre_mat :
                    #    tempre_mat = re.search("(m(?P<mat>\d+)_)",elem.InList[0].Label)

                    # tempre_dil = re.search("(_d(?P<dil>\d*\.\d*)_)",elem.Label)
                    # if not tempre_dil :
                    #    tempre_dil = re.search("(_d(?P<dil>\d*\.\d*)_)",elem.InList[0].Label)

                    # Paco modifications
                    # Search for material definition in tree
                    xelem = [elem]
                    while xelem and not tempre_mat:
                        # MIO: Modification of label if required
                        temp_label = LF.GetLabel(xelem[0].Label)
                        tempre_mat = re.search("_m(?P<mat>\d+)_", "_" + temp_label)
                        xelem = xelem[0].InList

                    # Search for dilution definition in tree
                    xelem = [elem]
                    while xelem and not tempre_dil:
                        # MIO: Modification of label if required
                        temp_label = LF.GetLabel(xelem[0].Label)
                        tempre_dil = re.search("_d(?P<dil>\d*\.\d*)_", temp_label)
                        xelem = xelem[0].InList
                    # Paco end
                else:
                    encl_label = None
                    envel_label = None

                # compSolid Diferent solid of the same cell are stored in the same metaObject (compSolid)
                # enclosures and envelopes are always stored as compound
                if comp_solids or encl_label or envel_label:

                    init = i_solid
                    end = i_solid + len(elem.Shape.Solids)
                    LF.fuseMetaObj(meta_list, init, end)
                    n_solids = 1
                else:
                    n_solids = len(elem.Shape.Solids)

                for i in range(n_solids):
                    meta_list[i_solid].setComments("{}{}".format(comment, (i + 1)))
                    meta_list[i_solid].setCADSolid()

                    if tempre_mat:
                        mat_label = int(tempre_mat.group("mat"))
                        if mat_label in m_dict.keys():
                            meta_list[i_solid].setMaterial(
                                mat_label, m_dict[mat_label][0], m_dict[mat_label][1]
                            )
                        else:
                            if mat_label == 0:
                                meta_list[i_solid].setMaterial(mat_label, 0, 0)
                            else:
                                meta_list[i_solid].setMaterial(
                                    mat_label,
                                    -100,
                                    "Missing material density information",
                                )
                                missing_mat.add(mat_label)
                    else:
                        # print('Warning : No material label associated to solid {}.\nDefault material used instead.'.format(comment))
                        if default_mat:
                            meta_list[i_solid].setMaterial(*default_mat)
                    if tempre_dil:
                        meta_list[i_solid].setDilution(float(tempre_dil.group("dil")))

                    if encl_label is not None:
                        meta_list[i_solid].EnclosureID = int(encl_label.group("encl"))
                        meta_list[i_solid].ParentEnclosureID = int(
                            encl_label.group("parent")
                        )
                        meta_list[i_solid].IsEnclosure = True
                        meta_list[i_solid].CellType = "void"

                    if envel_label is not None:
                        meta_list[i_solid].EnclosureID = int(envel_label.group("env"))
                        meta_list[i_solid].ParentEnclosureID = int(
                            envel_label.group("parent")
                        )
                        meta_list[i_solid].IsEnclosure = True
                        meta_list[i_solid].CellType = "envelope"
                    i_solid += 1

    LF.joinEnvelopes(meta_list)
    if missing_mat:
        print(
            "Warning!! At least one material in the CAD model is not present in the material file"
        )
        print("List of not present materials:", missing_mat)

    enclosure_list = LF.setEnclosureSolidList(meta_list)
    if enclosure_list:
        LF.checkEnclosure(cad_simplificado_doc, enclosure_list)
        # LF.RemoveEnclosure(meta_list)
        return meta_list, enclosure_list
    else:
        return meta_list, []
