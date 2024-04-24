import Part

from .options import Options
from .split_function import split_solid, join_base, SplitBase
from .utils.booleanFunction import BoolSequence


def get_part(slist):
    sol = []
    for s in slist:
        if type(s) is list:
            sol.extend(get_part(s))
        else:
            sol.append(s)
    return sol


def build_solid(cell, boundBox, mode="oneByOne", simplify=False):

    cutCell = cell.makeBox(boundBox)
    # cell.definition = BoolSequence(cell.definition.str)
    cell.clean_undefined()

    celParts = build_depth(cell, SplitBase(cutCell), mode, True, simplify)

    celParts = get_part(celParts)
    # print('celparts',len(celParts))
    shapeParts = []
    for i, s in enumerate(celParts):
        shapeParts.append(s.base)
        # s.base.exportStep('solid{}.stp'.format(i))

    #   tt = fuse_solid(shapeParts)
    #   tt = tt.removeSplitter()
    # tt=Part.makeCompound(shapeParts)
    # tt.exportStep('cell{}.stp'.format(cell.name))
    return shapeParts
    # return fuse_solid(shapeParts)


def build_depth(cell, cutShape, mode, baseBox, simplify=False, loop=0):

    loop += 1
    seq = cell.definition
    if seq.level == 0:
        if baseBox:
            cutShape, cut = build_solid_parts(cell, cutShape, mode)
        else:
            cutShape, cut = build_solid_parts(cell, cutShape, "solids")
        return cutShape

    if type(cutShape) is not list:
        cutShape = [cutShape]
    newCutShape = []
    for i, CS in enumerate(cutShape):
        cbaseBox = baseBox
        # CS.base.exportStep('CS_{}_{}.stp'.format(i,str(cell.definition)))
        # CTable =build_c_table_from_solids(cell.makeBox(CS.base.BoundBox),cell.surfaces,option='full')
        # cell.definition.simplify(CTable)
        cell.definition.group_single()

        if type(cell.definition.elements) is not bool:
            if cell.definition.level == 0:
                tmp = BoolSequence(operator=cell.definition.operator)
                tmp.append(cell.definition)
                cell.definition = tmp

            if seq.operator == "AND":
                part = CS
                for e in cell.definition.elements:
                    part = build_depth(
                        cell.get_sub_cell(e), part, mode, cbaseBox, simplify, loop=loop
                    )
                    cbaseBox = False
                newCutShape.extend(part)
            else:
                cellParts = []
                for e in cell.definition.elements:
                    sub = cell.get_sub_cell(e)
                    part = build_depth(sub, CS, mode, baseBox, simplify, loop=loop)
                    cellParts.extend(part)

                JB = join_base(cellParts)
                if JB.base is not None:
                    newCutShape.append(JB)

        elif cell.definition.elements:
            newCutShape.append(CS)

    cutShape = newCutShape
    return cutShape


def build_solid_parts(cell, base, mode):

    # part if several base in input
    if type(base) is list or type(base) is tuple:
        fullPart = []
        cutPart = []
        for b in base:
            fullList, cutList = build_solid_parts(cell, b, mode)
            fullPart.extend(fullList)
            cutPart.extend(cutList)
        return fullPart, cutPart

    boundBox = base.base.BoundBox
    surfaces = tuple(cell.surfaces.values())
    # print('\nbuild Parts :',mode)
    # print(cell.definition)
    # print(boundBox)

    if mode == "solids":
        # boundBox.enlarge(10)
        if cell.definition.operator == "OR" and False:
            Def = cell.definition
            cell.definition = cell.definition.get_complementary()
            cell.build_shape(boundBox, force=False, simplify=False)
            cell.definition = Def
        else:
            cell.build_shape(boundBox, force=True, simplify=False, fuse=True)

        # print('export')
        # base.base.exportStep('base.stp')
        # name=''.join(str(cell.definition))
        # cell.shape.exportStep('sol{}.stp'.format(name))
    else:
        cell.build_surface_shape(boundBox)

    if not surfaces:
        print("not cutting surfaces")
        return tuple(base.base), tuple()
    if mode == "solids":
        full, cut = split_solid(
            base, surfaces, cell, solidTool=True, tolerance=Options.split_tolerance
        )
    elif mode == "allSurfaces":
        full, cut = split_solid(base, surfaces, cell, tolerance=Options.split_tolerance)

    elif mode == "planeFirst":
        planes = []
        others = []
        for s in surfaces:
            # s.build_shape( boundBox)
            # s.shape.exportStep('Tool_{}_{}.stp'.format(s.type,s.id))
            if s.type == "plane":
                planes.append(s)
            else:
                others.append(s)

        if planes:

            full, cut = split_solid(base, planes, cell, tolerance=Options.split_tolerance)
            # for i,s in enumerate(full):
            #    s.exportStep('fullplane_{}.stp'.format(i))
            # for i,s in enumerate(cut):
            #    s.base.exportStep('cutplane_{}.stp'.format(i))
            # print('planes',full)
            # print('planes',cut)
        else:
            full = []
            cut = base

        if others:
            newf, cut = split_solid(cut, others, cell, tolerance=Options.split_tolerance)
            # print('others',newf)
            # print('others',cut)
        else:
            newf = []
        full.extend(newf)

    elif mode == "otherFirst":
        planes = []
        others = []
        for s in surfaces:
            # s.build_shape( boundBox)
            if s.type == "plane":
                planes.append(s)
            else:
                others.append(s)

        if others:
            full, cut = split_solid(base, others, cell, tolerance=Options.split_tolerance)
            # print('others',full)
            # print('others',cut)
        else:
            full = []
            cut = base

        if planes:
            newf, cut = split_solid(cut, planes, cell, tolerance=Options.split_tolerance)
            # print('planes',newf)
            # print('planes',cut)

        else:
            newf = []

        full.extend(newf)

    elif mode == "oneByOne":
        planes = []
        others = []
        for s in surfaces:
            if s.type == "plane":
                planes.append(s)
            else:
                others.append(s)

        if planes:
            full, cut = split_solid(base, planes, cell, tolerance=Options.split_tolerance)
        else:
            full = []
            cut = base

        # cut[0].base.exportStep('cutPlane.stp')
        for surf in others:
            newf, cut = split_solid(cut, (surf,), cell, tolerance=Options.split_tolerance)
            full.extend(newf)

    elif mode == "otherOneByOne":
        planes = []
        others = []
        for s in surfaces:
            if s.type == "plane":
                planes.append(s)
            else:
                others.append(s)

        cut = base
        full = []
        for surf in others:
            newf, cut = split_solid(cut, (surf,), cell, tolerance=Options.split_tolerance)
            full.extend(newf)

        for surf in planes:
            newf, cut = split_solid(cut, (surf,), cell, tolerance=Options.split_tolerance)
            full.extend(newf)

    return full, cut


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


def noOR(Seq):
    if len(Seq.elements) == 1:
        # Seq.operator = 'AND'
        return Seq
    newOR = BoolSequence(operator="OR")
    neg = []
    for e in Seq.elements:
        AND = BoolSequence(operator="AND")
        AND.append(*neg, e)
        newOR.append(AND)
        neg.append(-e)
    return newOR
