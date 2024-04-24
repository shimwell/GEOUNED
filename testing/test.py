import sys
import os
import subprocess
from pathlib import Path

sys.path.append(str(Path(__file__).parents[1] / "src"))
sys.path.append("/usr/lib64/freecad/lib64/")

from geouned import GEOUNED


def set_input(inName, inpDir, outDir):

    if inName.endswith(".step"):
        filename = inName[0:-5]
    elif inName.endswith(".stp"):
        filename = inName[0:-4]
    else:
        filename = inName

    if inpDir == "":
        inpDir = "."
    if outDir == "":
        outDir = "."

    inName = "{}/{}".format(inpDir, inName)
    outName = "{}/{}".format(outDir, filename)

    template = """[Files]
title    = Input Test
step_file = {}
geometry_name = {}

[parameters]
comp_solids = False
vol_card    = False
vol_sdef    = True
void_gen    = True
dummy_mat    = True
min_void_size =  100
cell_summary_file = False
cell_comment_file = False
debug       = False
simplify   = full

[Options]
force_cylinder = False
split_tolerance = 0
new_split_plane = True
n_plane_reverse = 0
""".format(
        inName, outName
    )

    file = open("config.ini", "w")
    file.write(template)
    file.close()


def get_input_list(folder, ext=None):
    filenameList = []
    if ext is None:
        return os.listdir(folder)

    elif type(ext) is str:
        if ext[0] != ".":
            ext = "." + ext
        for f in os.listdir(folder):
            if f.endswith(ext):
                filenameList.append(f)

    elif type(ext) is list or type(ext) is tuple:
        filenameList = []
        extList = []
        for e in ext:
            if e[0] != ".":
                e = "." + e
            extList.append(e)

        for f in os.listdir(folder):
            for e in extList:
                if f.endswith(e):
                    filenameList.append(f)

    return filenameList


def run_mcnp(path, inpFile):
    code = "/opt/mcnp5/last/bin/mcnp5"
    xsdir = "/opt/mcnp.data/ascii/xsdir"
    inp = inpFile
    out = inpFile[0:-1] + "o"
    mctal = inpFile[0:-1] + "m"
    cmd = "cd {} && {} i={} o={} mctal={} xsdir={}".format(
        path, code, inp, out, mctal, xsdir
    )
    os.system(cmd)


def clean(folder):
    filename = "{}/runtpe".format(folder)
    if os.path.isfile(filename):
        os.remove(filename)


def clean_dir(folder):
    for file in get_input_list(folder, (".o", ".m")):
        filename = "{}/{}".format(folder, file)
        os.remove(filename)


def check_lost(outp):
    cmd = "grep 'particles got lost.' {}".format(outp)
    stdout = subprocess.getoutput(cmd)
    lineout = stdout.split("\n")
    lost = 0
    for l in lineout:
        if l == "" or "terminated" in l:
            continue
        lost = int(l[0 : l.index("particles")])
    return lost


def get_mctal_values(mctal):
    file = open(mctal, "r")
    while True:
        line = file.readline()
        if "vals" in line:
            break

    rawData = []
    line = file.readline()
    while "tfc" not in line:
        rawData.extend(line.split())
        line = file.readline()
    file.close()

    values = []
    for i in range(int(len(rawData) / 2)):
        val = float(rawData[2 * i])
        err = float(rawData[2 * i + 1])
        values.append((val, err))
    return values


def print_results(f, res, lost):

    line = "{} :\n".format(f)
    if lost != 0:
        line += "{} particles got lost.\n".format(lost)

    for i, vals in enumerate(res):
        line += "  {:<3d} : {:8.5f}  {:7.5f}\n".format(i + 1, vals[0], vals[1])
    print(line)


def mk_geo_inp(inpDir, outDir):
    for f in get_input_list(inpDir, ("stp", "step")):
        set_input(f, inpDir, outDir)
        GEO = GEOUNED(inifile)
        GEO.SetOptions()
        GEO.Start()
        del GEO


def process_mcnp_folder(outDir):
    clean_dir(outDir)
    for f in get_input_list(outDir, ".mcnp"):
        os.rename(f"{outDir}/{f}", f"{outDir}/{f[:-4]}i")

    for f in get_input_list(outDir, ".i"):
        run_mcnp(outDir, f)
        clean(outDir)


def post_process(folder):
    for fmct in get_input_list(folder, ".m"):
        fout = fmct[:-1] + "o"
        mctal = "{}/{}".format(folder, fmct)
        outp = "{}/{}".format(folder, fout)
        res = get_mctal_values(mctal)
        lost = check_lost(outp)
        print_results(fmct[:-1] + "i", res, lost)


# ****************************************************

misInp = "inputSTEP/Misc"
misOut = "outMCNP/Misc"
cylInp = "inputSTEP/DoubleCylinder"
cylOut = "outMCNP/DoubleCylinder"
torInp = "inputSTEP/Torus"
torOut = "outMCNP/Torus"
lrgInp = "inputSTEP/large"
lrgOut = "outMCNP/large"
ms0Inp = "inputSTEP/"
ms0Out = "outMCNP/"

inifile = "config.ini"
folder = {
    "misc": (misInp, misOut),
    "cyl": (cylInp, cylOut),
    "torus": (torInp, torOut),
    "large": (lrgInp, lrgOut),
    "misc0": (ms0Inp, ms0Out),
}

test = "all"

if test == "all":
    for inpDir, outDir in folder.values():
        mk_geo_inp(inpDir, outDir)

    for inpDir, outDir in folder.values():
        process_mcnp_folder(outDir)

    for inpDir, outDir in folder.values():
        post_process(outDir)

else:
    inpDir, outDir = folder[test]
    mk_geo_inp(inpDir, outDir)
    process_mcnp_folder(outDir)
    post_process(outDir)
