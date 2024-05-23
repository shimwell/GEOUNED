# openMC geometry script generated by GEOUNED
import openmc

###############################################################################
# Define problem geometry
###############################################################################


# Surface setup
S1 = openmc.XPlane(x0=0.0)
S2 = openmc.XPlane(x0=1.0)
S3 = openmc.YPlane(y0=1.0)
S4 = openmc.YPlane(y0=0.0)
S5 = openmc.ZPlane(z0=1.0)
S6 = openmc.ZPlane(z0=0.0)
S7 = openmc.Plane(a=0.7071067811865477,b=0.7071067811865474,c=-0.0,d=0.0)
S8 = openmc.ZCylinder(x0=0.2,y0=0.2,r=0.2)
S9 = openmc.ZCylinder(x0=0.13999999999999999,y0=0.13999999999999999,r=0.09000000000000001)
S10 = openmc.Plane(a=0.7071067811865477,b=0.7071067811865475,c=-0.0,d=0.14142135623730953)
S11 = openmc.XPlane(x0=-1.0)
S12 = openmc.XPlane(x0=2.0)
S13 = openmc.YPlane(y0=-1.0)
S14 = openmc.YPlane(y0=2.0)
S15 = openmc.ZPlane(z0=-1.0)
S16 = openmc.ZPlane(z0=2.0)
S17 = openmc.Sphere(x0=0.5,y0=0.5,z0=0.5,r=2.6500377355803826, boundary_type="vacuum")

# Cell definition 
C1 = openmc.Cell(name="", region=((+S1 & +S4 & +S6 & +S8 & +S10 & -S5 & -S3 & -S2) | (-S8 & +S9 & -S5 & +S6)))
C2 = openmc.Cell(name="Automatic Generated Void Cell. Enclosure(-1.000, 2.000, -1.000, 2.000, -1.000, 2.000). Enclosed cells : (1)", region=(+S11 & -S12 & +S13 & -S14 & +S15 & -S16 & (((+S8 | -S9) & (-S10 | -S4 | +S3 | +S2 | -S1 | -S8)) | (+S5 | -S6))))
C3 = openmc.Cell(name="", region=(-S17 & (-S11 | +S12 | -S13 | +S14 | -S15 | +S16)))

geometry = openmc.Geometry([C1, C2, C3])
geometry.export_to_xml()
