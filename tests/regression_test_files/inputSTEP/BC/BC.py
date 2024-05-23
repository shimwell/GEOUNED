# openMC geometry script generated by GEOUNED
import openmc

###############################################################################
# Define problem geometry
###############################################################################


# Surface setup
S1 = openmc.Plane(a=0.7660444332791261,b=0.6427876214132181,c=1.9624920150279804e-14,d=-172.499860606956)
S2 = openmc.Plane(a=-0.5339138119292677,b=0.8303835411774745,c=0.15937131477375774,d=1469.6767710058505)
S3 = openmc.Plane(a=-0.04132029623640592,b=-0.9991459518603555,c=2.619034079809351e-15,d=-1076.953201832466)
S4 = openmc.Plane(a=-0.7497144041038003,b=0.6581628837016963,c=0.0689197417054896,d=1574.4284838080853)
S5 = openmc.Plane(a=-0.7941079324212609,b=0.3574324877858954,c=0.4915634326725455,d=1056.0673152355382)
S6 = openmc.Plane(a=-0.23108078046756733,b=0.7531775709073565,c=0.6158938379141279,d=774.2707076603047)
S7 = openmc.Plane(a=-0.7602144209302707,b=0.6496722513773027,c=0.0,d=1614.3079387793098)
S8 = openmc.Plane(a=-0.6122159223828901,b=0.7648615982474182,c=0.20044550358457108,d=1460.6708255532078)
S9 = openmc.Plane(a=-0.5339138119292673,b=0.8303835411774748,c=0.1593713147737576,d=1466.7388593706023)
S10 = openmc.Plane(a=-0.09948299242783774,b=-0.9372011204276478,c=-0.3343010530745669,d=-751.9981927955782)
S11 = openmc.Plane(a=0.5871010243189861,b=-0.5216095555996603,c=-0.6190604645353491,d=-929.1090584764745)
S12 = openmc.Plane(a=0.09948299242783747,b=0.9372011204276478,c=0.3343010530745669,d=753.9601665908509)
S13 = openmc.Plane(a=0.04132029623640854,b=0.9991459518603554,c=0.0,d=1074.7367117739009)
S14 = openmc.Plane(a=-0.9355273727615511,b=0.24455850778955301,c=0.25491110427287605,d=1224.40513475057)
S15 = openmc.Plane(a=-0.5871010243189837,b=0.5216095555996624,c=0.6190604645353498,d=926.7781863743836)
S16 = openmc.Plane(a=-0.6774011929718843,b=0.7154768836078964,c=0.17093990986016308,d=1496.769904745614)
S17 = openmc.Plane(a=-0.79410793242126,b=0.3574324877858938,c=0.49156343267254793,d=1053.7145594858289)
S18 = openmc.Plane(a=0.7178273776416673,b=-0.6844730256836112,c=-0.12735985638971906,d=-1532.4305249637562)
S19 = openmc.Plane(a=-0.9355273727615511,b=0.24455850778955301,c=0.25491110427287605,d=1222.0390704559463)
S20 = openmc.Plane(a=0.7602144209302698,b=-0.6496722513773039,c=-1.3197854849262398e-15,d=-1611.3489360862638)
S21 = openmc.Plane(a=-0.7178273776416672,b=0.6844730256836112,c=0.12735985638971908,d=1535.3877668749064)
S22 = openmc.Plane(a=-0.612215922382885,b=0.7648615982474217,c=0.20044550358457303,d=1463.6091539115253)
S23 = openmc.Plane(a=-0.6774011929718843,b=0.7154768836078964,c=0.17093990986016308,d=1499.7204413018176)
S24 = openmc.Plane(a=-0.23108078046756625,b=0.7531775709073572,c=0.6158938379141273,d=776.4472177123853)
S25 = openmc.Plane(a=-0.7497144041038003,b=0.6581628837016963,c=0.0689197417054896,d=1571.4702163526724)
S26 = openmc.Plane(a=0.4582650131091048,b=0.3845300154357572,c=0.8013300474767486,d=-543.3970773072595)
S27 = openmc.Cylinder(x0=-1148.0653755431501,y0=1099.84885492107,z0=-553.0,r=3.0,dx=-0.6427876214132191,dy=0.7660444332791252,dz=1.5226600442908503e-15)
S28 = openmc.Cylinder(x0=-1151.48896314332,y0=1103.92892761604,z0=-553.0,r=14.0,dx=-0.642787621413092,dy=0.766044433279232,dz=-3.04019139794915e-13)
S29 = openmc.Plane(a=0.3810598832843514,b=0.3197472174608166,c=-0.8674993269607659,d=396.5215787711977)
S30 = openmc.Plane(a=0.7393529106114152,b=0.6203907738909,c=-0.2616726222620143,d=-18.88896452772267)
S31 = openmc.Plane(a=0.7538154089568867,b=0.6325262512957599,c=0.17796873500953647,d=-265.21069941271315)
S32 = openmc.Plane(a=0.6582946325044312,b=0.5523748005129488,c=0.5114002899628974,d=-428.0887601805737)
S33 = openmc.Plane(a=0.16251502549115546,b=0.13636630218981532,c=0.9772374829675153,d=-574.0761680903925)
S34 = openmc.Plane(a=0.6443417442618621,b=0.5406669367707714,c=0.5408354464016684,d=-441.47636864660655)
S35 = openmc.Plane(a=0.5913521440836286,b=0.4962033814226792,c=0.6356766835053764,d=-481.6924917170763)
S36 = openmc.Plane(a=0.7265201771115892,b=0.6096228316094395,c=0.3170621948927908,d=-335.95014087722404)
S37 = openmc.Plane(a=0.4159337689120087,b=0.34900988293843116,c=0.8397566322988788,d=-555.5271862645349)
S38 = openmc.Plane(a=0.6920555976811689,b=0.5807036148216403,c=-0.428768424037079,d=82.55637992581215)
S39 = openmc.Plane(a=0.7162365824226509,b=0.6009938708305133,c=0.3546992038667094,d=-355.29167928263496)
S40 = openmc.Plane(a=0.5262660190183471,b=0.4415896361877841,c=0.7266654460190886,d=-518.1719343630841)
S41 = openmc.XPlane(x0=-1171.2431706642155)
S42 = openmc.XPlane(x0=-1152.0021214666456)
S43 = openmc.YPlane(y0=1115.5849830427574)
S44 = openmc.YPlane(y0=1135.3998928624294)
S45 = openmc.ZPlane(z0=-568.000000000009)
S46 = openmc.ZPlane(z0=-538.0000000000058)
S47 = openmc.Sphere(x0=-1161.6226460654307,y0=1125.4924379525933,z0=-553.0000000000075,r=20.79680089491407, boundary_type="vacuum")

# Cell definition 
C1 = openmc.Cell(name="", region=((+S1 & +S2 & +S3 & +S27 & +S29) | (+S5 & +S2 & -S28 & -S4) | (-S28 & +S8 & +S9 & -S7 & -S6) | (+S8 & -S10 & -S9 & -S6 & +S27 & +S30) | (+S5 & +S7 & +S9 & -S14 & -S13 & -S4) | (-S28 & +S6 & +S9 & +S16 & -S15 & -S7) | (+S6 & +S10 & +S16 & -S15 & +S27 & +S31) | (+S32 & -S18 & +S15 & -S17 & -S6 & +S27) | (+S1 & +S33 & -S20 & -S17 & +S19 & +S27) | (+S2 & -S28 & -S21 & -S11 & +S22 & -S5) | (+S2 & -S28 & -S23 & +S11 & +S22) | (+S2 & -S28 & -S23 & -S22 & +S24) | (-S24 & -S22 & +S2 & -S28) | (-S28 & +S7 & -S23 & +S9 & +S11 & +S13 & +S24 & -S2) | (-S28 & +S5 & +S7 & +S13 & -S14 & -S4 & -S2) | (-S28 & +S7 & +S9 & -S21 & +S13 & -S11 & -S5 & -S2) | (+S34 & +S7 & +S27 & +S9 & -S21 & -S13 & -S5) | (+S35 & +S11 & -S13 & +S21 & +S27) | (+S1 & -S28 & +S9 & +S17 & -S7 & +S25) | (-S28 & +S9 & -S18 & +S15 & -S17 & -S7) | (+S1 & -S28 & -S20 & -S16 & +S19 & -S10) | (+S36 & +S6 & +S10 & -S18 & +S15 & -S17 & +S27) | (+S26 & -S19 & +S17 & +S25 & -S6 & +S27) | (+S1 & +S37 & -S20 & +S17 & +S19 & -S6 & +S27) | (+S1 & -S28 & +S38 & +S8 & +S9 & -S24 & +S11 & +S13 & -S12 & +S27 & -S2) | (-S28 & +S7 & +S8 & -S24 & -S22 & +S9 & +S12 & -S2) | (-S28 & +S6 & +S39 & +S27 & -S18 & +S16 & -S10 & -S9 & +S25 & -S5) | (+S16 & +S18 & +S19 & -S20) | (-S19 & +S16 & +S17 & +S18 & -S9 & +S25) | (-S18 & -S17 & -S16 & +S15 & -S10 & -S9) | (-S19 & -S16 & +S17 & -S10 & +S25) | (+S6 & +S10 & -S19 & +S17 & +S25) | (+S1 & +S6 & -S5 & +S40 & +S10 & -S20 & +S15 & +S17 & -S14 & +S19 & +S27) | (-S7 & -S5 & +S1 & +S14) | (+S1 & +S5 & +S10 & -S20 & -S14) | (+S1 & +S5 & +S10 & -S20 & +S14 & -S7) | (-S28 & +S5 & -S20 & -S18 & +S16 & -S14 & +S19 & -S10 & -S9 & -S4) | (+S1 & -S28 & +S14 & +S16 & -S10 & -S9 & -S7) | (-S28 & +S5 & -S19 & -S18 & +S16 & -S9 & +S25)))
C2 = openmc.Cell(name="Automatic Generated Void Cell. Enclosure(-1171.243, -1152.002, 1115.585, 1135.400, -568.000, -538.000). Enclosed cells : (1)", region=(+S41 & -S42 & +S43 & -S44 & +S45 & -S46 & ((+S28 | -S8 | -S9 | +S7 | +S6) & (+S28 | -S6 | -S9 | -S16 | +S15 | +S7) & (-S2 | +S28 | +S21 | +S11 | -S22 | +S5) & (+S24 | +S22 | -S2 | +S28) & (+S28 | -S7 | -S9 | +S21 | -S13 | +S11 | +S5 | +S2) & (+S28 | -S9 | +S18 | -S15 | +S17 | +S7) & (+S28 | -S7 | -S8 | +S24 | +S22 | -S9 | -S12 | +S2) & (-S16 | -S18 | -S19 | +S20) & (+S19 | -S16 | -S17 | -S18 | +S9 | -S25) & (+S18 | +S17 | +S16 | -S15 | +S10 | +S9) & (+S19 | +S16 | -S17 | +S10 | -S25) & (-S6 | -S10 | +S19 | -S17 | -S25) & (+S28 | -S5 | +S19 | +S18 | -S16 | +S9 | -S25) & ((-S1) | ((-S5 | -S10 | +S20 | +S14) & ((+S7) | (((-S14) | (((+S5) & (-S10 | +S20)) & (-S16 | +S28 | +S9 | +S10))) & (-S25 | -S17 | +S28 | -S9))) & ((((-S27) | ((+S17 | -S33) & (-S37 | -S17 | +S6) & (-S40 | -S15 | +S14 | +S5 | -S10 | -S6 | -S17))) & (+S28 | +S16 | +S10)) | (+S20 | -S19)) & (-S29 | -S3 | -S2 | -S27) & (-S38 | +S24 | -S13 | +S12 | -S11 | -S8 | +S28 | -S9 | -S27 | +S2))) & (((-S2 | +S28) & ((+S14) | (((-S7) | ((+S13 | -S9) & (+S2 | +S28 | -S13))) & (+S20 | -S19 | +S18 | -S16 | +S10 | +S28 | +S9)))) | (-S5 | +S4)) & (((-S2 | +S22 | -S24) & ((-S11) | ((-S2 | -S22) & (-S13 | -S9 | -S7 | -S24 | +S2)))) | (+S28 | +S23)) & (-S26 | +S19 | -S17 | -S25 | +S6 | -S27) & (-S30 | -S8 | +S10 | +S9 | +S6 | -S27) & (-S31 | -S6 | -S10 | -S16 | +S15 | -S27) & (-S32 | +S18 | -S15 | +S17 | +S6 | -S27) & (-S34 | -S7 | -S27 | -S9 | +S21 | +S13 | +S5) & (-S35 | -S11 | +S13 | -S21 | -S27) & (-S36 | -S6 | -S10 | +S18 | -S15 | +S17 | -S27) & (-S39 | +S28 | -S6 | -S27 | +S18 | -S16 | +S10 | +S9 | -S25 | +S5))))
C3 = openmc.Cell(name="", region=(-S47 & (-S41 | +S42 | -S43 | +S44 | -S45 | +S46)))

geometry = openmc.Geometry([C1, C2, C3])
geometry.export_to_xml()
