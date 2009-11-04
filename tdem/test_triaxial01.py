from mechsys import *

d = DEM_Domain()
d.AddVoronoiPacking(-1, 0.1, 6, 6, 6, 6, 6, 6, True, 1.0)
d.WritePOV("test", [0,35,3])
