from plotter import *

test = 1

if test==1:
    p = Plotter()
    p.plot ("test.res", fem_res=False, dem_res=True, div_by_p=False)
    p.show ()
