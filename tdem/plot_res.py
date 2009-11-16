from msys_plotter import *

test = 1

if test==1:
    p = Plotter()
    p.show_k = True
    #p.fc_ty  = 'MC'
    p.fc_p   = -0.3
    p.fc_phi = 25.0
    p.fc_cu  = 1
    p.log_p  = False
    #p.plot ("test_triaxial01a_walls.res", draw_fl=True)
    #p.plot ("test_triaxial01b_walls.res", draw_fl=True)
    p.plot ("test_triaxialc_walls.res", draw_fl=True)
    p.show ()

