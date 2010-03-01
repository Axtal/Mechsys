from msys_plotter import *
import mechsys as ms

qcamf      = 8.29 #28.15
pcamf      = 5.00 #14.38
p          = Plotter()
p.idx_max  = 58
p.show_k   = True
p.fc_ty    = ['MC', 'VM']
p.fc_p     = pcamf*sqrt(3.0)
p.fc_phi   = ms.M2Phi (qcamf/pcamf, "cam")
p.fc_cu    = ms.qf2cu (qcamf,       "cam")
p.fc_c     = 0.0
p.fc_np    = 40
p.log_p    = False
p.div_by_p = False
#p.plot ("test_triaxialc_walls.res", draw_fl=True)
p.plot ("ttt_c_walls.res", draw_fl=True)
p.show ()
