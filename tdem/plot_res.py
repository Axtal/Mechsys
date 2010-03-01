from msys_plotter import *
import mechsys as ms

voronoi = 1 # 0 => spheres

if voronoi: # voronoi
    qcamf      = 9.35
    pcamf      = 5.0
    p          = Plotter()
    p.show_k   = False
    p.fc_p     = pcamf*sqrt(3.0)
    p.fc_phi   = ms.M2Phi (qcamf/pcamf, "cam")
    p.fc_cu    = ms.qf2cu (qcamf,       "cam")
    p.fc_c     = 0.0
    p.fc_np    = 40
    p.log_p    = False
    p.div_by_p = False
    p.fc_ty    = ['MN', 'LD', 'MC', 'VM']
    p.pcte     = True
    p.justone  = -1
    p.fc_np    = 20
    p.maxed    = 10.0
    #p.maxev    = 4.0
    #p.maxidx    = 10
    p.plot ("ttt_c_walls.res",   clr='blue',   txtmax=True,  txtlst=True,  draw_fl=True,  draw_ros=True)
    p.show ()

else:
    qcamf      = 8.31
    pcamf      = 5.0
    p          = Plotter()
    p.idx_max  = -1
    p.show_k   = False
    p.fc_p     = pcamf*sqrt(3.0)
    p.fc_phi   = ms.M2Phi (qcamf/pcamf, "cam")
    p.fc_cu    = ms.qf2cu (qcamf,       "cam")
    p.fc_c     = 0.0
    p.fc_np    = 40
    p.log_p    = False
    p.div_by_p = False
    p.fc_ty    = ['MN', 'LD', 'MC', 'VM']
    p.pcte     = True
    p.justone  = -1
    p.fc_np    = 20
    p.plot ("ttt_c_walls.res",   clr='blue',   txtmax=True,  txtlst=True,  draw_fl=True,  draw_ros=True)
    p.show ()
