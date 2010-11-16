import optparse
from mechsys         import String, Dict, InpFile, ReadMaterial
from msys_plotter    import *
from msys_invariants import *
from msys_readdata   import *
from msys_fig        import *

# input
op = optparse.OptionParser()
op.add_option('--inp',  '-i', dest='inp',  default='driver.inp',    help='input filename, ex: driver.inp')
op.add_option('--mat',  '-m', dest='mat',  default='materials.inp', help='materials filename, ex: materials.inp')
op.add_option('--tst',  '-t', dest='tst',  default='1',             help='test number')
op.add_option('--fem',  '-f', dest='fem',  default='0')
op.add_option('--both', '-b', dest='both', default='1')
opts, args = op.parse_args()

# input file
inp = InpFile()
inp.Read (opts.inp)

# default initial values
inis = Dict()
inis.Set (-1, {'sx':-inp.pCam0, 'sy':-inp.pCam0, 'sz':-inp.pCam0, 'pw':inp.pw0})

# parse materials file
prms = Dict()
name = String()
ReadMaterial (-1, inp.MatID, opts.mat, name, prms, inis)
print inp
print "Material data:"
print "  name = ", name
print "  prms : ", prms
print "  inis : ", inis

if inp.O2:
    if inp.NDiv==3: out_nods = [48,51,60,68, 0,3,12,15]
    else:           out_nods = [4,5,6,7, 10,11,14,15]
else:
    if inp.NDiv==3: out_nods = [0,3,12,15, 48,51,60,63]
    else:           out_nods = [0,1,2,3, 4,5,6,7]

# res file
fem = 'driver_nod_%d.res'%out_nods[2]
pnt = 'driver.res'
res = fem if opts.fem=='1' else pnt

if opts.tst=='1':
    p = Plotter()
    #p.fc_c    = 0.1
    p.fc_phi   = M_calc_phi(1,'cam')
    p.fc_poct  = 150.0*sqrt(3.0)
    p.show_k   = True
    p.oct_sxyz = True
    p.fsz      = 14

    dat = inp.RefDat.PyStr()
    sim = inp.RefSim.PyStr()
    ana = inp.RefAna.PyStr()

    if not dat=='': p.plot (dat,marker='o',clr='k', label='Dat')
    if not sim=='': p.plot (sim,marker='-',clr='g', label='Sim')
    if not ana=='': p.plot (ana,marker='-',clr='c', label='Ana')

    if opts.both=='1':
        p.lwd=2; p.plot (pnt,clr='blue', markevery=10, draw_ros=True, label='Point (%s)'%name)
        p.lwd=2; p.plot (fem,clr='red',  markevery=10, draw_ros=True, label='FEM (%s)'%name)
    else:
        p.lwd=2; p.plot (res,clr='blue', markevery=10, draw_ros=True, label='%s'%name)
    subplot(2,3,3)
    l = legend(loc='upper left',prop=FontProperties(size=8))
    p.show()

elif opts.tst=='2':
    res = []
    for n in out_nods:
        res.append(read_table('driver_nod_%d.res'%n))

    subplot(2,2,1)
    for k, n in enumerate(out_nods):
        plot (res[k]['Time'], -res[k]['sz'], label='Node # %d'%n)
    xlabel ('Time')
    ylabel (r'$-\sigma_z$')
    grid   ()
    legend (loc='best',prop=FontProperties(size=8))

    subplot(2,2,2)
    for k, n in enumerate(out_nods):
        plot (res[k]['uz'], -res[k]['sz'], label='Node # %d'%n)
    xlabel (r'$u_z$')
    ylabel (r'$-\sigma_z$')
    grid   ()
    legend (loc='best',prop=FontProperties(size=8))

    subplot(2,2,3)
    if res[0].has_key('pw'):
        for k, n in enumerate(out_nods):
            plot (res[k]['Time'], res[k]['pw'], label='Node # %d'%n)
        xlabel ('Time')
        ylabel (r'$p_w$')
        Grid   ()

    subplot(2,2,4)
    if res[0].has_key('pc'):
        for k, n in enumerate(out_nods):
            #plot (log(1.0+res[k]['pc']), res[k]['Sw'], label='Node # %d'%n)
            plot (res[k]['pc'], res[k]['Sw'], label='Node # %d'%n)
        #xlabel (r'$\log{(1+p_c)}$')
        xlabel (r'$p_c$')
        ylabel (r'$S_w$')
        Grid   ()

    show ()
