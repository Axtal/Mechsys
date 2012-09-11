import optparse
from msys_fig import *

# input
op = optparse.OptionParser()
op.add_option('--tst',  '-t', dest='tst',  default='1', help='test number')
op.add_option('--twos', '-w', dest='twos', default='0', help='two stages')
opts, args = op.parse_args()
two_stages = int(opts.twos)

if opts.tst=='0':
    # 8, 12, 14, 18
    P   = 18
    dat = read_table("owen_hinton_02_P%d.dat"%P)
    plot(dat['r'],dat['st'],'r*',label="OH")

    if two_stages:
        rs1 = read_table("owen_hinton_02_stg1_P%d.res"%P)
        rs2 = read_table("owen_hinton_02_stg2_P%d.res"%P)
        plot (rs1['r'],rs1['st'],'b-',marker='o',lw=2)
        plot (rs2['r'],rs2['st'],'g-',marker='o',lw=2)
    else:
        res = read_table("owen_hinton_02_stg1_P%d.res"%P)
        plot (res['r'],res['st'],'b-',marker='o',lw=2,label="MechSys")

    Gll('r','st')
    show()

elif opts.tst=='1':
    rs1 = read_table("owen_hinton_02_stg1_n41.res")
    dat = read_table("owen_hinton_02_pu.dat")
    da0 = read_table("jiang.dat")
    if two_stages:
        rs2 = read_table("owen_hinton_02_stg2_n41.res")

    subplot(1,2,1)
    if two_stages:
        plot (rs2['ur'],rs2['fr_ext'],'b-',lw=2,marker='.', label='MechSys: Fext (stg2)')
        plot (rs2['ur'],rs2['fr_int'],'c-',lw=2,marker='+', label='MechSys: Fint (stg2)')
    plot   (rs1['ur'],rs1['fr_ext'],'r-',lw=2,marker='o', label='MechSys: Fext (stg1)')
    plot   (rs1['ur'],rs1['fr_int'],'m-',lw=2,marker='*', label='MechSys: Fint (stg1)')
    xlabel ('u')
    ylabel ('force')
    legend (loc='best')

    subplot(1,2,2)
    if two_stages:
        plot (rs2['ur'],rs2['P'],'c-',lw=2, label='MechSys: stg2')
    plot   (rs1['ur'],rs1['P'],'r-',lw=2,   label='MechSys: stg1')
    plot   (da0['u'],da0['p'],'gd',         label='Jiang')
    plot   (dat['u']/100,dat['p'],'bo',     label='OH')
    xlabel ('u')
    ylabel ('P')
    legend (loc='best')
    show   ()

elif opts.tst=='2':
    dat = read_table("owen_hinton_02_mesh.dat")
    for i in range(1,5):
        print (dat['x'][i]-dat['x'][i-1])/10.
