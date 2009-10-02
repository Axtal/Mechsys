from numpy import *
from pylab import *
from data_handler import *

#dat = read_table("owen_hinton_nod_2.res")
#uy = array(dat['uy'])
#fy = array(dat['fy'])

#plot(-uy,-fy,'r-',lw=2)
#plot(-uy,-fy,'bo',lw=2)


dat = read_table("owen_hinton.res")
u    = array(dat['u'])
fint = array(dat['fint'])
fext = array(dat['fext'])

plot(-u,-fext,'r-',lw=2)
plot(-u,-fext,'ro')

plot(-u,-fint,'b-',lw=2)
plot(-u,-fint,'b*')

grid()
show()
