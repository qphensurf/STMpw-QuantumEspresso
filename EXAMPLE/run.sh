#!/bin/bash

echo "this example illustrates how to run STMpw using Quantum Espresso v. 7.0"


cat> scf.in << EOF
&control
 calculation='scf'
 prefix='cu'
 pseudo_dir= '/home/chris/Work/QE/PP/ONCV_TM/ONCVPSP-master/abinit',
 outdir = './tmp',
/
&system
 ibrav= 6
 A    = 7.524323259
 C    = 21.0000
 nat  = 40
 ntyp = 3
 ecutwfc=50
 occupations='smearing'
 degauss=0.005
 nosym = .true.
/
&electrons
conv_thr=1e-6
 electron_maxstep=100
/


ATOMIC_SPECIES
Cu 63.546 Cu_ONCV_PBE_sr.upf
C  12.001  C_ONCV_PBE_sr.upf
H   1.000  H_ONCV_PBE_sr.upf

ATOMIC_POSITIONS (angstrom)
Cu       1.2540539      1.2540539      1.7745096
Cu       1.2540539      3.7621616      1.7745096
Cu       1.2540539      6.2702694      1.7745096
Cu       3.7621616      1.2540539      1.7745096
Cu       3.7621616      3.7621616      1.7745096
Cu       3.7621616      6.2702694      1.7745096
Cu       6.2702694      1.2540539      1.7745096
Cu       6.2702694      3.7621616      1.7745096
Cu       6.2702694      6.2702694      1.7745096
Cu      -0.0000000     -0.0000000      3.5471481
Cu      -0.0000000      2.5081078      3.5471481
Cu      -0.0000000      5.0162155      3.5471481
Cu       2.5081078     -0.0000000      3.5471481
Cu       2.5081078      2.5081078      3.5471481
Cu       2.5081078      5.0162155      3.5471481
Cu       5.0162155     -0.0000000      3.5471481
Cu       5.0162155      2.5081078      3.5471481
Cu       5.0162155      5.0162155      3.5471481
Cu       1.2173089      1.2216649      5.3502926
Cu       1.2643659      3.7354671      5.3721237
Cu       1.3051483      6.2233711      5.3868335
Cu       3.7319576      1.2687723      5.3759110
Cu       3.7592229      3.7620892      5.3179781
Cu       3.7919749      6.2563927      5.3743474
Cu       6.2185227      1.3008503      5.3856822
Cu       6.2531116      3.7887965      5.3721853
Cu       6.3060424      6.3027225      5.3519483
Cu       0.0028142      0.0019800      7.1553930
Cu       7.4086628      2.5105506      7.1394184
Cu       0.1031360      5.0154410      7.1379067
Cu       2.5117288      7.4267928      7.1386841
Cu       2.3653450      2.3837857      7.1780170
Cu       2.6196981      4.9083667      7.2923268
Cu       5.0178023      0.0979564      7.1392638
Cu       4.8908890      2.6154483      7.2963865
Cu       5.1468891      5.1390028      7.1761747
C        3.2668548      3.2827450      8.6647390
C        4.2405386      4.2457130      8.6641832
H        2.8673693      2.8897722      9.6116949
H        4.6374034      4.6444382      9.6099712
K_POINTS (automatic)
 3 3 1 0 0 0 
EOF

 
 mpirun -np 6 pw.x <scf.in >scf.out 
 
 cat >nscf.in <<EOF
 &control
 calculation='nscf'
 prefix='cu'
 pseudo_dir= '/home/chris/Work/QE/PP/ONCV_TM/ONCVPSP-master/abinit',
 outdir = './tmp',
/
&system
 ibrav= 6
 A    = 7.524323259
 C    = 21.0000
 nat  = 40
 ntyp = 3
 ecutwfc=50
 occupations='smearing'
 degauss=0.005
 nbnd = 500
 nosym = .true.
/
&electrons
 conv_thr=1e-8
 electron_maxstep=100
 diago_full_acc= .true.
/


ATOMIC_SPECIES
Cu 63.546 Cu_ONCV_PBE_sr.upf
C  12.001  C_ONCV_PBE_sr.upf
H   1.000  H_ONCV_PBE_sr.upf

ATOMIC_POSITIONS (angstrom)
Cu       1.2540539      1.2540539      1.7745096
Cu       1.2540539      3.7621616      1.7745096
Cu       1.2540539      6.2702694      1.7745096
Cu       3.7621616      1.2540539      1.7745096
Cu       3.7621616      3.7621616      1.7745096
Cu       3.7621616      6.2702694      1.7745096
Cu       6.2702694      1.2540539      1.7745096
Cu       6.2702694      3.7621616      1.7745096
Cu       6.2702694      6.2702694      1.7745096
Cu      -0.0000000     -0.0000000      3.5471481
Cu      -0.0000000      2.5081078      3.5471481
Cu      -0.0000000      5.0162155      3.5471481
Cu       2.5081078     -0.0000000      3.5471481
Cu       2.5081078      2.5081078      3.5471481
Cu       2.5081078      5.0162155      3.5471481
Cu       5.0162155     -0.0000000      3.5471481
Cu       5.0162155      2.5081078      3.5471481
Cu       5.0162155      5.0162155      3.5471481
Cu       1.2173089      1.2216649      5.3502926
Cu       1.2643659      3.7354671      5.3721237
Cu       1.3051483      6.2233711      5.3868335
Cu       3.7319576      1.2687723      5.3759110
Cu       3.7592229      3.7620892      5.3179781
Cu       3.7919749      6.2563927      5.3743474
Cu       6.2185227      1.3008503      5.3856822
Cu       6.2531116      3.7887965      5.3721853
Cu       6.3060424      6.3027225      5.3519483
Cu       0.0028142      0.0019800      7.1553930
Cu       7.4086628      2.5105506      7.1394184
Cu       0.1031360      5.0154410      7.1379067
Cu       2.5117288      7.4267928      7.1386841
Cu       2.3653450      2.3837857      7.1780170
Cu       2.6196981      4.9083667      7.2923268
Cu       5.0178023      0.0979564      7.1392638
Cu       4.8908890      2.6154483      7.2963865
Cu       5.1468891      5.1390028      7.1761747
C        3.2668548      3.2827450      8.6647390
C        4.2405386      4.2457130      8.6641832
H        2.8673693      2.8897722      9.6116949
H        4.6374034      4.6444382      9.6099712
K_POINTS (automatic)
 3 3 1 0 0 0 
EOF

 
 mpirun -np 6 pw.x <nscf.in >nscf.out
 
cat >pp.in <<EOF
 &inputpp
   prefix="cu"
   outdir='./tmp'
/
EOF

cat >input.STMpw <<EOF
4.30 ! phi in eV 
2 ! number of voltages to calculate
-1.200 1.200  ! tip-substrate voltage in V (tip to mass)
50 ! sampling points in z (perpendicular to the surface)
10  ! Zmax, maximum tip-surface distance (from 'Zsurf') in \AA
.342 ! Zsurf, origin of the surface in direct coordinates
.550 ! z_s, direct coordinates. Zsurf and z_s must have the same origin!
F ! whether we use Bardeen aprox. or just Tersoff-Hamman
F ! whether we calculate the dIdV curve
POSCAR ! POSCAR file
WAVECAR  ! wf file 
F ! Do we read a mapfile?
OUTCAR ! OUTCAR file
T ! .true. or .false. whether the k-point sampling contains the Gamma point or not.
T ! wsxm output?
1000 ! for WSxM: increase the isovalues
T ! gnuplot output?
F ! cube output?
EOF

mpirun -np 1 STMpw.out <pp.in >pp.out

## The following creates STM differential conductance (dI/dV) images; make sure you have Cond_gnu.out on the path
## note that "gnu" does not refer to anything related to gnuplot specifically but the output is conveniently formatted for gnuplot
## in this example I use python 3.8 to create the plots so make sure that you have a compatible python 3 installation with matplotlib installed

# Variables:
   voltages="-1.200 1.200"
#  heights to use
   heights="4 5 6 7 8"
#  Repetitions in x and y
   nx=2
   ny=2

for j in $voltages; do
cd V_${j}
for i in $heights
do
cat > cond.dat << !
1
$i
$nx
$ny
0
!
Cond_gnu.out < cond.dat
cat > stm.py << EOF
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.ndimage as ndimage
from matplotlib import cm

data=np.loadtxt('Conductance.gnu')

fs=16

X=data[:,0]
Y=data[:,1]
xdim = np.count_nonzero(X == 0)
ydim = np.count_nonzero(Y== 0)
Z=data[:,2].reshape(ydim,xdim)

Z2 = ndimage.gaussian_filter(Z, sigma=3.0, order=0)
Z2=Z2.T

ratio=xdim/ydim
fig, axs = plt.subplots(1,1, figsize=(8,8*ratio))

IM=axs.imshow(Z2, extent=[0,np.max(X)*0.529, 0, np.max(Y)*0.529],cmap='afmhot')

divider = make_axes_locatable(axs)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(IM, cax=cax)
cb.ax.tick_params(labelsize=fs) 

axs.tick_params(axis='both', labelsize = fs)
axs.set_xlabel(r'x ($\AA$)', size=fs)
axs.set_ylabel(r'y ($\AA$)', size=fs)
axs.text(0.90, np.max(Y)*0.529*0.98, '%0.2f V' %$j, color='white', size=fs*2, ha='left', va='top')

fig.tight_layout()
fig.savefig("Conductance-${i}A.png", dpi=400)
EOF

python stm.py

done

rm -f stm.gnu Conductance.gnu cond.dat


cd ..
done


## The following creates STM Topography images; make sure you have Imagen_Bardeen_gnu.out on the path
## in this example I use python 3.8 to create the plots so make sure that you have a compatible python 3 installation with matplotlib installed

# Variables:
   voltages="-1.200 1.200"
#  heights to use
   heights="4 5 6 7 8"
#  Repetitions in x and y
   nx=2
   ny=2

for j in $voltages; do
cd V_${j}
for i in $heights
do
cat > topo.dat << !
F
$nx
$ny
1
$i
0
!
Imagen_Bardeen_gnu.out < topo.dat
cat > stm.py << !
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.ndimage as ndimage
from matplotlib import cm

data=np.loadtxt('Topography.gnu')

fs=16

X=data[:,0]
Y=data[:,1]
xdim = np.count_nonzero(X == 0)
ydim = np.count_nonzero(Y== 0)
Z=data[:,2].reshape(ydim,xdim)

Z2 = ndimage.gaussian_filter(Z, sigma=1.0, order=0)
Z2=Z2.T

ratio=xdim/ydim
fig, axs = plt.subplots(1,1, figsize=(8,8*ratio))

IM=axs.imshow(Z2, extent=[0,np.max(X)*0.529, 0, np.max(Y)*0.529],cmap='afmhot')

divider = make_axes_locatable(axs)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(IM, cax=cax)
cb.ax.tick_params(labelsize=fs) 

axs.tick_params(axis='both', labelsize = fs)
axs.set_xlabel(r'x ($\AA$)', size=fs)
axs.set_ylabel(r'y ($\AA$)', size=fs)
axs.text(0.90, np.max(Y)*0.529*0.98, '%0.2f V' %$j, color='white', size=fs*2, ha='left', va='top')

fig.tight_layout()
fig.savefig("Topo-${i}A.png", dpi=400)

# plt.close(fig)
!
python stm.py
done

rm -f stm.gnu Topography.gnu cond.dat


cd ..
done

