# 16k systems and Lz/Lx = 5
# averaged mole fraction profiles

file1='rho0.99/a.tmass.align.avg'
file2='rho0.85/a.tmass.align.avg'
file3='rho0.81/a.tmass.align.avg'
file4='rho0.78/a.tmass.align.avg'

nbin1=200
nbin2=200
nbin3=200
nbin4=216

Rzx=5
Lx1=14.85407461
Lx2=15.65173005
Lx3=15.92995119
Lx4=16.11427993

# labeling for figures
#set key at screen 0.0,0.9 title "(c)"
unset key

# output terminal
# load '~/Utility/gnuplot/terminal/svg.plt # svg default file
load '~/Utility/gnuplot/terminal/eps.plt' # postscript

# Palette
#load '~/Utility/gnuplot/palette/moreland.pal' # for pm3d
load '~/Utility/gnuplot/palette/set1.pal' # for line arts

# axes
#load '~/Utility/gnuplot/3dmap/axes.plt' # for pm3d
load '~/Utility/gnuplot/2d/axes-x1y1.plt' # for x1y1
set xl "z / {L}_{z}"
set yl "{/Symbol \162}(z) / {/Symbol \163}^{-3}"
set xr[0:1]
set xtics 0.2
set mxtics 2
scale=1
set yr[0.70:1.05]
set ytics 0.70, 0.10
set mytics 2
#set cbr[0:1]
#set cbtics 0.5
#set cbl "{x}_{A}^{align}(z)" offset 0.5

# smooth splines
#set samples 1000

# gaussian fitting
PI=3.14159
# rho0.9998
gau1(x) = a1/(2*PI*s1**2)**0.5*exp(-(x-m1)**2/(2*s1**2)) +  a1/(2*PI*s1**2)**0.5*exp(-(x-1+m1)**2/(2*s1**2))+b1
a1=-0.02;s1=-0.02;m1=0.25;b1=1.02
fit gau1(x) file1 u ($0/nbin1):(nbin1*scale*$1/Lx1/Lx1/Lx1/Rzx) via a1, s1, m1, b1
fit gau1(x) file1 u ($0/nbin1):(nbin1*scale*$1/Lx1/Lx1/Lx1/Rzx):(nbin1*scale*$2/Lx1/Lx1/Lx1/Rzx) via a1, s1, m1, b1
# rho0.8546
gau2(x) = a2/(2*PI*s2**2)**0.5*exp(-(x-m2)**2/(2*s2**2)) +  a2/(2*PI*s2**2)**0.5*exp(-(x-1+m2)**2/(2*s2**2))+b2
a2=-0.02;s2=-0.02;m2=0.25;b2=0.87
fit gau2(x) file2 u ($0/nbin2):(nbin2*scale*$1/Lx2/Lx2/Lx2/Rzx) via a2, s2, m2, b2
fit gau2(x) file2 u ($0/nbin2):(nbin2*scale*$1/Lx2/Lx2/Lx2/Rzx):(nbin2*scale*$2/Lx2/Lx2/Lx2/Rzx) via a2, s2, m2, b2
# rho0.8106
gau3(x) = a3/(2*PI*s3**2)**0.5*exp(-(x-m3)**2/(2*s3**2)) +  a3/(2*PI*s3**2)**0.5*exp(-(x-1+m3)**2/(2*s3**2))+b3
a3=-0.02;s3=-0.02;m3=0.25;b3=0.82
fit gau3(x) file3 u ($0/nbin3):(nbin3*scale*$1/Lx3/Lx3/Lx3/Rzx) via a3, s3, m3, b3
fit gau3(x) file3 u ($0/nbin3):(nbin3*scale*$1/Lx3/Lx3/Lx3/Rzx):(nbin3*scale*$2/Lx3/Lx3/Lx3/Rzx) via a3, s3, m3, b3
# rho0.7831
gau4(x) = a4/(2*PI*s4**2)**0.5*exp(-(x-m4)**2/(2*s4**2)) +  a4/(2*PI*s4**2)**0.5*exp(-(x-1+m4)**2/(2*s4**2))+b4
a4=-0.02;s4=-0.02;m4=0.25;b4=0.74
fit gau4(x) file4 u ($0/nbin4):(nbin4*scale*$1/Lx4/Lx4/Lx4/Rzx) via a4, s4, m4, b4
fit gau4(x) file4 u ($0/nbin4):(nbin4*scale*$1/Lx4/Lx4/Lx4/Rzx):(nbin4*scale*$2/Lx4/Lx4/Lx4/Rzx) via a4, s4, m4, b4

# fitting
## rho0.9998
#f1(x) = 1-wr1+0.5*(2*wr1-1)*(tanh((x-b1+c1)/lamda1)-tanh((x-b1-c1)/lamda1))
#wr1=0.959;b1=0.5;c1=0.25;lamda1=0.1
#fit f1(x) file1 u ($0/nbin1):1 via wr1, b1, c1, lamda1
#fit f1(x) file1 u ($0/nbin1):1:2 via wr1, b1, c1, lamda1
## rho0.8546
#f2(x) = 1-wr2+0.5*(2*wr2-1)*(tanh((x-b2+c2)/lamda2)-tanh((x-b2-c2)/lamda2))
#wr2=0.959;b2=0.5;c2=0.25;lamda2=0.1
#fit f2(x) file2 u ($0/nbin2):1 via wr2, b2, c2, lamda2
#fit f2(x) file2 u ($0/nbin2):1:2 via wr2, b2, c2, lamda2
## rho0.8106
#f3(x) = 1-wr3+0.5*(2*wr3-1)*(tanh((x-b3+c3)/lamda3)-tanh((x-b3-c3)/lamda3))
#wr3=0.959;b3=0.5;c3=0.25;lamda3=0.1
#fit f3(x) file3 u ($0/nbin3):1 via wr3, b3, c3, lamda3
#fit f3(x) file3 u ($0/nbin3):1:2 via wr3, b3, c3, lamda3
## rho0.7831
#f4(x) = 1-wr4+0.5*(2*wr4-1)*(tanh((x-b4+c4)/lamda4)-tanh((x-b4-c4)/lamda4))
#wr4=0.959;b4=0.5;c4=0.25;lamda4=0.1
#fit f4(x) file4 u ($0/nbin4):1 via wr4, b4, c4, lamda4
#fit f4(x) file4 u ($0/nbin4):1:2 via wr4, b4, c4, lamda4

# plot for fitting
#plot file1 u ($0/nbin1):1:2 w ye ls 91 lw 1, f1(x) ls 1 lw 4, \
#file2 u ($0/nbin2):1:2 w ye ls 92 lw 1, f2(x) ls 2 lw 4, \
#file3 u ($0/nbin3):1:2 w ye ls 93 lw 1, f3(x) ls 3 lw 4, \
#file4 u ($0/nbin4):1:2 w ye ls 94 lw 1, f4(x) ls 4 lw 4

# plot for filled circles
#plot 'hjung.out' i 0 u 3:1:(0.0065) w circles linecolor rgb '#E41A1C' fill solid 0.9 border lt -1 title "This Work"

# plot for tmass
plot file1 u ($0/nbin1):(nbin1*scale*$1/Lx1/Lx1/Lx1/Rzx):(nbin1*scale*$2/Lx1/Lx1/Lx1/Rzx) w ye ls 91 lw 1, \
file2 u ($0/nbin2):(nbin2*scale*$1/Lx2/Lx2/Lx2/Rzx):(nbin2*scale*$2/Lx2/Lx2/Lx2/Rzx) w ye ls 92 lw 1,  \
file3 u ($0/nbin3):(nbin3*scale*$1/Lx3/Lx3/Lx3/Rzx):(nbin3*scale*$2/Lx3/Lx3/Lx3/Rzx) w ye ls 93 lw 1,  \
file4 u ($0/nbin4):(nbin4*scale*$1/Lx4/Lx4/Lx4/Rzx):(nbin4*scale*$2/Lx4/Lx4/Lx4/Rzx) w ye ls 94 lw 1, \
gau1(x) ls 1 lw 4, gau2(x) ls 2 lw 4,  gau3(x) ls 3 lw 4,  gau4(x) ls 4 lw 4
