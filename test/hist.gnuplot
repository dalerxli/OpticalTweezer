set multiplot layout 1,3
plot "comhistxy.dat" using 1:3:5 with image
plot "comhistxz.dat" using 1:3:5 with image
plot "comhistyz.dat" using 1:3:5 with image
