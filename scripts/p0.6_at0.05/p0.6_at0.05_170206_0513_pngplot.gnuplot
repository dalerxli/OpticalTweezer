set term png 

b = system('cat  p0.6_at0.05_170206_0513_parameters.txt') 

set title b

set output 'p0.6_at0.05_170206_0513_temperature_internal.png'

plot 'p0.6_at0.05_170206_0513_temperature_internal.dat' w l 

unset output 

