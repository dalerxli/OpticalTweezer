set term png 

b = system('cat  p0.6_at0.1_170208_1340_parameters.txt') 

set title b

set output 'p0.6_at0.1_170208_1340_temperature_internal.png'

plot 'p0.6_at0.1_170208_1340_temperature_internal.dat' w l 

unset output 

