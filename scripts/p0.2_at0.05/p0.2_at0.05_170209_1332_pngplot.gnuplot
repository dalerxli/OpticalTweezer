set term png 

b = system('cat  p0.2_at0.05_170209_1332_parameters.txt') 

set title b

set output 'p0.2_at0.05_170209_1332_temperature_internal.png'

plot 'p0.2_at0.05_170209_1332_temperature_internal.dat' w l 

unset output 

