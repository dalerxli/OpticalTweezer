set term png 

b = system('cat  p0.9_at0.1_170209_1846_parameters.txt') 

set title b

set output 'p0.9_at0.1_170209_1846_temperature_internal.png'

plot 'p0.9_at0.1_170209_1846_temperature_internal.dat' w l 

unset output 

