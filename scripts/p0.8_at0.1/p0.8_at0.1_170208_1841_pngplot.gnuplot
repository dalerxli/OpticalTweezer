set term png 

b = system('cat  p0.8_at0.1_170208_1841_parameters.txt') 

set title b

set output 'p0.8_at0.1_170208_1841_temperature_internal.png'

plot 'p0.8_at0.1_170208_1841_temperature_internal.dat' w l 

unset output 

