set term png 

b = system('cat  p0.3_at0.05_170209_1620_parameters.txt') 

set title b

set output 'p0.3_at0.05_170209_1620_temperature_internal.png'

plot 'p0.3_at0.05_170209_1620_temperature_internal.dat' w l 

unset output 

