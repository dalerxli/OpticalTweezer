set term png 

b = system('cat  p0.1_at0.1_170209_1954_parameters.txt') 

set title b

set output 'p0.1_at0.1_170209_1954_temperature_internal.png'

plot 'p0.1_at0.1_170209_1954_temperature_internal.dat' w l 

unset output 

