set term png 

b = system('cat  p0.1_at0.1_170209_1955_parameters.txt') 

set title b

set output 'p0.1_at0.1_170209_1955_temperature_internal.png'

plot 'p0.1_at0.1_170209_1955_temperature_internal.dat' w l 

unset output 

