set term png 

b = system('cat  p0.5_at0.1_170212_1717_parameters.txt') 

set title b

set output 'p0.5_at0.1_170212_1717_temperature_internal.png'

plot 'p0.5_at0.1_170212_1717_temperature_internal.dat' w l 

unset output 

