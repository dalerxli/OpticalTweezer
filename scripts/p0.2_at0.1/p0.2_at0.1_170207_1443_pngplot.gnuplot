set term png 

b = system('cat  p0.2_at0.1_170207_1443_parameters.txt') 

set title b

set output 'p0.2_at0.1_170207_1443_temperature_internal.png'

plot 'p0.2_at0.1_170207_1443_temperature_internal.dat' w l 

unset output 

