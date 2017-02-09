set term png 

b = system('cat  p0.6_at0.1_170207_2212_parameters.txt') 

set title b

set output 'p0.6_at0.1_170207_2212_temperature_internal.png'

plot 'p0.6_at0.1_170207_2212_temperature_internal.dat' w l 

unset output 

