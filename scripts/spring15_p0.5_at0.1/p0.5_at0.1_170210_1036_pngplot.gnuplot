set term png 

b = system('cat  p0.5_at0.1_170210_1036_parameters.txt') 

set title b

set output 'p0.5_at0.1_170210_1036_temperature_internal.png'

plot 'p0.5_at0.1_170210_1036_temperature_internal.dat' w l 

unset output 

