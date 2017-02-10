set term png 

b = system('cat  p0.1_at0.05_170210_1222_parameters.txt') 

set title b

set output 'p0.1_at0.05_170210_1222_temperature_internal.png'

plot 'p0.1_at0.05_170210_1222_temperature_internal.dat' w l 

unset output 

