set term png 

b = system('cat  p0.7_at0.05_170207_1447_parameters.txt') 

set title b

set output 'p0.7_at0.05_170207_1447_temperature_internal.png'

plot 'p0.7_at0.05_170207_1447_temperature_internal.dat' w l 

unset output 

