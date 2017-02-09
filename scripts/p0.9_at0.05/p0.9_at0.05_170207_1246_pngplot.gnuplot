set term png 

b = system('cat  p0.9_at0.05_170207_1246_parameters.txt') 

set title b

set output 'p0.9_at0.05_170207_1246_temperature_internal.png'

plot 'p0.9_at0.05_170207_1246_temperature_internal.dat' w l 

unset output 

