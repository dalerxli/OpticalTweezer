set term png 

b = system('cat  p0.9_at0.05_170207_1844_parameters.txt') 

set title b

set output 'p0.9_at0.05_170207_1844_temperature_internal.png'

plot 'p0.9_at0.05_170207_1844_temperature_internal.dat' w l 

unset output 

