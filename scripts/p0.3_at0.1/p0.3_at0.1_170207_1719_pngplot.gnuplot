set term png 

b = system('cat  p0.3_at0.1_170207_1719_parameters.txt') 

set title b

set output 'p0.3_at0.1_170207_1719_temperature_internal.png'

plot 'p0.3_at0.1_170207_1719_temperature_internal.dat' w l 

unset output 

