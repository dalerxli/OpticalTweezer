set term png 

b = system('cat  p0.7_at0.1_170208_1515_parameters.txt') 

set title b

set output 'p0.7_at0.1_170208_1515_temperature_internal.png'

plot 'p0.7_at0.1_170208_1515_temperature_internal.dat' w l 

unset output 

