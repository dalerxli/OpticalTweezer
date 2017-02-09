set term png 

b = system('cat  p1_at0.05_170208_0529_parameters.txt') 

set title b

set output 'p1_at0.05_170208_0529_temperature_internal.png'

plot 'p1_at0.05_170208_0529_temperature_internal.dat' w l 

unset output 

