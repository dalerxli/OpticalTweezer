set term png 

b = system('cat  p1_at0.1_170211_1858_parameters.txt') 

set title b

set output 'p1_at0.1_170211_1858_temperature_internal.png'

plot 'p1_at0.1_170211_1858_temperature_internal.dat' w l 

unset output 

