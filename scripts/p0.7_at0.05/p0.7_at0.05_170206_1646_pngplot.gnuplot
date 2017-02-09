set term png 

b = system('cat  p0.7_at0.05_170206_1646_parameters.txt') 

set title b

set output 'p0.7_at0.05_170206_1646_temperature_internal.png'

plot 'p0.7_at0.05_170206_1646_temperature_internal.dat' w l 

unset output 

