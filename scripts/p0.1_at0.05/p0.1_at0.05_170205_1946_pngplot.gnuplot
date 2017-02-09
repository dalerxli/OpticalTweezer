set term png 

b = system('cat  p0.1_at0.05_170205_1946_parameters.txt') 

set title b

set output 'p0.1_at0.05_170205_1946_temperature_internal.png'

plot 'p0.1_at0.05_170205_1946_temperature_internal.dat' w l 

unset output 

