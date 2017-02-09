set term png 

b = system('cat  p0.4_at0.05_170205_2232_parameters.txt') 

set title b

set output 'p0.4_at0.05_170205_2232_temperature_internal.png'

plot 'p0.4_at0.05_170205_2232_temperature_internal.dat' w l 

unset output 

