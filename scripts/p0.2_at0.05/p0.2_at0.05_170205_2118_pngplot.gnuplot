set term png 

b = system('cat  p0.2_at0.05_170205_2118_parameters.txt') 

set title b

set output 'p0.2_at0.05_170205_2118_temperature_internal.png'

plot 'p0.2_at0.05_170205_2118_temperature_internal.dat' w l 

unset output 

