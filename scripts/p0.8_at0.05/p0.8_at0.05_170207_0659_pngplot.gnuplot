set term png 

b = system('cat  p0.8_at0.05_170207_0659_parameters.txt') 

set title b

set output 'p0.8_at0.05_170207_0659_temperature_internal.png'

plot 'p0.8_at0.05_170207_0659_temperature_internal.dat' w l 

unset output 

