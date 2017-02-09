set term png 

b = system('cat  p1_at0.05_170208_0732_parameters.txt') 

set title b

set output 'p1_at0.05_170208_0732_temperature_internal.png'

plot 'p1_at0.05_170208_0732_temperature_internal.dat' w l 

unset output 

