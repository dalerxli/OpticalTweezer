Tin=0.100496285714
Tout=0.108879642857
T_COM=0.205265643202
MassNano = 32
MassGas = 0.1

plot "170201_1701_histogram_mass_vx.dat" with boxes, \
    sqrt(MassNano/(2*pi*T_COM))*exp(-x**2*MassNano/(2*T_COM)) lt rgb "black", \
     sqrt(MassGas/(2*pi*Tin))*exp(-x**2*MassGas/(2*Tin)) lt rgb "red",\
     sqrt(MassGas/(2*pi*Tout))*exp(-x**2*MassGas/(2*Tout)) lt rgb "blue" 

