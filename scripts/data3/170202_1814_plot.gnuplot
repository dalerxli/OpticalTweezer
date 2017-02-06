Tin=system("./get_number.sh GasIn 170202_1814_statistics_mass.dat")+0
Tout=system("./get_number.sh GasOut 170202_1814_statistics_mass.dat")+0
Tcm=system("./get_number.sh T_COM 170202_1814_statistics_mass.dat")+0
M=32
plot "170202_1814_histogram_mass_vx.dat" with boxes, \
    sqrt(M/(2.*pi*Tcm))*exp(-x**2*M/(2*Tcm)), \
    sqrt(M/(2.*pi*Tin))*exp(-x**2*M/(2*Tin)), \
    sqrt(M/(2.*pi*Tout))*exp(-x**2*M/(2*Tout))
