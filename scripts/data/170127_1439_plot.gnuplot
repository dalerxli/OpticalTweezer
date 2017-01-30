Tin=system("./get_number.sh GasIn 170127_1439_statistics.dat")+0
Tout=system("./get_number.sh GasOut 170127_1439_statistics.dat")+0
Tcm=system("./get_number.sh T_COM 170127_1439_statistics.dat")+0
plot "170127_1439_histogram_vx.dat" with boxes, \
    sqrt(1./(2.*pi*Tcm))*exp(-x**2/(2*Tcm)), \
    sqrt(1./(2.*pi*Tin))*exp(-x**2/(2*Tin)), \
    sqrt(1./(2.*pi*Tout))*exp(-x**2/(2*Tout))
