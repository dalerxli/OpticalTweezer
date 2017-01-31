Tin=system("./get_number_osx.sh GasIn 170130_2359_statistics.dat")+0
Tout=system("./get_number_osx.sh GasOut 170130_2359_statistics.dat")+0
Tcm=system("./get_number_osx.sh T_COM 170130_2359_statistics.dat")+0
plot "170130_2359_histogram_vx.dat" with boxes, \
    sqrt(1./(2.*pi*Tcm))*exp(-x**2/(2*Tcm)), \
    sqrt(0.5/(2.*pi*Tin))*exp(-x**2*0.5/(2*Tin)), \
    sqrt(0.5/(2.*pi*Tout))*exp(-x**2*0.5/(2*Tout))
