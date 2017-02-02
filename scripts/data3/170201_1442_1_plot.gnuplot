Tin=system("./get_number_osx.sh GasIn 170201_1442_1_statistics.dat")+0
Tout=system("./get_number_osx.sh GasOut 170201_1442_1_statistics.dat")+0
Tcm=system("./get_number_osx.sh T_COM 170201_1442_1_statistics.dat")+0
Tin = Tin * 0.1
Tout = Tout * 0.1
plot "170201_1442_1_histogram_vx.dat" with boxes, \
    sqrt(1./(2.*pi*Tcm))*exp(-x**2/(2*Tcm)), \
    sqrt(1./(2.*pi*Tin))*exp(-x**2*1./(2*Tin)), \
    sqrt(1./(2.*pi*Tout))*exp(-x**2*1./(2*Tout))
