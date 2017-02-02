Tin=system("./get_number_osx.sh GasIn 170201_1303_1_statistics.dat")+0
Tout=system("./get_number_osx.sh GasOut 170201_1303_1_statistics.dat")+0
Tcm=system("./get_number_osx.sh T_COM 170201_1303_1_statistics.dat")+0
plot "170201_1303_1_histogram_vx.dat" with boxes, \
    sqrt(1./(2.*pi*Tcm))*exp(-x**2/(2*Tcm)), \
    sqrt(10./(2.*pi*Tin))*exp(-x**2*10./(2*Tin)), \
    sqrt(10./(2.*pi*Tout))*exp(-x**2*10./(2*Tout))
