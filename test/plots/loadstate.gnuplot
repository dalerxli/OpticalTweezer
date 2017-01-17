set term png
numlines(file)=system("cat ".file." | wc -l")
folder="~/codes/OpticalTweezer/output/runs/"
folders=folders=system("ls ".folder)
do for [i in folders] {
    set output i.".png"
    params=system("cat ".folder."".i."/parameters.txt")
    set title params
    if(numlines(folder."".i."/temperature_internal.dat") > 0) {
        plot folder."".i."/temperature_internal.dat" w l
}
#print numlines(folder."".i."/temperature_internal.dat")
#print folder."".i."/temperature_internal.dat
}

