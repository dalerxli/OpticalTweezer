#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstring>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_rng.h>
#include "../include/functions.hpp"
#include "../include/classes.hpp"
#include "../include/globals.hpp"


int main(int argc, char** argv)
{
    const gsl_rng_type* randtype = gsl_rng_taus;
    gsl_rng *randnum = gsl_rng_alloc(randtype);
    gsl_rng_set(randnum,98);
    gsl_histogram* h = gsl_histogram_alloc(3);
    gsl_histogram* g = gsl_histogram_alloc(3);
    double range[4] = { 1.0, 10.0, 100.0, 1000.0 };
    gsl_histogram_set_ranges(h,range,4);
    gsl_histogram_set_ranges_uniform(g,1.0,1000.0);
    for(int i=0;i<1000;i++)
    {
        gsl_histogram_increment(h,i);
        gsl_histogram_increment(g,i);
    }
    FILE* hist1 = fopen("histogram_unequal.dat","w");
    FILE* hist2 = fopen("histogram_equal.dat","w");
    gsl_histogram_fprintf(hist1,h,"%f","%f"); 
    gsl_histogram_fprintf(hist2,g,"%f","%f"); 
    fclose(hist1);
    fclose(hist2);

    gsl_histogram* k = gsl_histogram_alloc(100);
    gsl_histogram_set_ranges_uniform(k,0,10);
    for(int i=0;i<100000;i++)
    {
        double num = gsl_rng_uniform(randnum)*10.;
        gsl_histogram_increment(k,num);
    }
    int test = gsl_histogram_increment(k,1.);
    std::cout << test << std::endl;
    FILE* hist3 = fopen("histogram_test.dat","w");
    gsl_histogram_fprintf(hist3,k,"%.1f","%f"); 
    fclose(hist3);

	return 0;
}
