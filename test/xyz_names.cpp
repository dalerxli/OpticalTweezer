#include <iostream>
#include <cstdio>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <list>
#include <iterator>
#include <vector>
#include <sys/stat.h>
#include <unistd.h>
#include <fstream>
#include <iomanip>

int main()
{
    const gsl_rng_type* T = gsl_rng_taus;
    gsl_rng* RND = gsl_rng_alloc(T);
    gsl_rng_set(RND,10);
    FILE* output = fopen("xyz_test.xyz","w");
    char type1[] = "Hanspeter";
    //char type2[] = "null";
    std::string type2 = "null";
    fprintf(output,"20\nName Test\n");
    for(int i=0;i<10;i++)
        fprintf(output,"%s\t%lf\t%lf\t%lf\n",type1,gsl_rng_uniform(RND),gsl_rng_uniform(RND),gsl_rng_uniform(RND));
    for(int i=0;i<10;i++)
        fprintf(output,"%s\t%lf\t%lf\t%lf\n",type2.c_str(),gsl_rng_uniform(RND),gsl_rng_uniform(RND),gsl_rng_uniform(RND));
    
    fclose(output);

    return 0;
}
