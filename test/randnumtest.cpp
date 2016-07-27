#include "../include/globals.hpp"

int main()
{
    gsl_rng_set(r,98);
    std::cout << "L: " << L << std::endl;
    FILE* output = fopen("randnumtest.dat","w");
    for(int i=0;i<100000;i++)
        //fprintf(output,"%lf\n",2*(gsl_rng_uniform(r)-0.5));
        fprintf(output,"%lf\n",(2*(gsl_rng_uniform(r)-0.5)*L)+0.5*L);
    return 0;
}
