#include "../include/globals.hpp"
#include "../include/classes.hpp"
#include "../include/functions.hpp"

int main()
{
    //Particle* p1 = new Particle[N];
    //for(int i=0;i<N;i++)
        //p1[i] = new Particle();
/*
 *    int M=1;
 *    int n,x,y,z,k;
 *    n=0;
 *    
 *    while(4*M*M*M < N)
 *        ++M;
 *
 *    double a = L/M;
 *    double xCell[4] = {0.25,0.75,0.75,0.25};
 *    double yCell[4] = {0.25,0.75,0.25,0.75};
 *    double zCell[4] = {0.25,0.25,0.75,0.75};
 *    
 *    for(x=0;x<M;x++)
 *        for(y=0;y<M;y++)
 *            for(z=0;z<M;z++)
 *                for(k=0;k<4;k++)
 *                    if(n<N)
 *                    {
 *                        p1[n]->r[0] = (x + xCell[k])*a;
 *                        p1[n]->r[1] = (y + yCell[k])*a;
 *                        p1[n]->r[2] = (z + zCell[k])*a;
 *                        n++;
 *                    }
 */
    /*
     *InitPositions(p1);
     *FILE* p1data = fopen("p1data.dat","wb");
     *for(int i=0;i<N;i++)
     *{
     *    for(int j=0;j<3;j++)
     *        fwrite(&p1[i].r[j],sizeof(double),1,p1data);
     *    for(int j=0;j<3;j++)
     *        fwrite(&p1[i].v[j],sizeof(double),1,p1data);
     *    for(int j=0;j<3;j++)
     *        fwrite(&p1[i].a[j],sizeof(double),1,p1data);
     *}
     *FILE* testout1 = fopen("testread1.dat","w");
     *fclose(p1data);
     *
     *for(int n=0;n<N;n++)
     *{
     *    fprintf(testout1,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",p1[n].r[0],p1[n].r[1],p1[n].r[2],p1[n].v[0],p1[n].v[1],p1[n].v[2],p1[n].a[0],p1[n].a[1],p1[n].a[2]);
     *}
     */

    /*
     *FILE* p2data = fopen("p1data.dat","rb");
     *double testval = 0;
     *fread(&testval,sizeof(double),1,p2data); 
     *printf("%lf\n",testval);
     */
    Particle* p2 = new Particle[N];
    //for(int i=0;i<N;i++)
        //p2[i] = new Particle();
    FILE* p2data = fopen("p1data.dat","rb");
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<3;j++)
            fread(&(p2[i].r[j]),sizeof(double),1,p2data);
        for(int j=0;j<3;j++)
            fread(&(p2[i].v[j]),sizeof(double),1,p2data);
        for(int j=0;j<3;j++)
            fread(&(p2[i].a[j]),sizeof(double),1,p2data);
    }
    
    FILE* testout = fopen("testread.dat","w");
    
    for(int n=0;n<N;n++)
    {
        fprintf(testout,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",p2[n].r[0],p2[n].r[1],p2[n].r[2],p2[n].v[0],p2[n].v[1],p2[n].v[2],p2[n].a[0],p2[n].a[1],p2[n].a[2]);
    }

    fcloseall();

    return 0;
}
