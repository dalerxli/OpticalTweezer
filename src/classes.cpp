#include "../include/classes.hpp"


//constructors
Particle::Particle()
{
    r = new double[3];
    v = new double[3];      
    a = new double[3];
    rnew = new double[3];
    vnew = new double[3];
    vhalf = new double[3];
    r0 = new double[3];
    v0 = new double[3];
    m = 1.;
    potE = 0;
    

    for(int i=0;i<3;i++)
    {
            r[i]=0;
            v[i]=0;
            r0[i]=0;
            v0[i]=0;
            a[i]=0;
            rnew[i]=0;
            vnew[i]=0;
            vhalf[i]=0;
    }
    ID = g_ID;
    type = 0;
    g_ID++;
    name = "Glass";
}

Particle::Particle(double x,double y, double z, double v1, double v2, double v3)
{
    r = new double[3];
    v = new double[3];      
    a = new double[3];
    for(int i=0;i<3;i++)
            a[i]=0;
    r[0]=x;
    r[1]=y;
    r[2]=z;
    v[0]=v1;
    v[1]=v2;
    v[2]=v3;
    r0[0]=x;
    r0[1]=y;
    r0[2]=z;
    v0[0]=v1;
    v0[1]=v2;
    v0[2]=v3;
    ID = g_ID;
    type=2;
    g_ID++;
    m = 1.;
    potE = 0;
}

Particle::Particle(double x,double y, double z, double v1, double v2, double v3,std::string Name)
{
    r = new double[3];
    v = new double[3];      
    r0 = new double[3];
    v0 = new double[3];      
    a = new double[3];
    for(int i=0;i<3;i++)
            a[i]=0;
    r[0]=x;
    r[1]=y;
    r[2]=z;
    v[0]=v1;
    v[1]=v2;
    v[2]=v3;
    r0[0]=x;
    r0[1]=y;
    r0[2]=z;
    v0[0]=v1;
    v0[1]=v2;
    v0[2]=v3;
    ID = g_ID;
    type=2;
    g_ID++;
    name = Name;
    if(strcmp(name.c_str(),"GasIn")==0 || strcmp(name.c_str(),"GasOut")==0)
        m = 0.1;
    else
        m = 1.;
    potE = 0;
}

Particle::Particle(std::string Name)
{
    r = new double[3];
    v = new double[3];      
    a = new double[3];
    rnew = new double[3];
    vnew = new double[3];
    vhalf = new double[3];

    for(int i=0;i<3;i++)
    {
            r[i]=0;
            v[i]=0;
            a[i]=0;
            rnew[i]=0;
            vnew[i]=0;
            vhalf[i]=0;
    }
    ID = g_ID;
    type = 0;
    g_ID++;
    name = Name;
    if(strcmp(name.c_str(),"GasIn")==0 || strcmp(name.c_str(),"GasOut")==0)
        m = 0.1;
    else
        m = 1.;
    potE = 0;
}

