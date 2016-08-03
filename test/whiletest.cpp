#include <iostream>

int main()
{
    double dq = 0.05;
    double eps = 0.5;
    while(dq < 1)
    {
        while(eps < 2)
        {
            std::cout << dq << "\t" << eps << std::endl;
            eps += 0.01;
        }
        dq += 0.01;
        eps = 0.5;
    }
    return 0;
}
