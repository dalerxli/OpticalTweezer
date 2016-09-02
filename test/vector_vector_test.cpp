#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstring>
#include <vector>

int main(int argc, char** argv)
{
    std::vector<double> v1;
    std::vector<double> v3;

    typedef std::vector<double> dv;
    std::vector<dv> v2;
    v2.push_back(v1);
    v2[0].push_back(3.14);
    v2.push_back(v3);
    v2[1].push_back(5.23);
    std::cout << v2[0][0] << std::endl;
    return 0;
}
