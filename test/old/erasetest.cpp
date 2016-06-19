#include <iostream>
#include <cstdio>
#include <cmath>
#include <list>
#include <iterator>
#include <vector>

int main()
{
    std::list<int> list1;
    for(int i=0;i<100;i++)
    {
        list1.push_back(i);
    }

    std::list<int>::iterator iter1;
    std::list<int>::iterator iter2;
    for(iter1 = list1.begin();iter1 != list1.end();iter1++)
    {
        std::cout << *iter1 << ", ";
        //if(*iter1==50)
            //iter1=list1.erase(iter1);
    }
    std::cout << std::endl;

    for(iter1 = list1.begin();iter1 != list1.end();iter1++)
    {
        if(*iter1==50)
            iter1=list1.erase(iter1);
    }

    for(iter1 = list1.begin();iter1 != list1.end();iter1++)
    {
        std::cout << *iter1 << ", ";
    }
    std::cout << std::endl; 

    int a;
    int b = 5;
    
    for(a=10;a<b;a++)
    {
        std::cout << a << std::endl;
    }
    return 0;
}
