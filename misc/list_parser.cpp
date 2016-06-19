#include <iostream>
#include <cstdio>
#include <string>
#include <list>
#include <iterator>


void ManipulateList(std::list<int>& list)
{
    for(std::list<int>::iterator i1 = list.begin(); i1!=list.end();i1++)
    {
        (*i1) *= 2;
    }
}

void DontManipulateList(std::list<int> list)
{
    for(std::list<int>::iterator i1 = list.begin(); i1!=list.end();i1++)
    {
        (*i1) *= 2;
    }
}


void PrintList(std::list<int> list)
{
    for(std::list<int>::iterator i1 = list.begin(); i1!=list.end();i1++)
    {
        std::cout << (*i1) << std::endl;
    }
}


int main()
{
    std::list<int> list1;
    for(int i=0;i<20;i++)
    {
        list1.push_back(i);
    }



    PrintList(list1);

    DontManipulateList(list1);

    PrintList(list1);

    ManipulateList(list1);

    PrintList(list1);

    return 0;
}
