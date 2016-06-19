#include <iostream>
#include <vector>
#include <string>

int main()
{
    std::vector<int> numbers;
    for(int i=0;i<20;i++)
        numbers.push_back(i);
    std::cout << std::string(20,'-') << std::endl;
    std::cout << "Before erase" << std::endl; 
    std::cout << std::string(20,'-') << std::endl;
    for(int i=0;i<numbers.size();i++)
        std::cout << "numbers[" << i << "]: " << numbers[i] << std::endl;
    for(int i=0;i<numbers.size();i++)
        if(numbers[i] == 5)
            numbers.erase(numbers.begin()+i);
    std::cout << std::string(20,'-') << std::endl;
    std::cout << "After erase" << std::endl; 
    std::cout << std::string(20,'-') << std::endl;
    for(int i=0;i<numbers.size();i++)
        std::cout << "numbers[" << i << "]: " << numbers[i] << std::endl;
    return 0;
}
