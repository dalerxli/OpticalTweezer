#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>

char* input;

int main(int argc,char** argv)
{
    if(argc == 1)
        puts("There were no arguments given");
    if(argc > 1)
    {
        if(std::strcmp(argv[1],"init") == 0)
        {
            puts("Initializing system");
            input = argv[1];
        }
        else if(std::strcmp(argv[1],"heat") == 0)
        {
            puts("Heating system");
            input = argv[1];
        }
        else if(std::strcmp(argv[1],"baro") == 0)
        {
            puts("Running barostat");
            input = argv[1];
        }
        else
        {
            puts("Input not recognized. Options are 'init', 'heat' and 'baro'");
            return 1;
        }
    }
    if(input != NULL)
    {
        std::string output(input);
        std::cout << output << std::endl;
    }
    return 0;
}
