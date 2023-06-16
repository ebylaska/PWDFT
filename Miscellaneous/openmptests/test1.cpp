#include <iostream>
int main()
{
    #pragma omp parallel
    {
        std::cout << "Hello World" << std::endl;
    }
    return 0;
}
