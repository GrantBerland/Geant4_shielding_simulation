#include <iostream>
#include <ctime>
#include <cstdlib>

int main(int argc, char* argv[]){


time_t systime = time(NULL);

std::cout << systime << std::endl;


return 0;
}
