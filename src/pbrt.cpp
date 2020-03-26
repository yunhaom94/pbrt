#include <iostream>
#include <vector>
#include <Eigen/Core>
#include "core/ray.h"

void pbrtInit() { }
void pbrtCleanup() { }

int main() {

    std::vector<std::string> filenames;

    //TODO: Process command-line arguments
    pbrtInit();
    //TODO: Process scene description

    pbrtCleanup();
    std::cout << "Hello World!";
    return 0;
}