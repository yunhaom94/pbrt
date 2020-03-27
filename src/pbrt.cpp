#include <iostream>
#include <vector>
#include <Eigen/Core>
#include "core/bounding_boxes.h"

void pbrtInit() { }
void pbrtCleanup() { }

int main() {

    Eigen::Vector3i p1(0, 0, 0);

    std::cout << p1.x() << std::endl;

    Bounds3<int> b(p1, p1);

    std::vector<std::string> filenames;

    //TODO: Process command-line arguments
    pbrtInit();
    //TODO: Process scene description

    pbrtCleanup();
    std::cout << "Hello World!";
    return 0;
}