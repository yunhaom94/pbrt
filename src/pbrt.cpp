#include <iostream>
#include <vector>
#include <Eigen/Core>

#include "core/pbrt.h"
#include "core/spectrum.h"
#include "core/geometry.h"


void pbrtInit() { }
void pbrtCleanup() { }

int main() 
{

    Eigen::Matrix4d d;
    Eigen::Vector3d v(1, 2, 3);

    //Transform s;

    d << 1, 2, 3, 4,
        5, 6, 7, 8, 
        9, 10, 11, 12,
        13, 14, 15, 16;


    //std::cout << Eigen::Vector4d(v, 0);

    std::vector<std::string> filenames;

    Spectrum s = Spectrum(1);


    pbrtCleanup();
    std::cout << "Hello World!";
    return 0;
}

// TODO:
int main_real(int argc, char** argv)
{

    return 0;
}