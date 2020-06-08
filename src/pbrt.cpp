#include <iostream>
#include <vector>
#include <Eigen/Core>

#include "core/pbrt.h"
#include "core/spectrum.h"
#include "core/geometry.h"
#include "core/options.h"

int main()
{

    Matrix4x4 m;

    m << 1, 2, 3, 4,
        5, 6, 7, 8,
        9, 10, 11, 12,
        0, 0, 0, 1;

    Eigen::Vector4d v;
    v << 2, 3, 4, 0;


    std::cout << v;

    v.x() = 5;

    std::cout << v;

    return 0;
}

void pbrtInit() 
{
    // TODO: p1109
}
void pbrtCleanup() 
{
    // TODO: p1109
}

bool ParseFile(std::string f)
{
    
}

void Error(std::string s, std::string e)
{

}

// TODO:
int main_real(int argc, char** argv)
{
    Options options;
    std::vector<std::string> filenames;

    if (filenames.size() == 0) {
        ParseFile("-");
    }
    else {
    for (const std::string& f : filenames)
        if (!ParseFile(f))
            Error("Couldn¡¯t open scene file \"%s\"", f.c_str());
    }

    return 0;
}

