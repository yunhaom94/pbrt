#include <iostream>
#include <vector>
#include <Eigen/Core>

#include "core/pbrt.h"
#include "core/spectrum.h"
#include "core/options.h"

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
int main(int argc, char* argv[])
{
    Options options;
    std::vector<std::string> filenames;

    if (filenames.size() == 0) {
        ParseFile("-");
    }
    else
    {

        for (const std::string& f : filenames)
        {
            if (!ParseFile(f))
                Error("Couldn¡¯t open scene file \"%s\"", f.c_str());
        }
    }

    return 0;
}

