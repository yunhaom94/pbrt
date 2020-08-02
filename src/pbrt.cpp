#include "core/pbrt.h"
#include "core/spectrum.h"
#include <fstream>


#include "cameras/perspective.h"
#include "shapes/sphere.h"
#include "core/ray.h"
#include "core/interaction.h"
#include "core/transform.h"

bool ParseFile(std::string f)
{
    std::cout << f << std::endl;
    return false;
}

Vector3i ray_color(Ray ray, Shape *stuffs[100])
{
    // skyblue
    Vector3f color(135, 206, 235);


    SurfaceInteraction si, si_temp;
    bool ifhit = false , ifhit_temp;
    Float thit;
    Float tclosest = 99999999999;

    for (size_t i = 0; i < 2; i++)
    {
        Shape *o = stuffs[i];
        ifhit_temp = o->Intersect(ray, &thit, &si_temp, false);
     
        if (ifhit_temp)
        {
            if (thit < tclosest)
            {
                tclosest = thit;
                si = si_temp;
            }
            ifhit = true;
        }
    }


    if (ifhit)
    {
        color = Vector3f(1 * tclosest, 1, 1);
    }

    color.normalize();

    color *= 255.0;

    return (Vector3i)color;
}

int main(int argc, char* argv[])
{
    /*
    Options options;
    

    if (argc != 2) 
    {
        std::cout << "Must provide a file name" << std::endl;
    }
    else
    {
        std::string filename = argv[1];

        if (!ParseFile(filename))
        {
            std::cout << "File provide not correct" << std::endl;
        }
    }*/

    const auto aspect_ratio = 16.0 / 9.0;
    const int image_width = 400;
    const int image_height = static_cast<int>(image_width / aspect_ratio);

    // world
    Transform tm = Translate(Vector3f(0, 0, 0));
    Sphere sphere(&tm, &tm, false, 1, -1, 1, 360);

    Transform tm2 = Translate(Vector3f(0, 100, 0));
    Sphere sphere2(&tm2, &tm2, false, 100, -100, 100, 360);


    // camera
    auto viewport_height = 2.0;
    auto viewport_width = aspect_ratio * viewport_height;
    auto focal_length = 1.0;

    auto origin = Point3f(0, 0, 2);
    auto horizontal = Point3f(viewport_width, 0, 0);
    auto vertical = Point3f(0, viewport_height, 0);
    auto lower_left_corner = origin - horizontal / 2 - vertical / 2 - Point3f(0, 0, focal_length);

      

    std::ofstream outfile("out.ppm");

    outfile << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    for (int j = image_height - 1; j >= 0; --j) 
    {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) 
        {

            auto u = double(i) / (image_width - 1);
            auto v = double(j) / (image_height - 1);
            
            Ray ray(origin, lower_left_corner + u * horizontal + v * vertical - origin);
            
            Shape *o[100] = { &sphere, &sphere2 };
            Vector3i color = ray_color(ray, o);

            outfile << color[0] << ' ' << color[1] << ' ' << color[2] << '\n';
        }
    }

    return 0;
}

