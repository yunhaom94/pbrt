#include "core/pbrt.h"
#include "core/spectrum.h"
#include <fstream>


#include "cameras/perspective.h"
#include "shapes/sphere.h"
#include "core/ray.h"
#include "core/interaction.h"
#include "core/transform.h"
#include "core/bounding_boxes.h"
#include "core/film.h"

bool ParseFile(std::string f)
{
    std::cout << f << std::endl;
    return false;
}

Vector3i ray_color(Ray ray, Shape *stuffs[100])
{
    // skyblue
    Vector3f color(164, 217, 244);


    SurfaceInteraction si, si_temp;
    bool ifhit = false , ifhit_temp;
    Float thit;
    Float tclosest = 99999999999;
    int ihit = 0;

    for (size_t i = 0; i < 2; i++)
    {
        Shape *o = stuffs[i];
        ifhit_temp = o->Intersect(ray, &thit, &si_temp, false);
     
        if (ifhit_temp)
        {
            if (thit < tclosest)
            {
                ihit = i;
                tclosest = thit;
                si = si_temp;

            
            }
            ifhit = true;
        }
    }


    if (ifhit)
    {
        color = Vector3f(0, 0, 0);
        color[ihit] = 1;
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
    Transform tm = Translate(Vector3f(0, 0, 1));
    Sphere sphere(&tm, &tm, false, 0.5, 0.5, -0.5, 360);

    Transform tm2 = Translate(Vector3f(0, 100.5, 1));
    Sphere sphere2(&tm2, &tm2, false, 100, -100.5, 100.5, 360);

    // camera
    Transform tm3 = Translate(Vector3f(0, 0, 2));
    AnimatedTransform camtm(&tm3, 0, &tm3, 1);

    Bounds2f screenwindow = Bounds2f(Point2f(0, 0), Point2f(1, 1));
    std::unique_ptr<Filter> filter = std::unique_ptr<Filter>(new BoxFilter(Vector2f(0.5, 0.5)));;
    Film film(Point2i(image_width, image_height), 
                Bounds2f(Point2f(0, 0), Point2f(1, 1)), 
                std::move(filter),
                100.0, 
                std::string(""), 
                1.0);
    AnimatedTransform identity(new Transform, 0, new Transform, 1);
    PerspectiveCamera pcam = PerspectiveCamera(camtm, screenwindow, 0, 1, 0, 1, 90, &film, nullptr);

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
            
            CameraSample cms;
            cms.pFilm = Point2f(i, j);
            cms.pLens = Point2f(0, 0);
            cms.time = 0;

            Ray r;
            pcam.GenerateRay(cms, &r);

            Ray ray(origin, lower_left_corner + u * horizontal + v * vertical - origin);
            
            //std::cout << r.o << " +\n " << r.d << std::endl;
            //std::cout << "`````" << std::endl;
            //std::cout << ray.o << " +\n " << ray.d << std::endl;
            //std::cout << " ====== " << std::endl;

            Shape *o[100] = { &sphere, &sphere2 };
            Vector3i color = ray_color(r, o);

            outfile << color[0] << ' ' << color[1] << ' ' << color[2] << '\n';
        }
    }

    return 0;
}

