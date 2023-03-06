#include <iostream>
#include "fstream"
#include "rtweekend.h"
#include "color.h"
#include "hittable_list.h"
#include "sphere.h"
#include "camera.h"
#include "json.hpp"
#include "material.h"
#include "map"
using json = nlohmann::json;
color ray_color(const ray &r, const hittable &world, int depth)
{
    hit_record rec;
    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0)
        return color(0, 0, 0);
    if (world.hit(r, 0.001, infinity, rec))
    {
        ray scattered;
        color attenuation;
        if (rec.mat_ptr->scatter(r, rec, attenuation, scattered))
            return attenuation * ray_color(scattered, world, depth - 1);
        return color(0, 0, 0);
    }
    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5 * (unit_direction.y() + 1.0);
    return (1.0 - t) * color(1.0, 1.0, 1.0) + t * color(0.5, 0.7, 1.0);
}
camera createCamera(json data, double aspect_ratio)
{
    // defult values
    point3 lookfrom(13, 2, 3);
    point3 lookat(0, 0, 0);
    vec3 vup(0, 1, 0);
    auto focus_dist = 10.0;
    auto aperture = 0.1;
    auto dist_to_focus = 10.0;
    double vfov = 20;
    auto tmp_lookfrom = data["lookfrom"];
    if (!tmp_lookfrom.is_null())
    {
        lookfrom = point3(tmp_lookfrom[0], tmp_lookfrom[1], tmp_lookfrom[2]);
    }
    auto tmp_lookat = data["lookat"];
    if (!tmp_lookat.is_null())
    {
        lookat = point3(tmp_lookat[0], tmp_lookat[1], tmp_lookat[2]);
    }
    auto tmp_vup = data["vup"];
    if (!tmp_vup.is_null())
    {
        vup = vec3(tmp_vup[0], tmp_vup[1], tmp_vup[2]);
    }
    if (!data["vfov"].is_null())
    {
        vfov = data["vfov"]; // vertical field-of-view in degrees
    }
    if (!data["aspect_ratio"].is_null())
        aspect_ratio = data["aspect_ratio"];
    if (!data["aperture"].is_null())
        aperture = data["aperture"];
    if (!data["focus_dist"].is_null())
        focus_dist = data["focus_dist"];

    return camera(lookfrom, lookat, vup, vfov, aspect_ratio, aperture, focus_dist);
}
hittable_list createWorld(json data)
{
    hittable_list world;
    std::map<std::string, int> materials;
    materials["lambertian"]=1;
    materials["metal"]=2;
    materials["dielectric"]=3;
    int i=0;
    for (json::iterator it = data.begin(); it != data.end(); ++it)
    {  
        std::cout<<"sphare number "<<i<<std::endl; 
        // point3 cen, double r, shared_ptr<material> m
        point3 centerPoint(random_double(-11, 11) + 0.9 * random_double(), 0.2, random_double(-11, 11) + 0.9 * random_double());
        double radios = random_double(0.1, 0.7);
        auto center = (*it)["center-point"];
        if (!center.is_null())
        {
            centerPoint = point3(center[0], center[1], center[2]);
        }
        auto r = (*it)["radios"];
        if (!r.is_null())
        {
            radios = r;
        }
        shared_ptr<material> sphere_material;
        auto mat = (*it)["material"];
        if (!mat.is_null())
        {
            color c;
            std::string matType = mat["type"];
            
            if (materials[matType]==1){
                //lambertian
                auto col=mat["color"];
                sphere_material=make_shared<lambertian>(color(col[0],col[1],col[2]));
            }else if(materials[matType]==2){
             //metal
                auto col=mat["color"];
                auto albedo =color(col[0],col[1],col[2]);
                auto fuzz = mat["fuzz"];
                sphere_material = make_shared<metal>(albedo, fuzz);
            }else{
                //dielectric
                auto index_of_refraction=mat["index"];
                 sphere_material = make_shared<dielectric>(index_of_refraction);
            }
        }
        else
        {
            auto choose_mat = random_double();
            if (choose_mat < 0.8)
            {
                // diffuse
                auto albedo = color::random() * color::random();
                sphere_material = make_shared<lambertian>(albedo);
            }
            else if (choose_mat < 0.95)
            {
                // metal
                auto albedo = color::random(0.5, 1);
                auto fuzz = random_double(0, 0.5);
                sphere_material = make_shared<metal>(albedo, fuzz);
            }
            else
            {
                // glass
                sphere_material = make_shared<dielectric>(1.5);
            }
        }
        world.add(make_shared<sphere>(centerPoint, radios, sphere_material));
        i++;
    }
    return world;
}

hittable_list random_scene()
{
    hittable_list world;
    auto ground_material = make_shared<lambertian>(color(0.5, 0.5, 0.5));
    world.add(make_shared<sphere>(point3(0, -1000, 0), 1000, ground_material));

    for (int a = -11; a < 11; a++)
    {
        for (int b = -11; b < 11; b++)
        {
            auto choose_mat = random_double();
            point3 center(a + 0.9 * random_double(), 0.2, b + 0.9 * random_double());

            if ((center - point3(4, 0.2, 0)).length() > 0.9)
            {
                shared_ptr<material> sphere_material;

                if (choose_mat < 0.8)
                {
                    // diffuse
                    auto albedo = color::random() * color::random();
                    sphere_material = make_shared<lambertian>(albedo);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
                else if (choose_mat < 0.95)
                {
                    // metal
                    auto albedo = color::random(0.5, 1);
                    auto fuzz = random_double(0, 0.5);
                    sphere_material = make_shared<metal>(albedo, fuzz);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
                else
                {
                    // glass
                    sphere_material = make_shared<dielectric>(1.5);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
            }
        }
    }

    auto material1 = make_shared<dielectric>(1.5);
    world.add(make_shared<sphere>(point3(0, 1, 0), 1.0, material1));

    auto material2 = make_shared<lambertian>(color(0.4, 0.2, 0.1));
    world.add(make_shared<sphere>(point3(-4, 1, 0), 1.0, material2));

    auto material3 = make_shared<metal>(color(0.7, 0.6, 0.5), 0.0);
    world.add(make_shared<sphere>(point3(4, 1, 0), 1.0, material3));

    return world;
}
int main(int length, char **path)
{
    json jsonData;
    bool jsonExsit = false;
    std::string fileName = "image.ppm";
    if (length > 0)
    {
        std::ifstream dataFile(path[1]);
        
        if (dataFile)
        {
            jsonData = json::parse(dataFile);
            jsonExsit = true;
            dataFile.close();
        }
    }
    const auto aspect_ratio = 3.0 / 2.0;
    const int image_width = 1200;
    const int image_height = static_cast<int>(image_width / aspect_ratio);
    const int samples_per_pixel = 500;
    const int max_depth = 50;
    point3 lookfrom(13, 2, 3);
    point3 lookat(0, 0, 0);
    vec3 vup(0, 1, 0);
    auto dist_to_focus = 10.0;
    auto aperture = 0.1;
    auto cam = camera(lookfrom, lookat, vup, 20, aspect_ratio, aperture, dist_to_focus);
    auto world = random_scene();
    if (jsonExsit)
    {
        fileName = jsonData["file-name"];
        auto c=jsonData["camara"];
        if(!c.is_null()){
            cam = createCamera(c, aspect_ratio);
        }
        auto w=jsonData["world"];
        if(!w.is_null()){
            world = createWorld(jsonData["world"]);
        }
    }
    // Render
    std::ofstream file;
    std::string file_name = "image.ppm";
    file.open(file_name);
    file << "P3\n"
         << image_width << " " << image_height << "\n255\n";

    for (int j = image_height - 1; j >= 0; --j)
    {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i)
        {
            color pixel_color(0, 0, 0);
            for (int s = 0; s < samples_per_pixel; ++s)
            {
                auto u = (i + random_double()) / (image_width - 1);
                auto v = (j + random_double()) / (image_height - 1);
                ray r = cam.get_ray(u, v);
                pixel_color += ray_color(r, world, max_depth);
            }
            write_color(file, pixel_color, samples_per_pixel);
        }
    }

    std::cout << "\nDone.\n";
    file.close();
}