#include <iostream>
#include <vector>
#include "Common.h"
#include "Scene.h"
#include "Camera.h"
#include "Material.h"
#include "Hittable.h"
#include "Utils/lodepng.h"

const int kMaxTraceDepth = 5;

Color TraceRay(const Ray& ray,
               const std::vector<LightSource>& light_sources,
               const Hittable& scene,
               int trace_depth);

Color Shade(const std::vector<LightSource>& light_sources,
            const Hittable& hittable_collection,
            const HitRecord& hit_record,
            int trace_depth) {
    // TODO: Add your code here.
    Color color(0.f, 0.f, 0.f);

    HitRecord record = hit_record;
    Material material = record.material;

    color = material.ambient;

    for (auto& light : light_sources) {
        auto shadow_ray = light.position - record.position; // maybe wrong, here i use the [the light-position] minus [intersect position] to get [intersection-point to light vector]
        if (glm::dot(record.normal, glm::normalize(shadow_ray)) > 1e-5f) {
            // if not hit
            if (!hittable_collection.Hit(shadow_ray, &record)) {
                // fill in code                                                                                         
                Vec shadow_reflection_dir = glm::normalize(2.0f * glm::dot(glm::normalize(shadow_ray), record.normal) * record.normal - glm::normalize(shadow_ray));
                color = color + (1.0f * light.intensity * (material.k_d * material.diffuse * (glm::dot(record.normal, glm::normalize(shadow_ray))) + material.k_s * material.specular * (std::pow(glm::dot(shadow_reflection_dir, glm::normalize(-record.in_direction)), material.sh))));
            }
        }
    }

    if (trace_depth < kMaxTraceDepth) {
        if (material.k_s > 1e-5f) {
            auto reflected_ray = record.reflection;
            Color r_color = TraceRay(reflected_ray, light_sources, hittable_collection, trace_depth + 1);
            r_color = r_color * material.k_s;
            color = color + r_color;
        }
    }

    if (color.r > 1.0f) {
        color.r = 1.0f;
    }
    if (color.g > 1.0f) {
        color.g = 1.0f;
    }
    if (color.b > 1.0f) {
        color.b = 1.0f;
    }

    return color;
}

Color TraceRay(const Ray& ray,
               const std::vector<LightSource>& light_sources,
               const Hittable& hittable_collection,
               int trace_depth) {
    // TODO: Add your code here.
    HitRecord record;
    Color color(0.0f, 0.0f, 0.0f);
    if (hittable_collection.Hit(ray, &record)) {
        color = Shade(light_sources, hittable_collection, record, trace_depth);
        return color;
    }
    else
    {
        return color;
    }

    return color;
}

int main() {
    // TODO: Set your workdir (absolute path) here.
    const std::string work_dir("C:/Users/wei/Desktop/hku_course_material/COMP3271/COMP3271/PA3/Graphics_PA3_Release_2/Graphics_PA3_Release/");

    // Construct scene
    Scene scene(work_dir, "scene/spheres.toml");
    const Camera& camera = scene.camera_;
    int width = camera.width_;
    int height = camera.height_;

    std::vector<unsigned char> image(width * height * 4, 0);

    float progress = 0.f;

    // Traverse all pixels
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            Color color(0.f, 0.f, 0.f);
            int count = 0;
            for (float bias_x = 0.25f; bias_x < 1.f; bias_x += .5f) {
                for (float bias_y = 0.25f; bias_y < 1.f; bias_y += .5f) {
                    Ray ray = camera.RayAt(float(x) + bias_x, float(y) + bias_y);
                    color += TraceRay(ray, scene.light_sources_, scene.hittable_collection_, 1);
                    count++;
                }
            }
            color /= float(count);
            int idx = 4 * ((height - y - 1) * width + x);
            for (int i = 0; i < 3; i++) {
                image[idx + i] = (uint8_t) (glm::min(color[i], 1.f - 1e-5f) * 256.f);
            }
            image[idx + 3] = 255;

            float curr_progress = float(x * height + y) / float(height * width);
            if (curr_progress > progress + 0.05f) {
                progress += 0.05f;
                std::cout << "Progress: " << progress << std::endl;
            }
        }
    }

    // Save result as png file
    std::vector<unsigned char> png;
    unsigned error = lodepng::encode(png, image, width, height);
    lodepng::save_file(png, work_dir + "output.png");
}
