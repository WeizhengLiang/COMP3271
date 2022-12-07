#include "Hittable.h"

// Sphere
bool Sphere::Hit(const Ray& ray, HitRecord* hit_record) const {
    // TODO: Add your code here.
    bool ret = false;

    float A = glm::dot(ray.d, ray.d);
    float B = 2.0f * glm::dot(ray.o, ray.d);
    float C = glm::dot(ray.o, ray.o) - std::pow(r_, 2.f);

    std::cout << glm::dot(ray.o, ray.o) << std::endl;


    float t1 = (-B - std::sqrt(std::pow(B, 2.f) - 4.0f * A * C)) / (2.0f * A);
    float t2 = (-B + std::sqrt(std::pow(B, 2.f) - 4.0f * A * C)) / (2.0f * A);

    std::cout <<"t1: "<< t1 <<"t2: "<< t2 << std::endl;

    float t;

    if (t1 > 1e-5f && t2 > 1e-5f) {
        ret = true;
        t = std::min(t1, t2);
        std::cout << "1 t(sephere)" << std::endl;
    }
    else if (t1 == 0.f && t2 == 0.f) {
        ret = true;
        t = .0f;
        std::cout << "2 t(sephere)" << std::endl;
    }
    else {
        std::cout << "something wrong about the t(sephere)" << std::endl;
    }
    //if (determinant > 1e-5f) {
    //    if (t1 > 1e-5f && t2 > 1e-5f) {
    //        ret = true;
    //        t = std::min(t1, t2);
    //    }

    //    if (t1 == t2 && t1 >= -1e-5f && t1 <= 1e-5f) {
    //        ret = true;
    //        t = 0.f;
    //    }
    //}

    if (ret) {
        // hit_record->... = ...
        // hit_record->... = ...
        // ....

        Point intersectPoint = ray.At(t);
        Vec normal = glm::normalize(intersectPoint - o_);

        hit_record->position = intersectPoint;
        hit_record->normal = normal;
        hit_record->distance = std::sqrt(std::pow(intersectPoint.x - ray.o.x, 2.f) + std::pow(intersectPoint.y - ray.o.y, 2.f) + std::pow(intersectPoint.z - ray.o.z, 2.f));
        hit_record->in_direction = ray.d;
        hit_record->reflection = glm::normalize(ray.d - 2.0f * glm::dot(ray.d, normal) * normal);
        hit_record->material = material_;
    }
    return ret;
}

// Quadric
bool Quadric::Hit(const Ray& ray, HitRecord* hit_record) const {
    // TODO: Add your code here.
    bool ret = false;

    glm::vec4 D(ray.d[0], ray.d[1], ray.d[2], .0f);
    glm::vec4 O(ray.o[0], ray.o[1], ray.o[2], 1.f);

    float A = glm::dot(D, A_ * D);
    float B = 2.0f * (glm::dot(O, A_ * D));
    float C = glm::dot(O, A_ * O);

    float determinant = std::pow(B, 2.f) - 4.f * A * C;

    float t;

    if (determinant > 1e-5f) {
        ret = true;
        float t1 = (-B - std::sqrt(std::pow(B, 2.f) - 4.0f * A * C)) / (2.0f * A);
        float t2 = (-B + std::sqrt(std::pow(B, 2.f) - 4.0f * A * C)) / (2.0f * A);
        t = std::min(t1, t2);
    }
    else if (determinant == 1e-5f) {
        ret = true;
        t = -B / (2.f * A);
    }
    else {
        ret = false;
        /*std::cout << "no feasible t(quadric)" << std::endl;*/
    }

    if (ret) {
        // hit_record->... = ...
        // hit_record->... = ...
        // ....

        Point intersectPoint = ray.At(t);
        Vec normal = glm::normalize((A_ + glm::transpose(A_)) * (O + D * t));

        hit_record->position = intersectPoint;
        hit_record->normal = normal;
        hit_record->distance = std::sqrt(std::pow(intersectPoint.x - ray.o.x, 2.f) + std::pow(intersectPoint.y - ray.o.y, 2.f) + std::pow(intersectPoint.z - ray.o.z, 2.f));
        hit_record->in_direction = ray.d;
        hit_record->reflection = glm::normalize(ray.d - 2.0f * glm::dot(ray.d, normal) * normal);
        hit_record->material = material_;
    }
    return ret;
}

// Triangle
bool Triangle::Hit(const Ray& ray, HitRecord* hit_record) const {
    // TODO: Add your code here.
    bool ret = false;
    float t;
    Point o;
    Vec oa;
    Vec ob;
    Vec oc;

    Vec oa_x_ob;
    Vec ob_x_oc;
    Vec oc_x_oa;

    // ray pass through the plane the triangle lay at
    if (glm::dot(n_a_, ray.d) > 1e-5f || glm::dot(n_a_, ray.d) < -1e-5f) {
        t = glm::dot((a_ - ray.o), n_a_) / (glm::dot(ray.d, n_a_));

        o = ray.At(t);
        oa = a_ - o;
        ob = b_ - o;
        oc = c_ - o;

        oa_x_ob = glm::cross(oa, ob);
        ob_x_oc = glm::cross(ob, oc);
        oc_x_oa = glm::cross(oc, oa);

        // check whether intersection point inside triangle
        if ((glm::dot(oa_x_ob, ob_x_oc) < (1.0f+ 1e-5f) &&  glm::dot(oa_x_ob, ob_x_oc) > (1.0f - 1e-5f)) && (glm::dot(ob_x_oc, oc_x_oa) < (1.0f + 1e-5f) && glm::dot(ob_x_oc, oc_x_oa) > (1.0f - 1e-5f))) {
            // yes, the point is inside
            ret = true;
        }

    }



    if (ret) {
        // hit_record->... = ...
        // hit_record->... = ...
        // ....
        Vec normal;
        Point intersectPoint = ray.At(t);

        if (phong_interpolation_) {
            // hit_record->normal = ...

            Vec ab = b_ - a_;
            Vec ac = c_ - a_;
            Vec ba = a_ - b_;
            Vec bc = c_ - b_;
            Vec ca = a_ - c_;
            Vec cb = b_ - c_;

            Vec ao = o - a_;
            Vec bo = o - b_;
            Vec co = o - c_;

            float tri_abc = glm::length(glm::cross(ab, ac)) / 2.f;
            // the top vertex is b, corresponds to the base triangle
            float tri_b = glm::length(glm::cross(ac, ao)) / 2.f;
            // the top vertex is c, corresponds to the base triangle
            float tri_c = glm::length(glm::cross(ba, bo)) / 2.f;

            // alpha2 corresponds to b
            float alpha_2 = tri_b / tri_abc;
            float alpha_3 = tri_c / tri_abc;
            float alpha_1 = 1.f - alpha_2 - alpha_3;

            normal = glm::normalize(alpha_1 * n_a_ + alpha_2 * n_b_ + alpha_3 * n_c_);
        }
        else {
            // hit_record->normal = ...
            normal = glm::normalize(oa_x_ob);
        }

        hit_record->position = intersectPoint;
        hit_record->normal = normal;
        hit_record->distance = std::sqrt(std::pow(intersectPoint.x - ray.o.x, 2.f) + std::pow(intersectPoint.y - ray.o.y, 2.f) + std::pow(intersectPoint.z - ray.o.z, 2.f));
        hit_record->in_direction = ray.d;
        hit_record->reflection = glm::normalize(ray.d - 2.0f * glm::dot(ray.d, normal) * normal);
        // no need to set material in this function
    }
    return ret;
}

// ---------------------------------------------------------------------------------------------
// ------------------------------ no need to change --------------------------------------------
// ---------------------------------------------------------------------------------------------

// CompleteTriangle
bool CompleteTriangle::Hit(const Ray& ray, HitRecord* hit_record) const {
    bool ret = triangle_.Hit(ray, hit_record);
    if (ret) {
        hit_record->material = material_;
    }
    return ret;
}


// Mesh
Mesh::Mesh(const std::string& file_path,
    const Material& material,
    bool phong_interpolation) :
    ply_data_(file_path), material_(material), phong_interpolation_(phong_interpolation) {
    std::vector<std::array<double, 3>> v_pos = ply_data_.getVertexPositions();
    vertices_.resize(v_pos.size());

    for (int i = 0; i < vertices_.size(); i++) {
        vertices_[i] = Point(v_pos[i][0], v_pos[i][1], v_pos[i][2]);
    }

    f_ind_ = ply_data_.getFaceIndices();

    // Calc face normals
    for (const auto& face : f_ind_) {
        Vec normal = glm::normalize(glm::cross(vertices_[face[1]] - vertices_[face[0]], vertices_[face[2]] - vertices_[face[0]]));
        face_normals_.emplace_back(normal);
    }

    // Calc vertex normals
    vertex_normals_.resize(vertices_.size(), Vec(0.f, 0.f, 0.f));
    for (int i = 0; i < f_ind_.size(); i++) {
        for (int j = 0; j < 3; j++) {
            vertex_normals_[f_ind_[i][j]] += face_normals_[i];
        }
    }
    for (auto& vertex_normal : vertex_normals_) {
        vertex_normal = glm::normalize(vertex_normal);
    }

    // Construct hittable triangles
    for (const auto& face : f_ind_) {
        triangles_.emplace_back(vertices_[face[0]], vertices_[face[1]], vertices_[face[2]],
            vertex_normals_[face[0]], vertex_normals_[face[1]], vertex_normals_[face[2]],
            phong_interpolation_);
    }

    // Calc bounding box
    Point bbox_min(1e5f, 1e5f, 1e5f);
    Point bbox_max(-1e5f, -1e5f, -1e5f);
    for (const auto& vertex : vertices_) {
        bbox_min = glm::min(bbox_min, vertex - 1e-3f);
        bbox_max = glm::max(bbox_max, vertex + 1e-3f);
    }

    // Build Octree
    tree_nodes_.emplace_back(new OctreeNode());
    tree_nodes_.front()->bbox_min = bbox_min;
    tree_nodes_.front()->bbox_max = bbox_max;

    root_ = tree_nodes_.front().get();
    for (int i = 0; i < f_ind_.size(); i++) {
        InsertFace(root_, i);
    }
}

bool Mesh::Hit(const Ray& ray, HitRecord* hit_record) const {
    const bool brute_force = false;
    if (brute_force) {
        // Naive hit algorithm
        float min_dist = 1e5f;
        for (const auto& triangle : triangles_) {
            HitRecord curr_hit_record;
            if (triangle.Hit(ray, &curr_hit_record)) {
                if (curr_hit_record.distance < min_dist) {
                    *hit_record = curr_hit_record;
                    min_dist = curr_hit_record.distance;
                }
            }
        }
        if (min_dist + 1.0 < 1e5f) {
            hit_record->material = material_;
            return true;
        }
        return false;
    }
    else {
        bool ret = OctreeHit(root_, ray, hit_record);
        if (ret) {
            hit_record->material = material_;
        }
        return ret;
    }
}

bool Mesh::IsFaceInsideBox(const std::vector<size_t>& face, const Point& bbox_min, const Point& bbox_max) const {
    for (size_t idx : face) {
        const auto& pt = vertices_[idx];
        for (int i = 0; i < 3; i++) {
            if (pt[i] < bbox_min[i] + 1e-6f) return false;
            if (pt[i] > bbox_max[i] - 1e-6f) return false;
        }
    }
    return true;
}

bool Mesh::IsRayIntersectBox(const Ray& ray, const Point& bbox_min, const Point& bbox_max) const {
    float t_min = -1e5f;
    float t_max = 1e5f;

    for (int i = 0; i < 3; i++) {
        if (glm::abs(ray.d[i]) < 1e-6f) {
            if (ray.o[i] < bbox_min[i] + 1e-6f || ray.o[i] > bbox_max[i] - 1e-6f) {
                t_min = 1e5f;
                t_max = -1e5f;
            }
        }
        else {
            if (ray.d[i] > 0.f) {
                t_min = glm::max(t_min, (bbox_min[i] - ray.o[i]) / ray.d[i]);
                t_max = glm::min(t_max, (bbox_max[i] - ray.o[i]) / ray.d[i]);
            }
            else {
                t_min = glm::max(t_min, (bbox_max[i] - ray.o[i]) / ray.d[i]);
                t_max = glm::min(t_max, (bbox_min[i] - ray.o[i]) / ray.d[i]);
            }
        }
    }

    return t_min + 1e-6f < t_max;
}

void Mesh::InsertFace(OctreeNode* u, size_t face_idx) {
    const Point& bbox_min = u->bbox_min;
    const Point& bbox_max = u->bbox_max;

    Vec bias = bbox_max - bbox_min;
    Vec half_bias = bias * 0.5f;

    bool inside_childs = false;

    for (size_t a = 0; a < 2; a++) {
        for (size_t b = 0; b < 2; b++) {
            for (size_t c = 0; c < 2; c++) {
                size_t child_idx = ((a << 2) | (b << 1) | c);
                Point curr_bbox_min = bbox_min + half_bias * Vec(float(a), float(b), float(c));
                Point curr_bbox_max = curr_bbox_min + half_bias;
                if (IsFaceInsideBox(f_ind_[face_idx], curr_bbox_min, curr_bbox_max)) {
                    if (u->childs[child_idx] == nullptr) {
                        tree_nodes_.emplace_back(new OctreeNode());
                        OctreeNode* child = tree_nodes_.back().get();
                        u->childs[child_idx] = tree_nodes_.back().get();
                        child->bbox_min = curr_bbox_min;
                        child->bbox_max = curr_bbox_max;
                    }
                    InsertFace(u->childs[child_idx], face_idx);
                    inside_childs = true;
                }
            }
        }
    }

    if (!inside_childs) {
        u->face_index.push_back(face_idx);
    }
}

bool Mesh::OctreeHit(OctreeNode* u, const Ray& ray, HitRecord* hit_record) const {
    if (!IsRayIntersectBox(ray, u->bbox_min, u->bbox_max)) {
        return false;
    }
    float distance = 1e5f;
    for (const auto& face_idx : u->face_index) {
        HitRecord curr_hit_record;
        if (triangles_[face_idx].Hit(ray, &curr_hit_record)) {
            if (curr_hit_record.distance < distance) {
                distance = curr_hit_record.distance;
                *hit_record = curr_hit_record;
            }
        }
    }

    for (const auto& child : u->childs) {
        if (child == nullptr) {
            continue;
        }
        HitRecord curr_hit_record;
        if (OctreeHit(child, ray, &curr_hit_record)) {
            if (curr_hit_record.distance < distance) {
                distance = curr_hit_record.distance;
                *hit_record = curr_hit_record;
            }
        }
    }
    return distance + 1 < 1e5f;
}


// Hittable list
void HittableList::PushHittable(const Hittable& hittable) {
    hittable_list_.push_back(&hittable);
}

bool HittableList::Hit(const Ray& ray, HitRecord* hit_record) const {
    float min_dist = 1e5f;
    for (const auto& hittable : hittable_list_) {
        HitRecord curr_hit_record;
        if (hittable->Hit(ray, &curr_hit_record)) {
            if (curr_hit_record.distance < min_dist) {
                *hit_record = curr_hit_record;
                min_dist = curr_hit_record.distance;
            }
        }
    }
    return min_dist + 1.0 < 1e4f;
}