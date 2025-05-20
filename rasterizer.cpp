#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>

using namespace std;

class Color {
public:
    int r, g, b;
    Color(int r = 0, int g = 0, int b = 0) : r(r), g(g), b(b) {}
};

class Point {
public:
    float x, y;
    Point(float x = 0, float y = 0) : x(x), y(y) {}
};

class Canvas {
private:
    int width, height;
    vector<vector<Color>> pixels;
    vector<vector<float>> depthBuffer;
    Color background;
    const float INITIAL_DEPTH = numeric_limits<float>::max();
    
public:
    Canvas(int w, int h, Color bg = Color(0, 0, 0)) : width(w), height(h), background(bg) {
        pixels.resize(height, vector<Color>(width, background));
        depthBuffer.resize(height, vector<float>(width, INITIAL_DEPTH));
    }
    
    void setPixel(int x, int y, float z, const Color& color) {
        if (x >= 0 && x < width && y >= 0 && y < height) {
            if (z < depthBuffer[y][x]) {
                pixels[y][x] = color;
                depthBuffer[y][x] = z;
            }
        }
    }
    
    void clear() {
        for (auto& row : pixels) {
            fill(row.begin(), row.end(), background);
        }
        for (auto& row : depthBuffer) {
            fill(row.begin(), row.end(), INITIAL_DEPTH);
        }
    }
    
    void writePPM(const string& filename) {
        ofstream file(filename);
        file << "P3\n" << width << " " << height << "\n255\n";
        
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                file << pixels[y][x].r << " " 
                     << pixels[y][x].g << " " 
                     << pixels[y][x].b << " ";
            }
            file << "\n";
        }
        file.close();
    }
};

vector<float> Interpolate(int i0, float d0, int i1, float d1) {
    vector<float> values;
    if (i0 == i1) {
        values.push_back(d0);
        return values;
    }
    float a = (d1 - d0) / (i1 - i0);
    float d = d0;

    for (int i = i0; i <= i1; i++) {
        values.push_back(d);
        d += a;
    }
    return values;
}

void DrawLine(Canvas& canvas, const Point& P0, float h0, const Point& P1, float h1, const Color& color) {
    int x0 = static_cast<int>(P0.x);
    int y0 = static_cast<int>(P0.y);
    int x1 = static_cast<int>(P1.x);
    int y1 = static_cast<int>(P1.y);

    if (abs(x1 - x0) > abs(y1 - y0)) {
        if (x0 > x1) {
            swap(x0, x1);
            swap(y0, y1);
            swap(h0, h1);
        }

        auto ys = Interpolate(x0, static_cast<float>(y0), x1, static_cast<float>(y1));
        auto hs = Interpolate(x0, h0, x1, h1);
        for (int x = x0; x <= x1; x++) {
            canvas.setPixel(x, static_cast<int>(ys[x - x0]), hs[x - x0], color);
        }
    } else {
        if (y0 > y1) {
            swap(x0, x1);
            swap(y0, y1);
            swap(h0, h1);
        }

        auto xs = Interpolate(y0, static_cast<float>(x0), y1, static_cast<float>(x1));
        auto hs = Interpolate(y0, h0, y1, h1);
        for (int y = y0; y <= y1; y++) {
            canvas.setPixel(static_cast<int>(xs[y - y0]), y, hs[y - y0], color);
        }
    }
}

void DrawTriangle(Canvas& canvas, 
                 const Point& P0, float h0,
                 const Point& P1, float h1,
                 const Point& P2, float h2,
                 const Color& color) {
    DrawLine(canvas, P0, h0, P1, h1, color);
    DrawLine(canvas, P1, h1, P2, h2, color);
    DrawLine(canvas, P2, h2, P0, h0, color);
}

void DrawFilledTriangle(Canvas& canvas, 
                const Point& P0, float h0, const Point& P1, float h1,
                const Point& P2, float h2, const Color& color) {

    int x0 = static_cast<int>(P0.x);
    int y0 = static_cast<int>(P0.y);
    int x1 = static_cast<int>(P1.x);
    int y1 = static_cast<int>(P1.y);
    int x2 = static_cast<int>(P2.x);
    int y2 = static_cast<int>(P2.y);

    // Sort vertices by y
    if (y1 < y0) { swap(x0, x1); swap(y0, y1); swap(h0, h1); }
    if (y2 < y0) { swap(x0, x2); swap(y0, y2); swap(h0, h2); }
    if (y2 < y1) { swap(x1, x2); swap(y1, y2); swap(h1, h2); }

    // Compute x and h values for edges
    auto x01 = Interpolate(y0, static_cast<float>(x0), y1, static_cast<float>(x1));
    auto h01 = Interpolate(y0, h0, y1, h1);
    
    auto x12 = Interpolate(y1, static_cast<float>(x1), y2, static_cast<float>(x2));
    auto h12 = Interpolate(y1, h1, y2, h2);
    
    auto x02 = Interpolate(y0, static_cast<float>(x0), y2, static_cast<float>(x2));
    auto h02 = Interpolate(y0, h0, y2, h2);

    // Concatenate the short sides
    x01.pop_back();
    auto x012 = x01;
    x012.insert(x012.end(), x12.begin(), x12.end());
    
    h01.pop_back();
    auto h012 = h01;
    h012.insert(h012.end(), h12.begin(), h12.end());

    // Determine which is left and which is right
    int m = x012.size() / 2;
    
    vector<float> x_left, x_right, h_left, h_right;
    
    if (x02[m] < x012[m]) {
        x_left = x02;
        h_left = h02;

        x_right = x012;
        h_right = h012;
    } else {
        x_left = x012;
        h_left = h012;

        x_right = x02;
        h_right = h02;
    }

    // Draw the horizontal segments 
    for (int y = y0; y <= y2; y++) {
        int y_idx = y - y0;
        int x_l = static_cast<int>(x_left[y_idx]);
        int x_r = static_cast<int>(x_right[y_idx]);
        
        auto h_segment = Interpolate(x_l, h_left[y_idx], x_r, h_right[y_idx]);
        
        for (int x = x_l; x <= x_r; x++) {
            float z = h_segment[x - x_l];
            canvas.setPixel(x, y, z, color);
        }
    }
}

Point ViewportToCanvas(float x, float y, int Cw, int Ch) {
    const float Vw = 1.0f, Vh = 1.0f;
    return Point(x * Cw / Vw, y * Ch / Vh);
}

Point ProjectVertex(const vector<float>& v, int Cw, int Ch) {
    const float d = 1.0f;
    Point point = ViewportToCanvas(v[0] * d / v[2], v[1] * d / v[2], Cw, Ch);
    return Point(point.x + Cw/2, Ch/2 - point.y);
}

struct Transform {
    vector<float> translation;
    vector<float> rotation;  // Euler angles (x, y, z)
    vector<float> scale;
    
    Transform(const vector<float>& t = {0,0,0}, 
              const vector<float>& r = {0,0,0}, 
              const vector<float>& s = {1,1,1})
        : translation(t), rotation(r), scale(s) {}
};

vector<float> ApplyTransform(const vector<float>& vertex, const Transform& transform) {
    // Scale first
    vector<float> scaled = {
        vertex[0] * transform.scale[0],
        vertex[1] * transform.scale[1],
        vertex[2] * transform.scale[2]
    };
    
  
    float x = scaled[0];
    float y = scaled[1] * cos(transform.rotation[0]) - scaled[2] * sin(transform.rotation[0]);
    float z = scaled[1] * sin(transform.rotation[0]) + scaled[2] * cos(transform.rotation[0]);
    
    // Y rotation
    float x2 = x * cos(transform.rotation[1]) + z * sin(transform.rotation[1]);
    float y2 = y;
    float z2 = -x * sin(transform.rotation[1]) + z * cos(transform.rotation[1]);
    
    // Z rotation
    float x3 = x2 * cos(transform.rotation[2]) - y2 * sin(transform.rotation[2]);
    float y3 = x2 * sin(transform.rotation[2]) + y2 * cos(transform.rotation[2]);
    float z3 = z2;
    
    // Finally translate
    return {
        x3 + transform.translation[0],
        y3 + transform.translation[1],
        z3 + transform.translation[2]
    };
}

struct Plane {
    float a, b, c, d;
    
    Plane(float a = 0, float b = 0, float c = 0, float d = 0) 
        : a(a), b(b), c(c), d(d) {}
};

struct BoundingSphere {
    vector<float> center;
    float radius;
    
    BoundingSphere(const vector<float>& c = {0,0,0}, float r = 1.0f) 
        : center(c), radius(r) {}
};

struct Triangle {
    vector<int> vertices;
    int color_index;
    
    Triangle(const vector<int>& v, int ci) : vertices(v), color_index(ci) {}
};

struct Model {
    vector<vector<float>> vertices;
    vector<Triangle> triangles;
};

struct Camera {
    vector<float> position = {0, 0, 0};
    vector<float> orientation = {0, 0, 0};
};

struct Instance {
    Model model;
    Transform transform;
    BoundingSphere bounding_sphere;
    
    Instance(const Model& m = Model(), 
             const Transform& t = Transform(),
             const BoundingSphere& bs = BoundingSphere())
        : model(m), transform(t), bounding_sphere(bs) {}
};

struct Scene {
    vector<Instance> instances;
};



float SignedDistance(const Plane& plane, const vector<float>& vertex) {
    return plane.a * vertex[0] + 
           plane.b * vertex[1] + 
           plane.c * vertex[2] + 
           plane.d;
}

vector<float> IntersectLinePlane(const vector<float>& v0, const vector<float>& v1, const Plane& plane) {
    float d0 = SignedDistance(plane, v0);
    float d1 = SignedDistance(plane, v1);
    float t = d0 / (d0 - d1);
    
    return {
        v0[0] + t * (v1[0] - v0[0]),
        v0[1] + t * (v1[1] - v0[1]),
        v0[2] + t * (v1[2] - v0[2])
    };
}

vector<Triangle> ClipTriangle(const Triangle& triangle, 
                             const vector<vector<float>>& vertices, 
                             const Plane& plane) {
    const auto& v0 = vertices[triangle.vertices[0]];
    const auto& v1 = vertices[triangle.vertices[1]];
    const auto& v2 = vertices[triangle.vertices[2]];
    
    float d0 = SignedDistance(plane, v0);
    float d1 = SignedDistance(plane, v1);
    float d2 = SignedDistance(plane, v2);
    
    bool inside0 = d0 >= 0;
    bool inside1 = d1 >= 0;
    bool inside2 = d2 >= 0;
    
    int inside_count = (inside0 ? 1 : 0) + (inside1 ? 1 : 0) + (inside2 ? 1 : 0);
    
    if (inside_count == 3) {
        return {triangle};
    } else if (inside_count == 0) {
        return {};
    } else if (inside_count == 1) {
        int inside_idx = inside0 ? 0 : (inside1 ? 1 : 2);
        int outside_idx1 = (inside_idx + 1) % 3;
        int outside_idx2 = (inside_idx + 2) % 3;
        
        const auto& A = vertices[triangle.vertices[inside_idx]];
        const auto& B = vertices[triangle.vertices[outside_idx1]];
        const auto& C = vertices[triangle.vertices[outside_idx2]];
        
        auto B_prime = IntersectLinePlane(A, B, plane);
        auto C_prime = IntersectLinePlane(A, C, plane);
        
        vector<int> new_vertices = {
            triangle.vertices[inside_idx],
            static_cast<int>(vertices.size()),
            static_cast<int>(vertices.size() + 1)
        };
        
        Triangle new_tri(new_vertices, triangle.color_index);
        return {new_tri};
    } else {
        int outside_idx = !inside0 ? 0 : (!inside1 ? 1 : 2);
        int inside_idx1 = (outside_idx + 1) % 3;
        int inside_idx2 = (outside_idx + 2) % 3;
        
        const auto& A = vertices[triangle.vertices[inside_idx1]];
        const auto& B = vertices[triangle.vertices[inside_idx2]];
        const auto& C = vertices[triangle.vertices[outside_idx]];
        
        auto A_prime = IntersectLinePlane(A, C, plane);
        auto B_prime = IntersectLinePlane(B, C, plane);
        
        vector<int> new_vertices1 = {
            triangle.vertices[inside_idx1],
            triangle.vertices[inside_idx2],
            static_cast<int>(vertices.size())
        };
        
        vector<int> new_vertices2 = {
            static_cast<int>(vertices.size()),
            triangle.vertices[inside_idx2],
            static_cast<int>(vertices.size() + 1)
        };
        
        Triangle new_tri1(new_vertices1, triangle.color_index);
        Triangle new_tri2(new_vertices2, triangle.color_index);
        return {new_tri1, new_tri2};
    }
}

Instance ClipTrianglesAgainstPlane(const Instance& instance, const Plane& plane) {
    Instance clipped_instance = instance;
    clipped_instance.model.triangles.clear();
    
    for (const auto& triangle : instance.model.triangles) {
        auto clipped_tris = ClipTriangle(triangle, instance.model.vertices, plane);
        for (const auto& tri : clipped_tris) {
            clipped_instance.model.triangles.push_back(tri);
        }
    }
    
    return clipped_instance;
}

Instance* ClipInstanceAgainstPlane(Instance& instance, const Plane& plane) {
    float d = SignedDistance(plane, instance.bounding_sphere.center);
    float r = instance.bounding_sphere.radius;
    
    if (d > r) {
        return &instance;
    } else if (d < -r) {
        return nullptr;
    } else {
        Instance* clipped_instance = new Instance(ClipTrianglesAgainstPlane(instance, plane));
        return clipped_instance;
    }
}

Instance* ClipInstance(Instance& instance, const vector<Plane>& planes) {
    Instance* current = &instance;
    
    for (const auto& plane : planes) {
        current = ClipInstanceAgainstPlane(*current, plane);
        if (current == nullptr) {
            return nullptr;
        }
    }
    
    return current;
}

Scene ClipScene(const Scene& scene, const vector<Plane>& planes) {
    Scene clipped_scene;
    
    for (const auto& instance : scene.instances) {
        Instance inst_copy = instance;
        Instance* clipped_instance = ClipInstance(inst_copy, planes);
        if (clipped_instance != nullptr) {
            clipped_scene.instances.push_back(*clipped_instance);
            delete clipped_instance;
        }
    }
    
    return clipped_scene;
}

void RenderTriangle(Canvas& canvas, 
                   const Triangle& triangle, 
                   const vector<vector<float>>& transformedVertices,
                   const vector<Point>& projected) {

    const auto& v0 = transformedVertices[triangle.vertices[0]];
    const auto& v1 = transformedVertices[triangle.vertices[1]];
    const auto& v2 = transformedVertices[triangle.vertices[2]];

    Color color;
    switch(triangle.color_index) {
        case 0: color = Color(0, 0, 0); break;
        case 1: color = Color(0, 0, 255); break;
        case 2: color = Color(255, 0, 0); break;
        case 3: color = Color(0, 255, 0); break;
        case 4: color = Color(128, 0, 128); break;
        case 5: color = Color(255, 255, 0); break;
        case 6: color = Color(0, 255, 255); break;
        default: color = Color(255, 255, 255); break;
    }


    DrawFilledTriangle(canvas, 
                         projected[triangle.vertices[0]], v0[2],
                         projected[triangle.vertices[1]], v1[2],
                         projected[triangle.vertices[2]], v2[2],
                         color);
    
}

void RenderInstance(Canvas& canvas, 
                  const Instance& instance, 
                  const Camera& camera,
                  int canvasWidth, int canvasHeight,
                  bool filled = true) {
    vector<Point> projected;
    vector<vector<float>> transformedVertices;
    const Model& model = instance.model;
    
    
    Transform cameraTransform;
    cameraTransform.translation = {
        -camera.position[0],
        -camera.position[1],
        -camera.position[2]
    };
    
    for (const auto& vertex : model.vertices) {
       
        auto transformed = ApplyTransform(vertex, instance.transform);
        transformed = ApplyTransform(transformed, cameraTransform);
        
        transformedVertices.push_back(transformed);
        projected.push_back(ProjectVertex(transformed, canvasWidth, canvasHeight));
    }
    
    for (const auto& triangle : model.triangles) {
        RenderTriangle(canvas, triangle, transformedVertices, projected);
    }
}

void RenderScene(Canvas& canvas, const Scene& scene, const Camera& camera, int canvasWidth, int canvasHeight, bool filled = true) {
    canvas.clear();
    
    vector<Plane> clipping_planes = {
        Plane(0, 0, 1, 1),    // Near plane
        Plane(0, 0, -1, 10),   // Far plane
        Plane(1, 0, 0, 5),     // Right plane
        Plane(-1, 0, 0, 5),    // Left plane
        Plane(0, 1, 0, 5),     // Top plane
        Plane(0, -1, 0, 5)     // Bottom plane
    };
    
    Scene clipped_scene = ClipScene(scene, clipping_planes);
    
    for (const auto& instance : clipped_scene.instances) {
        RenderInstance(canvas, instance, camera, canvasWidth, canvasHeight, filled);
    }
}

int main() {
    const int Cw = 800;
    const int Ch = 600;
    
    Canvas canvas(Cw, Ch, Color(255, 255, 255));
    
    // Define cube model
    Model cube;
    cube.vertices = {
        { 1,  1,  1}, {-1,  1,  1}, {-1, -1,  1}, { 1, -1,  1},
        { 1,  1, -1}, {-1,  1, -1}, {-1, -1, -1}, { 1, -1, -1}
    };
    
    cube.triangles = {
        Triangle({0, 1, 2}, 2), Triangle({0, 2, 3}, 2),
        Triangle({4, 0, 3}, 3), Triangle({4, 3, 7}, 3),
        Triangle({5, 4, 7}, 1), Triangle({5, 7, 6}, 1),
        Triangle({1, 5, 6}, 5), Triangle({1, 6, 2}, 5),
        Triangle({4, 5, 1}, 4), Triangle({4, 1, 0}, 4),
        Triangle({2, 6, 7}, 6), Triangle({2, 7, 3}, 6)
    };

    // Create scene with instances
    Scene scene;
    
    // Original cube (centered)
    Instance originalCube;
    originalCube.model = cube;
    originalCube.transform = Transform({0, 0, 10.0f});
    originalCube.bounding_sphere = BoundingSphere({0, 0, 8}, sqrt(3));
    scene.instances.push_back(originalCube);
    
    // Translated cube (left)
    Instance translatedCube;
    translatedCube.model = cube;
    translatedCube.transform = Transform({-3.0f, 0.0f, 10.0f});
    translatedCube.bounding_sphere = BoundingSphere({-3, -1, 7}, sqrt(3));
    scene.instances.push_back(translatedCube);
    
    // Scaled cube (right)
    Instance scaledCube;
    scaledCube.model = cube;
    scaledCube.transform = Transform({3.0f, 0.0f, 10.0f}, {0, 0, 0}, {0.5f, 0.5f, 0.5f});
    scaledCube.bounding_sphere = BoundingSphere({3, 1, 9}, sqrt(3) * 0.5f);
    scene.instances.push_back(scaledCube);
    

    // Set up camera
    Camera camera;
    camera.position = {0, 0, 0};

    // Render scene with filled triangles
    RenderScene(canvas, scene, camera, Cw, Ch, true);
    
    canvas.writePPM("rasterizer_output.ppm");
    
    return 0;
}
