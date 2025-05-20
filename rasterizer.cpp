
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

void DrawLine(Canvas& canvas, const Point& P0, float z0, const Point& P1, float z1, const Color& color) {
    int x0 = static_cast<int>(P0.x);
    int y0 = static_cast<int>(P0.y);
    int x1 = static_cast<int>(P1.x);
    int y1 = static_cast<int>(P1.y);

    if (abs(x1 - x0) > abs(y1 - y0)) {
        if (x0 > x1) {
            swap(x0, x1);
            swap(y0, y1);
            swap(z0, z1);
        }

        auto ys = Interpolate(x0, static_cast<float>(y0), x1, static_cast<float>(y1));
        auto zs = Interpolate(x0, z0, x1, z1);
        for (int x = x0; x <= x1; x++) {
            canvas.setPixel(x, static_cast<int>(ys[x - x0]), zs[x - x0], color);
        }
    } else {
        if (y0 > y1) {
            swap(x0, x1);
            swap(y0, y1);
            swap(z0, z1);
        }

        auto xs = Interpolate(y0, static_cast<float>(x0), y1, static_cast<float>(x1));
        auto zs = Interpolate(y0, z0, y1, z1);
        for (int y = y0; y <= y1; y++) {
            canvas.setPixel(static_cast<int>(xs[y - y0]), y, zs[y - y0], color);
        }
    }
}

void DrawTriangle(Canvas& canvas, 
                 const Point& P0, float z0,
                 const Point& P1, float z1,
                 const Point& P2, float z2,
                 const Color& color) {
    DrawLine(canvas, P0, z0, P1, z1, color);
    DrawLine(canvas, P1, z1, P2, z2, color);
    DrawLine(canvas, P2, z2, P0, z0, color);
}
/*
void DrawFilledTriangle(Canvas& canvas, 
                       const Point& P0, float z0,
                       const Point& P1, float z1,
                       const Point& P2, float z2,
                       const Color& color) {
    int x0 = static_cast<int>(P0.x);
    int y0 = static_cast<int>(P0.y);
    int x1 = static_cast<int>(P1.x);
    int y1 = static_cast<int>(P1.y);
    int x2 = static_cast<int>(P2.x);
    int y2 = static_cast<int>(P2.y);

    // Sort vertices by y
    if (y1 < y0) { swap(x0, x1); swap(y0, y1); swap(z0, z1); }
    if (y2 < y0) { swap(x0, x2); swap(y0, y2); swap(z0, z2); }
    if (y2 < y1) { swap(x1, x2); swap(y1, y2); swap(z1, z2); }

    // Compute x and z values for edges
    auto x01 = Interpolate(y0, static_cast<float>(x0), y1, static_cast<float>(x1));
    auto z01 = Interpolate(y0, z0, y1, z1);
    
    auto x12 = Interpolate(y1, static_cast<float>(x1), y2, static_cast<float>(x2));
    auto z12 = Interpolate(y1, z1, y2, z2);
    
    auto x02 = Interpolate(y0, static_cast<float>(x0), y2, static_cast<float>(x2));
    auto z02 = Interpolate(y0, z0, y2, z2);

    // Concatenate the short sides
    x01.pop_back();
    auto x012 = x01;
    x012.insert(x012.end(), x12.begin(), x12.end());
    
    z01.pop_back();
    auto z012 = z01;
    z012.insert(z012.end(), z12.begin(), z12.end());

    // Determine which is left and which is right
    int m = x012.size() / 2;
    vector<float> x_left, x_right, z_left, z_right;
    
    if (x02[m] < x012[m]) {
        x_left = x02;
        z_left = z02;
        x_right = x012;
        z_right = z012;
    } else {
        x_left = x012;
        z_left = z012;
        x_right = x02;
        z_right = z02;
    }

    // Draw the horizontal segments
    for (int y = y0; y <= y2; y++) {
        int y_idx = y - y0;
        int x_l = static_cast<int>(x_left[y_idx]);
        int x_r = static_cast<int>(x_right[y_idx]);
        
        auto z_segment = Interpolate(x_l, z_left[y_idx], x_r, z_right[y_idx]);
        
        for (int x = x_l; x <= x_r; x++) {
            float z = z_segment[x - x_l];
            canvas.setPixel(x, y, z, color);
        }
    }
}
*/
Point ViewportToCanvas(float x, float y, int canvasWidth, int canvasHeight) {
    const float Vw = 1.0f, Vh = 1.0f;
    return Point(x * canvasWidth / Vw, y * canvasHeight / Vh);
}

Point ProjectVertex(const vector<float>& v, int canvasWidth, int canvasHeight) {
    const float d = 1.0f;
    Point point = ViewportToCanvas(v[0] * d / v[2], v[1] * d / v[2], canvasWidth, canvasHeight);
    return Point(point.x + canvasWidth/2, canvasHeight/2 - point.y);
}

vector<vector<float>> MultiplyMatrices(const vector<vector<float>>& a, const vector<vector<float>>& b) {
    vector<vector<float>> result(4, vector<float>(4, 0));
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            for (int k = 0; k < 4; ++k) {
                result[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return result;
}

vector<vector<float>> MakeTranslationMatrix(float tx, float ty, float tz) {
    return {
        {1, 0, 0, tx},
        {0, 1, 0, ty},
        {0, 0, 1, tz},
        {0, 0, 0, 1}
    };
}

vector<vector<float>> MakeCameraMatrix(const vector<float>& position, const vector<float>& orientation) {
    return {
        {1, 0, 0, -position[0]},
        {0, 1, 0, -position[1]},
        {0, 0, 1, -position[2]},
        {0, 0, 0, 1}
    };
}

vector<float> TransformVertex(const vector<float>& vertex, const vector<vector<float>>& transform) {
    vector<float> v = {vertex[0], vertex[1], vertex[2], 1.0f};
    vector<float> result(4, 0);
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            result[i] += transform[i][j] * v[j];
        }
    }
    return {result[0], result[1], result[2]};
}

struct Triangle {
    vector<int> vertices;
    int color_index;
    
    Triangle(const vector<int>& v, int ci) : vertices(v), color_index(ci) {}
};

struct Model {
    vector<vector<float>> vertices;
    vector<Triangle> triangles;
};

struct Instance {
    Model model;
    vector<vector<float>> transform;
};

struct Camera {
    vector<float> position = {0, 0, 0};
    vector<float> orientation = {0, 0, 0};
};

struct Scene {
    vector<Instance> instances;
};

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

    DrawTriangle(canvas, 
                      projected[triangle.vertices[0]], v0[2],
                      projected[triangle.vertices[1]], v1[2],
                      projected[triangle.vertices[2]], v2[2],
                      color);
}

// ... (previous code remains the same until the RenderModel function)

vector<float> ApplyTransform(const vector<float>& vertex, const vector<vector<float>>& transform) {
    // This function combines scaling, rotation, and translation from the transform matrix
    // Since we're already using matrix multiplication for transformations, we can just use TransformVertex
    return TransformVertex(vertex, transform);
}

void RenderInstance(Canvas& canvas, 
                   const Instance& instance, 
                   const vector<vector<float>>& cameraMatrix,
                   int canvasWidth, int canvasHeight) {
    vector<Point> projected;
    vector<vector<float>> transformedVertices;
    const Model& model = instance.model;
    
    // Combine camera and instance transform
    auto M = MultiplyMatrices(cameraMatrix, instance.transform);
    
    // Transform and project vertices
    for (const auto& vertex : model.vertices) {
        auto transformed = ApplyTransform(vertex, M);
        transformedVertices.push_back(transformed);
        projected.push_back(ProjectVertex(transformed, canvasWidth, canvasHeight));
    }
    
    // Render each triangle with depth information
    for (const auto& triangle : model.triangles) {
        RenderTriangle(canvas, triangle, transformedVertices, projected);
    }
}

void RenderScene(Canvas& canvas, const Scene& scene, const Camera& camera, int canvasWidth, int canvasHeight) {
    canvas.clear();
    auto M_camera = MakeCameraMatrix(camera.position, camera.orientation);
    
    for (const auto& instance : scene.instances) {
        RenderInstance(canvas, instance, M_camera, canvasWidth, canvasHeight);
    }
}



int main() {
    const int Cw = 500;
    const int Ch = 500;
    
    Canvas canvas(Cw, Ch, Color(50, 50, 50));
    
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

    // Create scene with instance
    Scene scene;

    //Original cube instance
    Instance originalCube;
    originalCube.model = cube;
    originalCube.transform = MakeTranslationMatrix(0, 0, 10.0f);
    scene.instances.push_back(originalCube);

    //Translated cube instance
    Instance translatedCube;
    translatedCube.model = cube;
    translatedCube.transform = MakeTranslationMatrix(-3.0f, -1.0, 10.0);
    scene.instances.push_back(translatedCube);

    //Scaled cube instance
    Instance scaledCube;
    scaledCube.model = cube;
    scaledCube.transform = MultiplyMatrices(MakeTranslationMatrix(3.0f, 1.0f, 10.0), 
                                             MakeTranslationMatrix(0.5f, 0.5f, 0.5f));
    scene.instances.push_back(scaledCube);
    
    // Set up camera
    Camera camera;
    camera.position = {0, 0, 0};

    // Render scene
    RenderScene(canvas, scene, camera, Cw, Ch);
    
    canvas.writePPM("output.ppm");
    cout << "Rendered image saved to output.ppm" << endl;
    return 0;
}