#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>

#include "lighting.h"

// Function prototypes.
Color cwise_min(Color a, Color b);
Color cwise_mult(Color a, Color b);
Color normalize(Color a);
// Helper function for barycentric coordinates
double compute_f(int *X, int *Y, int x, int y, int i, int j);

double compute_f(int *X, int *Y, int x, int y, int i, int j) {
    return (Y[i] - Y[j]) * x + (X[j] - X[i]) * y + X[i] * Y[j] - X[j] * Y[i];
}

// Compute barycentric coordinate alpha.
double compute_alpha(int xa, int ya, int xb, int yb, int xc, int yc, int x, int y) {
    int X[] = {xa, xb, xc};
    int Y[] = {ya, yb, yc};

    return compute_f(X, Y, x, y, 1, 2) / compute_f(X, Y, xa, ya, 1, 2);
}

// Compute barycentric coordinate beta.
double compute_beta(int xa, int ya, int xb, int yb, int xc, int yc, int x, int y) {
    int X[] = {xa, xb, xc};
    int Y[] = {ya, yb, yc};

    return compute_f(X, Y, x, y, 0, 2) / compute_f(X, Y, xb, yb, 0, 2);
}


// Compute barycentric coordinate gamma.
double compute_gamma(int xa, int ya, int xb, int yb, int xc, int yc, int x, int y) {
    int X[] = {xa, xb, xc};
    int Y[] = {ya, yb, yc};

    return compute_f(X, Y, x, y, 0, 1) / compute_f(X, Y, xc, yc, 0, 1);
}

// Implement lighting model (Algorithm 3 in the notes).
Color lighting(Vertex *P, Vertex *n, LightProps *lp, 
               std::vector<PointLight*> lights, Vertex *e) {

    Color n_c(n->x, n->y, n->z);

    Color c_d = *(lp->props["diffuse"]);
    Color c_a = *(lp->props["ambient"]);
    Color c_s = *(lp->props["specular"]);
    double phong = lp->phong;

    Color diffuse_sum = Color(0, 0, 0);
    Color specular_sum = Color(0, 0, 0);

    Color e_direction = normalize(Color(e->x, e->y, e->z) - Color(P->x, P->y, P->z));

    for (int i = 0; i < lights.size(); i++) {
        Vertex position = lights[i]->v;
        Color l_p(position.x, position.y, position.z);
        // Account for attenuation.
        double attenuation = 1.0 / (1 + lights[i]->k * square_dist(*P, position));
        Color l_c = lights[i]->c * attenuation;
        
        Color l_direction = normalize(l_p - Color(P->x, P->y, P->z));

        Color l_diffuse = l_c * std::max(0.0, n_c.dot(l_direction));
        diffuse_sum += l_diffuse;
        

        double dp = n_c.dot(normalize(e_direction + l_direction));
        Color l_specular = l_c * std::pow(std::max(0.0, dp), phong);
        specular_sum += l_specular;

    }

    Color c = cwise_min(Color(1, 1, 1), c_a + cwise_mult(diffuse_sum, c_d) +
                        cwise_mult(specular_sum, c_s));
    return c;
}

// Compute component-wise min of two colors.
Color cwise_min(Color a, Color b) {
    return Color(std::min(a(0), b(0)), std::min(a(1), b(1)),
                 std::min(a(2), b(2)));
}

// Compute component-wise multiplication of two colors.
Color cwise_mult(Color a, Color b) {
    return Color(a(0) * b(0), a(1) * b(1), a(2) * b(2));
}

// Normalize the given vector
Color normalize(Color a) {
    Color result(a(0), a(1), a(2));
    double norm = sqrt(a(0) * a(0) + a(1) * a(1) + a(2) * a(2));
    result /= norm;
    return result;
}

LightProps::~LightProps() {
    for (std::map<std::string, Color*>::iterator it = props.begin(); it != props.end(); it++) {
        delete it->second;
    }
}

// Add material light property.
void LightProps::add_prop(const std::string& line) {
    std::stringstream ss(line);
    std::string buf, property;

    ss >> property;

    if (property == "ambient" || property == "diffuse" || property == "specular") {
        double r, g, b;
        ss >> buf;
        r = std::stod(buf);
        ss >> buf;
        g = std::stod(buf);
        ss >> buf;
        b = std::stod(buf);
        props[property] = new Color(r, g, b);
    }
    else if (property == "shininess") {
        ss >> buf;
        phong = std::stod(buf);
    }
}

float *LightProps::get_gl_prop(std::string prop) {
    Color c = *(props[prop]);
    
    if (prop == "ambient") {
        for (int i = 0; i < 3; i++) { 
            ambient[i] = c(i);
        }
        return ambient;
    }
    else if (prop == "diffuse") {
        for (int i = 0; i < 3; i++) { 
            diffuse[i] = c(i);
        }
        return diffuse;
    }
    else if (prop == "specular") { 
        for (int i = 0; i < 3; i++) { 
            specular[i] = c(i);
        }
        return specular;
    }
}

PointLight::PointLight(const std::string& line) {
    std::stringstream ss(line);
    std::string buf;

    // Read in vertex
    Vertex v;
    ss >> buf;
    v.x = std::stod(buf);
    ss >> buf;
    v.y = std::stod(buf);
    ss >> buf;
    v.z = std::stod(buf);
    v.w = 1;
    this->v = v;

    // Read in comma
    ss >> buf;

    // Read in color
    double r, g, b;
    ss >> buf;
    r = std::stod(buf);
    ss >> buf;
    g = std::stod(buf);
    ss >> buf;
    b = std::stod(buf);
    c = Color(r, g, b);

    // Read in comma
    ss >> buf; 

    // Read in attenuation
    ss >> buf;
    k = std::stod(buf);
}


void PointLight::debug() {
    std::cerr << "Position: " << v.x << " "
              << v.y << " " << v.z << " " << v.w << "\n";
    
    std::cerr << "Color:\n" << c << "\n";

    std::cerr << "Attenuation: " << k << "\n";
}

float *PointLight::get_gl_position() {
    gl_position[0] = v.x;
    gl_position[1] = v.y;
    gl_position[2] = v.z;
    gl_position[3] = v.w;
    return gl_position;
}

// Return the color in the array format needed by OpenGL.
float *PointLight::get_gl_color() {
    gl_color[0] = c(0);
    gl_color[1] = c(1);
    gl_color[2] = c(2);
    return gl_color;
}


