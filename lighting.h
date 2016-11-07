#ifndef LIGHTING_H
#define LIGHTING_H

#include <map>
#include <string>

#include "obj_parser.h"

// Represent color with a 3 dimensional vector.
typedef Vector3d Color;

class LightProps {
    public:
        LightProps() {};
        ~LightProps();
        void add_prop(const std::string& line);
        // OpenGL
        float *get_gl_prop(std::string prop);

        std::map<std::string, Color*> props;
        double phong;
    private:
        float ambient[3];
        float diffuse[3];
        float specular[3];
};

class PointLight {
    public:
        PointLight(const std::string& line);
        ~PointLight() {};

        void debug();

        // OpenGL
        float *get_gl_position();
        float *get_gl_color();

        Vertex v;
        Color c;
        double k; // Attenuation parameter.

    private:
        float gl_position[4];
        float gl_color[3];
};

Color lighting(Vertex* P, Vertex* n, LightProps* lp, 
               std::vector<PointLight*> lights, Vertex* e);
double compute_alpha(int xa, int ya, int xb, int yb, int xc, int yc, int x, int y);
double compute_beta(int xa, int ya, int xb, int yb, int xc, int yc, int x, int y);
double compute_gamma(int xa, int ya, int xb, int yb, int xc, int yc, int x, int y);


#endif // LIGHTING_H
