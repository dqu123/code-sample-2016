#ifndef MULT_OBJ_PARSER_H
#define MULT_OBJ_PARSER_H

#include <map>
#include <vector>
#include <string>

#include "obj_parser.h"
#include "transform.h"
#include "lighting.h"

// Camera details.
struct Camera {
    double near;
    double far;
    double left;
    double right;
    double top;
    double bottom;
};

// Contains metadata for a given copy of an object.
class ObjectCopy {
    public:
        ObjectCopy(ObjParser *op, Transformation *t, LightProps *lp, std::string name);
        ~ObjectCopy();
        void print();
        std::string name;
        ObjParser *op;
        Transformation *t;
        LightProps *lp;
};

// Parser for the multi object file format described in the assignment.
class MultObjParser {
    public:
        MultObjParser(const std::string& file_str, const std::string& data_path = "",
                      int xres = 800, int yres = 800);
        ~MultObjParser();
        void print();
        void to_ppm(int xres, int yres);
        bool do_bresenham(int x0, int y0, int x1, int y1, int *grid, int xres, int yres, int color);
        void shade(int shade_mode);
        void debug_camera();

        std::vector<ObjectCopy*> objects;
        std::vector<PointLight*> lights;
        Camera cam;
        Transformation T_c;
        Transformation R_c;
    
    private:
        Vertex *world_to_NDC(Vertex *v);
        Vertex *NDC_to_screen(Vertex *v);
        void raster_colored_triangle(Vertex *a, Vertex *b, Vertex *c, 
                Color& c_a, Color& c_b, Color& c_c);
        void gouraud_shading(LightProps *lp, Vertex *a, Vertex *b, Vertex *c, 
                Vertex *n_a, Vertex *n_b, Vertex *n_c);
        void do_gouraud();
        void phong_shading(LightProps *lp, Vertex *a, Vertex *b, Vertex *c,
                Vertex *n_a, Vertex *n_b, Vertex *n_c);
        void do_phong();

        // Fields
        std::map<std::string, ObjParser*> obj_map;
        std::map<std::string, int> obj_count;
        MatrixXd C_inverse;
        MatrixXd PPM;
        int xres;
        int yres;
        Color **grid;
        double *buffer;
};


#endif // MULT_OBJ_PARSER_H
