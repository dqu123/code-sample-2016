#ifndef OBJ_PARSER_H
#define OBJ_PARSER_H

#include <vector>
#include <string>

#include "transform.h"

struct Vertex {
    double x;
    double y;
    double z;
    double w;
};

// Compatible with the GL demo.
struct VertexGL { 
    float x;
    float y;
    float z;
};

// 3-D distance between points.
double square_dist(Vertex a, Vertex b);
// Convert vertex to 3-D vector.
Vector3d to_vector(Vertex *v);

VertexGL to_gl_vertex(Vertex v);

struct Face {
    // Vertices
    int v1;
    int v2;
    int v3;
    // Surface normal
    int n1;
    int n2;
    int n3;
};

// Parser for the .obj file format.
class ObjParser {
    public:
        ObjParser() : vertices(), normals(), faces() {};
        ObjParser(const std::string& file_str, const std::string& data_path = "");
        ~ObjParser();
        void print(); 
        void update(const MatrixXd& m);
        void transform(Transformation *t);
        ObjParser *new_copy();

        // OpenGL
        void *get_gl_vertices();
        void *get_gl_normals();
        int get_gl_buffer_size();

        std::vector<Vertex*> vertices;
        std::vector<Vertex*> normals;
        std::vector<Face*> faces;

    private:
        std::vector<VertexGL> gl_vertices;
        std::vector<VertexGL> gl_normals;
};

#endif // OBJ_PARSER_H
