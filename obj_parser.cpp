#include <iostream>
#include <fstream>
#include <sstream>

#include <Eigen/Dense>

#include "obj_parser.h"

using namespace Eigen;

double square_dist(Vertex a, Vertex b) {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    double dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
}

// Convert to Vector3d.
Vector3d to_vector(Vertex *v) {
    return Vector3d(v->x, v->y, v->z);
}

// Construct ObjParser from the string name of a file.
ObjParser::ObjParser(const std::string& file_str, const std::string& data_path) {
    this->vertices.push_back(NULL); // 1-index
    this->normals.push_back(NULL); // 1-index

    std::ifstream file(data_path + file_str, std::ifstream::in);
    std::string line;

    if (file.is_open()) {
        // Process file line by line.
        while (std::getline(file, line)) {
            std::stringstream ss(line);
            std::string buf;

            ss >> buf;

            // Determine line type based on first space delimited token.
            if (buf == "v") {
                Vertex *v = new Vertex;
                ss >> buf;
                v->x = std::stod(buf);
                ss >> buf;
                v->y = std::stod(buf);
                ss >> buf;
                v->z = std::stod(buf);
                v->w = 1;
                this->vertices.push_back(v);
            }
            else if (buf == "vn") {
                Vertex *v = new Vertex;
                ss >> buf;
                v->x = std::stod(buf);
                ss >> buf;
                v->y = std::stod(buf);
                ss >> buf;
                v->z = std::stod(buf);
                v->w = 1;
                this->normals.push_back(v);
            }
            else if (buf == "f") {
                Face *f = new Face;
                ss >> buf;
                // Location of "//" divider
                std::size_t pos = buf.find("/");
                f->v1 = std::stoi(buf.substr(0, pos));
                f->n1 = std::stoi(buf.substr(pos + 2));
                ss >> buf;
                pos = buf.find("/");
                f->v2 = std::stoi(buf.substr(0, pos));
                f->n2 = std::stoi(buf.substr(pos + 2));
                ss >> buf;
                pos = buf.find("/");
                f->v3 = std::stoi(buf.substr(0, pos));
                f->n3 = std::stoi(buf.substr(pos + 2));
                this->faces.push_back(f);
            }
        }
    }

    file.close();

}

// Clean up vertices and faces.
ObjParser::~ObjParser() {
    for (int i = 0; i < vertices.size(); i++) {
        delete vertices[i];
    }

    for (int i = 0; i < normals.size(); i++) {
        delete normals[i];
    }

    for (int i = 0; i < faces.size(); i++) {
        delete faces[i];
    }
}

// Outputs contents of .obj file.
void ObjParser::print() {
    // Print vertices in order.
    for (int i = 1; i < vertices.size(); i++) {
        if (vertices[i] == NULL) { 
            std::cout << "NULL vertex\n";
        }
        else {
            std::cout << "v " << vertices[i]->x << " ";
            std::cout << vertices[i]->y << " " << vertices[i]->z << "\n";
        }
    }

    // Print normals in order.
    for (int i = 1; i < normals.size(); i++) {
        if (normals[i] == NULL) { 
            std::cout << "NULL normal\n";
        }
        else {
            std::cout << "vn " << normals[i]->x << " ";
            std::cout << normals[i]->y << " " << normals[i]->z << "\n";
        }
    }

    // Print faces in order.
    for (int i = 0; i < faces.size(); i++) {
        if (faces[i] == NULL) {
            std::cout << "NULL face\n";
        }
        else {
            std::cout << "f " << faces[i]->v1 << "//" << faces[i]->n1 << " "
                      << faces[i]->v2 << "//" << faces[i]->n2 << " "
                      << faces[i]->v3 << "//" << faces[i]->n3 << "\n";
        }
    }
}

// Applys matrix to the vertices
void ObjParser::update(const MatrixXd& m) {

    // Update each vertice using the matrix.
    for (int i = 1; i < vertices.size(); i++) {
        VectorXd temp(4);
        temp << vertices[i]->x, vertices[i]->y, vertices[i]->z, 1;
        temp = m * temp;
        vertices[i]->x = temp(0);
        vertices[i]->y = temp(1);
        vertices[i]->z = temp(2);
        vertices[i]->w = temp(3);
    }
    
}

// Transforms vertices and normals.
void ObjParser::transform(Transformation *t) {
    update(t->m);
    // Transform each normal using the matrix.
    for (int i = 1; i < normals.size(); i++) {
        VectorXd temp(4);
        temp << normals[i]->x, normals[i]->y, normals[i]->z, 1;
        temp = t->n_transform() * temp;

        double x = temp(0);
        double y = temp(1);
        double z = temp(2);
        double w = temp(3);
        normals[i]->x = x;
        normals[i]->y = y;
        normals[i]->z = z;
        normals[i]->w = w;

        // Normalize vector
        double norm = sqrt(x * x + y * y + z * z);

        normals[i]->x /= norm;
        normals[i]->y /= norm;
        normals[i]->z /= norm;
    }

}



// Allocate a new copy of a ObjParser with the same values for
// vertices and faces. Caller is responsible for freeing used
// memory.
ObjParser *ObjParser::new_copy() {
    ObjParser *op = new ObjParser();

    // Copy vertices.
    op->vertices.push_back(NULL);
    for (int i = 1; i < this->vertices.size(); i++) {
        Vertex *v = new Vertex();
        v->x = this->vertices[i]->x;
        v->y = this->vertices[i]->y;
        v->z = this->vertices[i]->z;
        v->w = this->vertices[i]->w;
        op->vertices.push_back(v);
    }

    // Copy normals.
    op->normals.push_back(NULL);
    for (int i = 1; i < this->normals.size(); i++) {
        Vertex *v = new Vertex();
        v->x = this->normals[i]->x;
        v->y = this->normals[i]->y;
        v->z = this->normals[i]->z;
        v->w = this->normals[i]->w;
        op->normals.push_back(v);
    }

    // Copy faces.
    for (int i = 0; i < this->faces.size(); i++) {
        Face *f = new Face();
        f->v1 = this->faces[i]->v1;
        f->v2 = this->faces[i]->v2;
        f->v3 = this->faces[i]->v3;
        f->n1 = this->faces[i]->n1;
        f->n2 = this->faces[i]->n2;
        f->n3 = this->faces[i]->n3;
        op->faces.push_back(f);
    }

    return op;
}


VertexGL to_gl_vertex(Vertex v) { 
    VertexGL w;
    w.x = v.x;
    w.y = v.y;
    w.z = v.z;

    //std::cerr << "Next vertex: ";
    //std::cerr << w.x << " " << w.y << " " << w.z << "\n";
    return w;
}

void *ObjParser::get_gl_vertices() {
    gl_vertices = std::vector<VertexGL>();

    //std::cerr << "Printing vertices:\n";
    for (int i = 0; i < faces.size(); i++) {
        Face *f = faces[i];
        gl_vertices.push_back(to_gl_vertex(*vertices[f->v1])); 
        gl_vertices.push_back(to_gl_vertex(*vertices[f->v2])); 
        gl_vertices.push_back(to_gl_vertex(*vertices[f->v3])); 
    }
    return &gl_vertices[0];
}

void *ObjParser::get_gl_normals() {
    gl_normals = std::vector<VertexGL>();

    //std::cerr << "Printing normals:\n";
    for (int i = 0; i < faces.size(); i++) {
        Face *f = faces[i];
        gl_normals.push_back(to_gl_vertex(*normals[f->n1])); 
        gl_normals.push_back(to_gl_vertex(*normals[f->n2])); 
        gl_normals.push_back(to_gl_vertex(*normals[f->n3])); 
    }
    return &gl_normals[0];
}

int ObjParser::get_gl_buffer_size() {
    return gl_vertices.size();
}
