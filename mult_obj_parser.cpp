#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>

#include "mult_obj_parser.h"

#define BG_COLOR "0 0 0"
#define COLOR_1 "255 255 255"
#define OUT_OF_NDC 2
#define COLOR_RESOLUTION 255

// Create object copy and apply the given transformation.
ObjectCopy::ObjectCopy(ObjParser *op, Transformation *t, LightProps *lp, std::string name) {
    this->op = op;
    this->t = t;
    this->lp = lp;
    this->name = name;

    // Transform vertices and normals according to t.
    op->transform(t);
}

// Clean up object copy.
ObjectCopy::~ObjectCopy() {
    delete op;
    delete t;
    delete lp;
}

// Print out all vertices for each copy.
void ObjectCopy::print() { 
    std::cout << name << "\n";
    for (int i = 1; i < op->vertices.size(); i++) {
        std::cout << op->vertices[i]->x << " " <<
            op->vertices[i]->y << " " << op->vertices[i]->z << "\n";
    }
}

// Parse the given file.
MultObjParser::MultObjParser(const std::string& file_str, const std::string& data_path,
                             int xres, int yres) : C_inverse(4, 4), PPM(4, 4) {
    this->xres = xres;
    this->yres = yres;

    // Grid representing picture.
    this->grid = new Color*[xres * yres]();
    // Buffer for depth buffering.
    this->buffer = new double[xres * yres]();
    for (int i = 0; i < xres * yres; i++) {
        this->buffer[i] = OUT_OF_NDC;
        this->grid[i] = NULL;
    }

    std::ifstream file(data_path + file_str, std::ifstream::in);

    // Camera transformation
    double near, far, left, right, top, bottom;

    if (file.is_open()) {
        std::string line;

        // Read in camera details
        std::getline(file, line); // camera: 

        while (std::getline(file, line) && line != "objects:") {
            std::stringstream ss(line);
            std::string buf;

            ss >> buf;

            if (buf == "position") {
                // The remaining portion is a translation vector
                std::getline(ss, line); 
                T_c.add("t " + line);

            }
            else if (buf == "orientation") {
                // The remaining portion is a rotation vector
                std::getline(ss, line);
                R_c.add("r " + line);
            }
            else if (buf == "near") {
                ss >> buf;
                near = std::stof(buf);
            }
            else if (buf == "far") {
                ss >> buf;
                far = std::stof(buf);
            }
            else if (buf == "left") {
                ss >> buf;
                left = std::stof(buf);
            }
            else if (buf == "right") {
                ss >> buf;
                right = std::stof(buf);
            }
            else if (buf == "top") {
                ss >> buf;
                top = std::stof(buf);
            }
            else if (buf == "bottom") {
                ss >> buf;
                bottom = std::stof(buf);
            }
            else if (buf == "light") {
                std::getline(ss, buf);
                lights.push_back(new PointLight(buf));
            }
        }

        cam.near = near; 
        cam.far = far;
        cam.left = left; 
        cam.right = right; 
        cam.top = top;
        cam.bottom = bottom;
                    
        // First read object files
        while (std::getline(file, line) && line != "") { 
            std::stringstream ss(line);
            std::string buf;
            std::string obj_file_str;
            
            ss >> buf;
            ss >> obj_file_str;

            obj_map[buf] = new ObjParser(obj_file_str, data_path);
            obj_count[buf] = 0;
        }

        ObjParser *op = NULL;
        Transformation *t = NULL;
        LightProps *lp = NULL;
        std::string obj_name = "";

        // Create copies of objects based on ObjParsers in map
        while (!file.eof()) {

            std::getline(file, obj_name);

            op = obj_map[obj_name]->new_copy();
            t = new Transformation();
            lp = new LightProps();
            
            while (std::getline(file, line)) {
                t->add(line);
                lp->add_prop(line);

                if (line == "" && op != NULL && t != NULL && lp != NULL && obj_name != "") {
                    obj_count[obj_name]++;
                    objects.push_back(new ObjectCopy(op, t, lp,
                                obj_name + "_copy" + std::to_string(obj_count[obj_name])));
                    op = NULL;
                    t = NULL;
                    lp = NULL;
                    obj_name = "";
                    break;
                }
            }
        }

        if (op != NULL && t != NULL && lp != NULL && obj_name != "") { 
            obj_count[obj_name]++;
            objects.push_back(new ObjectCopy(op, t, lp,
                        obj_name + "_copy" + std::to_string(obj_count[obj_name])));
        }
        
        // Compute camera transformation.
        C_inverse = (T_c.m * R_c.m).inverse();

        // Compute Perspective projection matrix
        PPM << 2 * near / (right - left),                         0,  (right + left) / (right - left), 0,
                                       0, 2 * near / (top - bottom),  (top + bottom) / (top - bottom), 0,
                                       0,                         0, -(far + near) / (far - near), -2 * far * near / (far - near),
                                       0,                         0,                                -1, 0;


    }

    file.close();
}

void MultObjParser::shade(int shade_mode) {
    // Shade depending on shade_mode.
    switch (shade_mode) {
    case 0:
        do_gouraud();
        break;
    case 1:
        do_phong();
        break;
    default:
        break;
    }

    // Print in PPM format
    std::cout << "P3\n" << xres << " " << yres << "\n"
              << COLOR_RESOLUTION << "\n";
    for (int j = 0; j < yres; j++) {
        for (int i = 0; i < xres; i++) {
            if (grid[i + j * xres] != NULL) {
                Color c = *(grid[i + j * xres]);
                c *= COLOR_RESOLUTION;
                // Round to integer color values.
                std::cout << std::round(c(0)) << " " << std::round(c(1)) << 
                          " " << std::round(c(2)) << "\n";
            }
            else {
                std::cout << "0 0 0\n";
            }
        }
    }
}

Vertex *MultObjParser::world_to_NDC(Vertex *v) {
    Vertex *new_v = NULL;

    // Convert to NDC coordinates.
    if (v != NULL) {
        VectorXd temp(4);
        temp << v->x, v->y, v->z, 1;
        temp = PPM * C_inverse * temp;

        new_v = new Vertex;
        new_v->x = temp(0);
        new_v->y = temp(1);
        new_v->z = temp(2);
        new_v->w = temp(3);

        double w = new_v->w;
        new_v->x /= w;
        new_v->y /= w;
        new_v->z /= w;
        new_v->w = 1;
    }

    return new_v;

}

Vertex *MultObjParser::NDC_to_screen(Vertex *v) {
    Vertex *new_v = NULL;
    
    // Apply only to NDC coordinates in screen.
    if (v != NULL // && v->x >= -1 && v->x <= 1 && 
        // v->y >= -1 && v->y <= 1
        ) {
        new_v = new Vertex;
        // Map coordinates from NDC to screen coordinates.
        new_v->x = xres / 2.0 * (v->x + 1);
        new_v->y = yres / -2.0 * (v->y - 1);
        new_v->z = v->z;
    }
    
    return new_v;
}


// Rasterize triangle in color.
void MultObjParser::raster_colored_triangle(Vertex *a, Vertex *b, Vertex *c,
        Color& c_a, Color& c_b, Color& c_c) {
    Vertex *NDC_a = world_to_NDC(a);
    Vertex *NDC_b = world_to_NDC(b);
    Vertex *NDC_c = world_to_NDC(c);

    Vector3d cross = (to_vector(NDC_c) - to_vector(NDC_b)).cross(
                      to_vector(NDC_a) - to_vector(NDC_b));

    // Backface cull
    if (cross(2) < 0) {
        delete NDC_a;
        delete NDC_b;
        delete NDC_c;
        return;
    }

    // Get screen coordinates;
    Vertex *sa = NDC_to_screen(NDC_a);
    Vertex *sb = NDC_to_screen(NDC_b);
    Vertex *sc = NDC_to_screen(NDC_c);

    // Only raster triangles in view.
    if (sa == NULL || sb == NULL || sc == NULL) {
        delete NDC_a;
        delete NDC_b;
        delete NDC_c;
        delete sa;
        delete sb;
        delete sc;
        return;
    }

    // Compute bounding box.
    int x_min = std::max((int) std::min({sa->x, sb->x, sc->x}), 0);
    int x_max = std::min((int) std::max({sa->x, sb->x, sc->x}), xres - 1);
    int y_min = std::max((int) std::min({sa->y, sb->y, sc->y}), 0);
    int y_max = std::min((int) std::max({sa->y, sb->y, sc->y}), yres - 1);

    // Iterate through bounding box.
    for (int x = x_min; x <= x_max; x++) {
        for (int y = y_min; y <= y_max; y++) {
            double alpha = compute_alpha(sa->x, sa->y, sb->x, sb->y, sc->x, sc->y, x, y);
            double beta = compute_beta(sa->x, sa->y, sb->x, sb->y, sc->x, sc->y, x, y);
            double gamma = compute_gamma(sa->x, sa->y, sb->x, sb->y, sc->x, sc->y, x, y);

            // The coefficients are all in [0,1] iff the point is in the triangle.
            if (alpha >= 0 && alpha <= 1 &&
                    beta >= 0 && beta <= 1 &&
                    gamma >= 0 && gamma <= 1) {
                Vector3d NDC = alpha * to_vector(NDC_a) + beta * 
                    to_vector(NDC_b) + gamma * to_vector(NDC_c);

                // Make sure the point is visible.
                if (NDC(0) >= -1 && NDC(0) <= 1 &&
                    NDC(1) >= -1 && NDC(1) <= 1 &&
                    NDC(2) >= -1 && NDC(2) <= 1 &&
                    !(NDC(2) > buffer[x + y * xres])) {
                    
                    // Update line of sight.
                    buffer[x + y * xres] = NDC(2);

                    double R = alpha * c_a(0) + beta * c_b(0) + gamma * c_c(0);
                    double G = alpha * c_a(1) + beta * c_b(1) + gamma * c_c(1);
                    double B = alpha * c_a(2) + beta * c_b(2) + gamma * c_c(2);
                    
                    Color *old = grid[x + y * xres];
                    delete old;
                    grid[x + y * xres] = new Color(R, G, B);
                }

            }
        }
    }

    delete NDC_a;
    delete NDC_b;
    delete NDC_c;
    delete sa;
    delete sb;
    delete sc;
}

// The Gouraud shading algorithm. Requires less computation than Phong, but
// is more pixelated.
void MultObjParser::gouraud_shading(LightProps *lp, Vertex *a, Vertex *b, Vertex *c,
        Vertex *n_a, Vertex *n_b, Vertex *n_c) {
    Vertex *e = new Vertex;
    e->x = T_c.m(0,3);
    e->y = T_c.m(1,3);
    e->z = T_c.m(2,3);

    Color c_a = lighting(a, n_a, lp, lights, e);
    Color c_b = lighting(b, n_b, lp, lights, e);
    Color c_c = lighting(c, n_c, lp, lights, e);

    raster_colored_triangle(a, b, c, c_a, c_b, c_c); 
    delete e;
}

// Apply Gouraud algorithm to all faces.
void MultObjParser::do_gouraud() {
    for (int i = 0; i < objects.size(); i++) {
        std::vector<Vertex*> vertices = objects[i]->op->vertices;
        std::vector<Vertex*> normals = objects[i]->op->normals;
        std::vector<Face*> faces = objects[i]->op->faces;
        LightProps *lp = objects[i]->lp;

        for (int j = 0; j < faces.size(); j++) {
            Vertex *a = vertices[faces[j]->v1];
            Vertex *b = vertices[faces[j]->v2];
            Vertex *c = vertices[faces[j]->v3];
            Vertex *n_a = normals[faces[j]->n1];
            Vertex *n_b = normals[faces[j]->n2];
            Vertex *n_c = normals[faces[j]->n3];

            gouraud_shading(lp, a, b, c, n_a, n_b, n_c);
        }

    }
}

// The more computationally intensive algorithm.
void MultObjParser::phong_shading(LightProps *lp, Vertex *a, Vertex *b, Vertex *c,
        Vertex *n_a, Vertex *n_b, Vertex *n_c) {
    // Create
    Vertex *e = new Vertex;
    e->x = T_c.m(0,3);
    e->y = T_c.m(1,3);
    e->z = T_c.m(2,3);

    Vertex *NDC_a = world_to_NDC(a);
    Vertex *NDC_b = world_to_NDC(b);
    Vertex *NDC_c = world_to_NDC(c);

    Vector3d cross = (to_vector(NDC_c) - to_vector(NDC_b)).cross(
                      to_vector(NDC_a) - to_vector(NDC_b));

    // Backface cull
    if (cross(2) < 0) {
        delete NDC_a;
        delete NDC_b;
        delete NDC_c;
        delete e;
        return;
    }

    // Get screen coordinates;
    Vertex *sa = NDC_to_screen(NDC_a);
    Vertex *sb = NDC_to_screen(NDC_b);
    Vertex *sc = NDC_to_screen(NDC_c);

    // Only raster triangles in view.
    if (sa == NULL || sb == NULL || sc == NULL) {
        delete e;
        delete NDC_a;
        delete NDC_b;
        delete NDC_c;
        delete sa;
        delete sb;
        delete sc;
        return;
    }

    // Compute bounding box.
    int x_min = std::max((int) std::min({sa->x, sb->x, sc->x}), 0);
    int x_max = std::min((int) std::max({sa->x, sb->x, sc->x}), xres - 1);
    int y_min = std::max((int) std::min({sa->y, sb->y, sc->y}), 0);
    int y_max = std::min((int) std::max({sa->y, sb->y, sc->y}), yres - 1);

    for (int x = x_min; x <= x_max; x++) {
        for (int y = y_min; y <= y_max; y++) {
            double alpha = compute_alpha(sa->x, sa->y, sb->x, sb->y, sc->x, sc->y, x, y);
            double beta = compute_beta(sa->x, sa->y, sb->x, sb->y, sc->x, sc->y, x, y);
            double gamma = compute_gamma(sa->x, sa->y, sb->x, sb->y, sc->x, sc->y, x, y);

            if (alpha >= 0 && alpha <= 1 &&
                    beta >= 0 && beta <= 1 &&
                    gamma >= 0 && gamma <= 1) {
                // Compute camera coordinates based on vertices.
                Vector3d NDC = alpha * to_vector(NDC_a) + beta * 
                    to_vector(NDC_b) + gamma * to_vector(NDC_c);

                // Make sure NDC is a visible coordinate.
                if (NDC(0) >= -1 && NDC(0) <= 1 &&
                    NDC(1) >= -1 && NDC(1) <= 1 &&
                    NDC(2) >= -1 && NDC(2) <= 1 &&
                    !(NDC(2) > buffer[x + y * xres])) {
                    
                    // Update depth buffer representing line of sight.
                    buffer[x + y * xres] = NDC(2);

                    Vertex v;
                    v.x = alpha * a->x + beta * b->x + gamma * c->x;
                    v.y = alpha * a->y + beta * b->y + gamma * c->y;
                    v.z = alpha * a->z + beta * b->z + gamma * c->z;
                    Vertex n;
                    n.x = alpha * n_a->x + beta * n_b->x + gamma * n_c->x;
                    n.y = alpha * n_a->y + beta * n_b->y + gamma * n_c->y;
                    n.z = alpha * n_a->z + beta * n_b->z + gamma * n_c->z;
                    double norm = n.x * n.x + n.y * n.y + n.z * n.z;
                    n.x /= norm;
                    n.y /= norm;
                    n.z /= norm;

                    // Compute new color.
                    Color c = lighting(&v, &n, lp, lights, e);
                    
                    Color *old = grid[x + y * xres];
                    delete old;
                    grid[x + y * xres] = new Color(c(0), c(1), c(2));
                }

            }
        }
    }

    delete NDC_a;
    delete NDC_b;
    delete NDC_c;
    delete sa;
    delete sb;
    delete sc;
    delete e;
}

// Apply Phong shading to all surfaces.
void MultObjParser::do_phong() {
    // Iterate over all objects.
    for (int i = 0; i < objects.size(); i++) {
        std::vector<Vertex*> vertices = objects[i]->op->vertices;
        std::vector<Vertex*> normals = objects[i]->op->normals;
        std::vector<Face*> faces = objects[i]->op->faces;
        LightProps *lp = objects[i]->lp;

        // Shade each face.
        for (int j = 0; j < faces.size(); j++) {
            Vertex *a = vertices[faces[j]->v1];
            Vertex *b = vertices[faces[j]->v2];
            Vertex *c = vertices[faces[j]->v3];
            Vertex *n_a = normals[faces[j]->n1];
            Vertex *n_b = normals[faces[j]->n2];
            Vertex *n_c = normals[faces[j]->n3];

            // Use the Phong algorithm.
            phong_shading(lp, a, b, c, n_a, n_b, n_c);
        }

    }
}

// Clean up memory.
MultObjParser::~MultObjParser() {
    for (int i = 0; i < xres * yres; i++) { 
        delete grid[i];
    }
    delete[] grid;
    delete[] buffer;

    for (int i = 0; i < objects.size(); i++) { 
        delete objects[i];
    }

    for (int i = 0; i < lights.size(); i++) {
        delete lights[i];
    }


    for (std::map<std::string, ObjParser*>::iterator it=obj_map.begin(); it!=obj_map.end(); it++) {
        delete it->second;
    }
}

// Print out objects.
void MultObjParser::print() {
    for (int i = 0; i < objects.size(); i++) {
        objects[i]->print();
        if (i != objects.size() - 1) {
            std::cout << "\n";
        }
    }
}

void MultObjParser::debug_camera() {
    std::cerr << "Near: " << cam.near
              << " Far: " << cam.far 
              << " Left: " << cam.left 
              << " Right: " << cam.right 
              << " Top: " << cam.top 
              << " Bottom: " << cam.bottom << std::endl;

    std::cerr << "T_c: " << T_c.v[0] << " "
              << T_c.v[1] << " " << T_c.v[2] << std::endl;

    std::cerr << "R_c: " << R_c.v[0] << " " << R_c.v[1] 
              << " " << R_c.v[2] << " " << R_c.angle << std::endl;

}

// Print out coordinates in ppm format.
void MultObjParser::to_ppm(int xres, int yres) {
    // File header.
    std::cout << "P3\n" << xres << " " << yres << "\n255\n";

    // Grid representing picture.
    int *grid = new int[xres * yres]();
    
    for (int i = 0; i < objects.size(); i++) {
        ObjParser *op = objects[i]->op;

        std::vector<Vertex*> vertices = op->vertices;

        std::vector<Face*> faces = op->faces;
        for (int j = 0; j < faces.size(); j++) {
            Vertex *v1 = vertices[faces[j]->v1];
            Vertex *v2 = vertices[faces[j]->v2];
            Vertex *v3 = vertices[faces[j]->v3];
            Vertex *ndc;

            ndc = world_to_NDC(v1);
            v1 = NDC_to_screen(ndc);
            delete ndc;

            ndc = world_to_NDC(v2);
            v2 = NDC_to_screen(ndc);
            delete ndc;
            
            ndc = world_to_NDC(v3);
            v3 = NDC_to_screen(ndc);
            delete ndc;

            int x0, y0, x1, y1;
            if (v1 != NULL && v2 != NULL) {
                x0 = (int) std::round(v1->x);
                y0 = (int) std::round(v1->y);
                x1 = (int) std::round(v2->x);
                y1 = (int) std::round(v2->y);
                do_bresenham(x0, y0, x1, y1, grid, xres, yres, 1);
            }

            if (v2 != NULL && v3 != NULL) {
                x0 = (int) std::round(v3->x);
                y0 = (int) std::round(v3->y);
                x1 = (int) std::round(v2->x);
                y1 = (int) std::round(v2->y);
                do_bresenham(x0, y0, x1, y1, grid, xres, yres, 1);
            }

            if (v3 != NULL && v1 != NULL) {
                x0 = (int) std::round(v1->x);
                y0 = (int) std::round(v1->y);
                x1 = (int) std::round(v3->x);
                y1 = (int) std::round(v3->y);
                do_bresenham(x0, y0, x1, y1, grid, xres, yres, 1);
            }

            delete v1;
            delete v2;
            delete v3;
        }
    }

    for (int j = 0; j < yres; j++) {
        for (int i = 0; i < xres; i++) {
            switch (grid[i + j * xres]) {
            case 0:
                std::cout << BG_COLOR;
                break;
            case 1:
                std::cout << COLOR_1;
                break;
            default:
                std::cout << BG_COLOR;
                break;
            }
            std::cout << "\n";
        }
    }

    delete[] grid;
}


// Do generalized Bresenham line drawing algorithm.
bool MultObjParser::do_bresenham(int x0, int y0, int x1, int y1, 
        int *grid, int xres, int yres, int color) {
    int e, x, y, dx, dy, x_min, x_max, y_min, y_max, update;
    e = 0;
    dx = x1 - x0;
    dy = y1 - y0;

    
    // Increment or decrement depending on positive or negative slope.
    update = 1;
    if (dx * dy < 0) {
        update = -1;
    }

    // Convert to absolute dx and dy.
    dx = std::abs(dx);
    dy = std::abs(dy);

    // We should iterate by the coordinate that grows faster, so there
    // won't be two pixels in that coordinate with the same value.
    bool iter_x = true;
    if (dy > dx) { 
        iter_x = false;
    }

    // |slope| <= 1, so we can iterate by x without worrying
    // about having a x coordinate with multiple filled pixels.
    if (iter_x) {
        // Determine which point is first in terms of x coordinate.
        if (x0 < x1) {
            x_min = x0;
            x_max = x1;
            y = y0;
        }
        // x1 comes before x0, so we should use x1 for x_min
        // and start with y = y1.
        else {
            x_min = x1;
            x_max = x0;
            y = y1;
        }

        // Iterate from x_min to x_max, incrementing / decrementing
        // y depending on the slope.
        for (x = x_min; x < x_max; x++) {
            grid[x + y * xres] = color;
            if (2 * (e + dy) < dx) {
                e += dy;
            }
            else {
                e += dy - dx;
                y += update;
            }
        }
    }
    // |slope| > 1, so we can iterate by y without worrying
    // about having a y coordinate with multiple filled pixels.
    else {
        // Determine which point is first in terms of y coordinate.
        if (y0 < y1) {
            y_min = y0;
            y_max = y1;
            x = x0;
        }
        // y1 comes before y0, so we should use y1 for y_min
        // and start x = x1.
        else {
            y_min = y1;
            y_max = y0;
            x = x1;
        }

        // Iterate from y_min to y_max, incrementing / decrementing
        // x depending on the slope.
        for (y = y_min; y < y_max; y++) {
            grid[x + y * xres] = color;
            if (2 * (e + dx) < dy) {
                e += dx;
            }
            else {
                e += dx - dy;
                x += update;
            }
        }
    }

    // Succesfully completed.
    return true; 
}
