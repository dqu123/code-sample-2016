#include <cmath>
#include <sstream>
#include <iostream>

#include "transform.h" // Initialize m as identity matrix.
Transformation::Transformation() : m(4, 4), _normal(4, 4) {
    m << 1.0, 0.0, 0.0, 0.0,
         0.0, 1.0, 0.0, 0.0,
         0.0, 0.0, 1.0, 0.0,
         0.0, 0.0, 0.0, 1.0;

    _normal << 1.0, 0.0, 0.0, 0.0,
         0.0, 1.0, 0.0, 0.0,
         0.0, 0.0, 1.0, 0.0,
         0.0, 0.0, 0.0, 1.0;
}

// Adds a transformation to the current matrix.
bool Transformation::add(const std::string& vector) {
    std::stringstream ss(vector);
    std::string buf;

    if (vector == "") {
        return false;
    }
    
    ss >> buf;
    MatrixXd temp(4, 4);
    temp << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1;
    
    // Translation
    if (buf == "t") {
        ss >> buf;
        double tx = std::stod(buf);
        ss >> buf;
        double ty = std::stod(buf);
        ss >> buf;
        double tz = std::stod(buf);
        temp << 1, 0, 0, tx,
                0, 1, 0, ty,
                0, 0, 1, tz,
                0, 0, 0, 1;

        v[0] = tx;
        v[1] = ty;
        v[2] = tz;

    }
    // Rotation
    else if (buf == "r") {
        ss >> buf;
        double x = std::stod(buf);
        ss >> buf;
        double y = std::stod(buf);
        ss >> buf;
        double z = std::stod(buf);
        ss >> buf;
        double theta = std::stod(buf);

        // Normalize vector.
        double norm = sqrt(x * x + y * y + z * z);
        x /= norm;
        y /= norm;
        z /= norm;

        double c = cos(theta);
        double s = sin(theta);

        temp << x * x + (1 - x * x) * c, x * y * (1 - c) - z * s, x * z * (1 - c) + y * s, 0,
                y * x * (1 - c) + z * s, y * y + (1 - y * y) * c, y * z * (1 - c) - x * s, 0,
                z * x * (1 - c) - y * s, z * y * (1 - c) + x * s, z * z + (1 - z * z) * c, 0,
                                      0,                       0,                       0, 1;
        v[0] = x;
        v[1] = y;
        v[2] = z;
        angle = theta;
     }
    // Scale
    else if (buf == "s") {
        ss >> buf;
        double sx = std::stod(buf);
        ss >> buf;
        double sy = std::stod(buf);
        ss >> buf;
        double sz = std::stod(buf);
        temp << sx,  0,  0, 0,
                 0, sy,  0, 0,
                 0,  0, sz, 0,
                 0,  0,  0, 1;

        v[0] = sx;
        v[1] = sy; 
        v[2] = sz;

    }

    m = temp * m;
    
    // Apply rotations and scaling to normal transformation.
    if (buf != "t") {
        _normal = temp * _normal;
    }

    return true;
}

MatrixXd Transformation::n_transform() {
    return _normal.inverse().transpose();
}

// Print the inverse of matrix m.
void Transformation::print() {
    MatrixXd inv_m = m.inverse();
    std::cout << inv_m << "\n";
}

// Get the matrix for OpenGL.
float *Transformation::get_gl_matrix() {
    // Get transpose matrix.
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            gl_m[j + 4 * i] = m(j, i);
        }
    }

    return gl_m;
}
