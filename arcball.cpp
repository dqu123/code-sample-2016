#include <cmath>
#include <iostream>

#include "arcball.h"

// Compute Arcball rotation matrix.
Quatern compute_rotation_matrix(int x_cur, int y_cur, int x, int y, int xres, int yres) { 
    Matrix4d m;

    // Transform x from screen coordinates to NDC
    // [0, xres] -> [-1, 1]
    // (x_cur - (xres / 2.0)) -> [-xres/2.0, xres/2.0]
    // * 2.0 / xres -> [-1, 1] NDC coordinates
    double x_cur_ndc = 2.0 * x_cur / xres - 1;
    double x_ndc = 2.0 * x / xres - 1;
    
    // Transform y from screen coordinates to NDC
    // [0, yres] -> [1, -1]
    // (y_cur - (yres / 2.0)) -> [-yres/2.0, yres/2.0]
    // * -2.0 / yres -> [1, -1]
    double y_cur_ndc = -2.0 * y_cur / yres + 1;
    double y_ndc = -2.0 * y / yres + 1;

    // Compute z coordinate based on sphere equation.
    double cur_square_sum = x_cur_ndc * x_cur_ndc + y_cur_ndc * y_cur_ndc;
    double z_cur_ndc = (cur_square_sum <= 1) ? std::sqrt(1 - cur_square_sum) : 0;
    
    double square_sum = x_ndc * x_ndc + y_ndc * y_ndc;
    double z_ndc = (square_sum <= 1) ? std::sqrt(1 - square_sum) : 0;

    Vector3d v_cur(x_cur_ndc, y_cur_ndc, z_cur_ndc);
    Vector3d v(x_ndc, y_ndc, z_ndc);

    // Generate rotation angle and axis.
    double theta = acos(std::min(1.0, v.dot(v_cur) / (v.norm() * v_cur.norm())));
    Vector3d u = v.cross(v_cur);
    
    // Normalize vector.
    double norm = sqrt(u(0) * u(0) + u(1) * u(1) + u(2) * u(2));
    for (int i = 0; i < 3; i++) { 
        u(i) /= norm;
    }

    Quatern q(theta, u);
    q.normalize();
    return q;

}

// Construct Quatern from components.
Quatern::Quatern(double s, double x, double y, double z) : v(x, y, z) { 
    this->s = s;
}

// Construct Quatern from rotation.
Quatern::Quatern(double theta, Vector3d& u) {
    s = cos(theta / 2);
    for (int i = 0; i < 3; i++)
        v(i) = u(i) * sin(theta / 2);
}

// Normalize Quatern.
void Quatern::normalize() {
    double norm = sqrt(s * s + v.norm() * v.norm());
    s /= norm;
    v /= norm;
}

// Multiply and update current Quatern. 
Quatern& Quatern::operator*=(const Quatern& rhs) { 
    double old_s = s;
    s = s * rhs.s - v.dot(rhs.v);
    v = old_s * rhs.v + rhs.s * v + v.cross(rhs.v);

    return *this;
}

// Multiply two Quaterns.
Quatern Quatern::operator*(const Quatern& other) { 
    Quatern q = *this;
    q *= other;
    return q;
}

// Get array for OpenGL.
double *Quatern::get_gl_rotation() {
    update_rotation_matrix();
    return &rotation_matrix[0];
}

// Compute rotation in matrix form.
void Quatern::update_rotation_matrix() { 
    double x = v(0);
    double y = v(1);
    double z = v(2);

    // Transpose matrix (since OpenGL uses column major format).
    // Converts from quaternion to rotation matrix.
    rotation_matrix = {1 - 2 * y * y - 2 * z * z, 2 * (x * y + z * s), 2 * (x * z - y * s), 0,
                       2 * (x * y - z * s), 1 - 2 * x * x - 2 * z * z, 2 * (y * z + x * s), 0,
                       2 * (x * z + y * s), 2 * (y * z - x * s), 1 - 2 * x * x - 2 * y * y, 0,
                                         0,                   0,                         0, 1};

}

// Pipe debug messages to cerr.
void Quatern::debug() {
    std::cerr << "s: " << s << "\nv:\n" << v;

    double norm = sqrt(s * s + v.norm() * v.norm());
    std::cerr << "\nNorm: " << norm << "\n";
    update_rotation_matrix();
    for (int i = 0; i < 4; i++) { 
        for (int j = 0; j < 4; j++) { 
            std::cerr << rotation_matrix[j + 4 * i] << " ";
        }
        std::cerr << "\n";
    }

}
