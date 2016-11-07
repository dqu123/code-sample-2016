#ifndef ARCBALL_H
#define ARCBALL_H

#include <Eigen/Dense>

#include <vector>

using namespace Eigen;

class Quatern;

// Function prototypes
Quatern compute_rotation_matrix(int x_cur, int y_cur, int x, int y, int xres, int yres); 


class Quatern {
    public:
        Quatern(double s, double x, double y, double z);
        Quatern() : Quatern(1, 0, 0, 0) {};
        Quatern(double theta, Vector3d& u);
        ~Quatern() {};

        void normalize();

        // Operators
        Quatern & operator*=(const Quatern& rhs);
        Quatern operator*(const Quatern& rhs);

        // OpenGL
        double *get_gl_rotation();

        // Debug
        void debug();

        double s;
        Vector3d v;
    private:
        // Sync rotation matrix with Quaternion values.
        void update_rotation_matrix();

        std::vector<double> rotation_matrix;
};

#endif // ARCBALL_H
