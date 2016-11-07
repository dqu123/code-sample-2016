#ifndef TRANSFORM_H
#define TRANSFORM_H

#include <string>
#include <Eigen/Dense>

using namespace Eigen;

class Transformation {
    public:
        Transformation();
        ~Transformation() {};
        bool add(const std::string& vector);
        void print();
        // Matrix used to transform normal vectors.
        MatrixXd n_transform();
        // OpenGL
        float *get_gl_matrix();

        double v[3];
        double angle;

        MatrixXd m;
    private:
        // Intermediate matrix to compute normal transform matrix.
        MatrixXd _normal;

        // OpenGL
        float gl_m[16];
};

#endif // TRANSFORM_H
