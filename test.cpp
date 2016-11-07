#include <assert.h>
#include <iostream>
#include <string>
#include <vector>


#include "arcball.h"
#include "obj_parser.h"
#include "mult_obj_parser.h"

#define TOLERANCE 0.00001

int main(int argc, char **argv) {
    // Regression tests for new format.
    if (argc != 2) {
        std::cerr << "Usage: ./test [file.txt]";
        return 1;
    }

    std::cout << "Printing cube.obj\n";
    ObjParser op("cube.obj", "data/");
    ObjParser *new_op = op.new_copy();
    new_op->print();
    delete new_op;

    std::string file_str(argv[1]);
    
    std::cout << "Printing " << file_str << "\n";
    MultObjParser mop(file_str, "data/");
    mop.print();
    //mop.shade(1);

    for (int i = 0; i < mop.lights.size(); i++) {
        mop.lights[i]->debug();
    }

    // Test quatern multiplication.
    // Example taken from testing with mouse.
    Quatern q1(0.884939, -0.363, 0.290843, 0.0228829);
    Quatern q2(0.73952, -0.601555, -0.0706038, 0.293695);
    Quatern q3 = q1 * q2;
    q3.debug();
    
    // Expected value:
    assert ((q3.s - 0.44988) < TOLERANCE);
    double expected_v[] = {-0.71375, 0.24545, 0.477412};
    for (int i = 0; i < 3; i++)
        assert ((q3.v[i] - expected_v[i]) < TOLERANCE); 

    // mop.to_ppm(800, 800);
    return 0;
}
