// Class for storing data

#include <iostream>

#ifndef VEC_3D
#define VEC_3D

class vec{
    private:
        double array[3];
    public:
        vec(double x, double y, double z){
            array[0] = x;
            array[1] = y;
            array[2] = z;
        }

        inline void print(){  // Use inline as overhead ~ function call
            std::cout << array[0] << " ";
            std::cout << array[1] << " ";
            std::cout << array[2] << std::endl;
        }

        double &operator[](int i){
            return array[i];
        }
};

#endif