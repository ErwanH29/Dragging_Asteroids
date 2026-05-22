// Class for storing data

#include <iostream>

#ifndef VEC_3D
#define VEC_3D

class vec{
    private:
        double array[3];  // The actual variable.

    public:
        vec(double x, double y, double z){  // Constructor
            array[0] = x;
            array[1] = y;
            array[2] = z;
        }

        inline void print(){  // Use inline as overhead ~ function call
            std::cout << array[0] << " ";
            std::cout << array[1] << " ";
            std::cout << array[2] << std::endl;
        }
        
        // Overloading operators for vector operations
        // return_type operator<symbol>(arguments) const_optional {
        //     implementation
        //  }
        inline double& operator[](int i){ // Modification operator
            return array[i];
        }

        inline const double& operator[](int i) const { // Read-only operator
            return array[i];
        }

        double operator*(const vec& other) {
            return array[0]*other[0] 
              + array[1]*other[1] 
              + array[2]*other[2];
        }

        inline vec operator-(const vec& other) const {
            return vec(
                array[0] - other[0],
                array[1] - other[1],
                array[2] - other[2]
            );
        }

        inline vec operator*(double scalar) const {
            return vec(
                scalar * array[0],
                scalar * array[1],
                scalar * array[2]
            );
        }

};

#endif