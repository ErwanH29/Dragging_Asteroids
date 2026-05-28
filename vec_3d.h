// Class for storing data

#ifndef VEC_3D
#define VEC_3D

#include <iostream>

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
        inline double& operator[](int i){
            return array[i];
        }
        
        // First const is return type, second const is so function doesn't modify object
        inline const double& operator[](int i) const {
            return array[i];
        }

        double operator*(const vec& other) const {
            return array[0]*other[0] 
              + array[1]*other[1] 
              + array[2]*other[2];
        }

        vec operator*(double scalar) const {
            return vec(
                scalar * array[0],
                scalar * array[1],
                scalar * array[2]
            );
        }

        vec operator-(const vec& other) const {
            return vec(
                array[0] - other[0],
                array[1] - other[1],
                array[2] - other[2]
            );
        }

        vec operator+(double scalar) const {
            return vec(
                scalar + array[0],
                scalar + array[1],
                scalar + array[2]
            );
        }

        vec& operator+=(const vec& other){
            array[0] += other[0];
            array[1] += other[1];
            array[2] += other[2];
            return *this;  // Return pointer since modifying reference
        }

        vec& operator-=(const vec& other){
            array[0] -= other[0];
            array[1] -= other[1];
            array[2] -= other[2];
            return *this;
        }

        inline vec& operator*=(double scalar){
            array[0] *= scalar;
            array[1] *= scalar;
            array[2] *= scalar;
            return *this;
        }

};

#endif