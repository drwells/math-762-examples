#include <iostream>

// This second example was a quick overview of how objects work in C++. We have
// two definitions of what objects are:
//
// 1. functions and data
// 2. code that models part of the program
//
// In this case, we created a Point class in 3D that models points in Cartesian
// space: we can multiply them by scalars, add them, and do dot products with
// standard mathematical notation.
class Point
{
public:
    // A member function: the first argument is implicitly the current object.
    Point
    operator+(Point other_point)
    {
        // this is a pointer to the object on which we are calling the method.
        std::cout << "location of this = " << this << '\n';
        Point result;
        // We can get values out of this object by dereferencing the pointer
        // and then taking the value of a class member with ->
        result.x = this->x + other_point.x;
        result.y = this->y + other_point.y;
        result.z = this->z + other_point.z;

        // You don't always need to specify this - the compiler can figure out
        // which x it is. The following commented-out code would also work:
        // result.x = x + other_point.x;
        // result.y = y + other_point.y;
        // result.z = z + other_point.z;
        // Since, in the context of this method call, x has to refer to this->x.

        return result;
    }

    // inner product.
    double
    operator*(Point other_point)
    {
        return x * other_point.x + y * other_point.y + z * other_point.z;
    }

    // scalar multiplication. Note that we can have operator* defined in two ways:
    // the compiler will call the correct function based on the type of the input
    // argument.
    Point
    operator*(double scalar)
    {
        Point new_point;

        new_point.x = x * scalar;
        new_point.y = y * scalar;
        new_point.z = z * scalar;

        return new_point;
    }

    // TODO - implement other dimensions
    double x;
    double y;
    double z;
};

int main()
{
    Point point_1;
    point_1.x = 1.0;
    point_1.y = 2.0;
    point_1.z = 3.0;

    Point point_2;
    point_2.x = 0.0;
    point_2.y = 1.0;
    point_2.z = 0.0;
    std::cout << "location of point_1 = " << &point_1 << '\n';
    std::cout << "location of point_2 = " << &point_2 << '\n';

    // Point point_3 = point_2;
    Point point_3 = point_1 + point_2;
    Point point_4 = point_2 + point_1;
    std::cout << "x = " << point_3.x << '\n';
    std::cout << "y = " << point_3.y << '\n';
    std::cout << "z = " << point_3.z << '\n';

    std::cout << point_3 * point_4 << '\n';
    std::cout << (point_3 * (point_4 * 5.0)) << '\n';
}
