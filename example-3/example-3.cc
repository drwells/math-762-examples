#include <iostream>

#include <algorithm> // std::sort
#include <array>  // std::array
#include <vector> // std::vector

// This example showed how templates containers, and iterators work in C++. 
// Consider the point class we used in example-2. It might be useful to use
// floats, rather than doubles, to decrease memory usage. Similarly, we might
// want to express coordinates with integers instead of floating-point values.
// We can achieve both by defining Point as a class template:
template <typename T>
class Point
{
    T x;
    T y;
    T z;
};

int main()
{
    // The class template named vector in namespace std:
    std::vector<int> ints;
    std::vector<double> doubles;

    // for a vector:
    // 1. We want the number of entries in the vector
    // 2. We want to append to the vector
    // 3. We want to pop entries off the back of the vector
    // etc etc etc
    //
    // C++ has template classes, where the type is an argument.
    //
    // The motivation for this is that we usually want to do the same types of
    // things with arrays independent of their underlying type: add things,
    // remove things, compute the size, check if its empty, etc. Implementing
    // std::vector as a template lets us use one piece of code for every conceivable
    // type that we could put in an array.

    // Templates:
    Point<int> point_1; // contains 3 integers
    Point<double> point_2; // contains 3 doubles
    Point<float> point_3; // contains 3 floats
    Point<Point<Point<int>>> point_4; // contains 3 Point<Point<int>>

    {
        std::vector<float> my_vector;
        std::cout << "size = " << my_vector.size() << std::endl;
        my_vector.push_back(3);
        my_vector.push_back(2);
        my_vector.push_back(1);
        std::cout << "size = " << my_vector.size() << std::endl;

        // in C++, standard containers offer similar interfaces.
        std::cout << "first entry in my_vector = " << my_vector.front() << '\n';
        std::cout << "last entry in my_vector = " << my_vector.back() << '\n';

        // All containers offer a similar interface. std::array is a type of 
        // array that has a fixed number of entries (i.e., push_back doesn't 
        // make sense here). The number of entries is encoded in the type with 
        // its second template argument:
        std::array<float, 3> my_array;
        my_array[0] = 3;
        my_array[1] = 2;
        my_array[2] = 1;
        // contains exactly 3 floats

        // same for both containers:
        std::sort(my_vector.begin(), my_vector.end());
        std::cout << "first entry in my_vector = " << my_vector.front() << '\n';
        std::cout << "last entry in my_vector = " << my_vector.back() << '\n';

        std::sort(my_array.begin(), my_array.end());
        std::cout << "first entry in my_array = " << my_array.front() << '\n';
        std::cout << "last entry in my_array = " << my_array.back() << '\n';

        // We access elements of a container with iterators
        std::cout << "elements of my_vector:\n";
        for (std::vector<float>::iterator it = my_vector.begin();
             it < my_vector.end();
             ++it)
        {
            // Much like this is the memory address of the current object, an
            // iterator is a generalization of a memory address of a container
            // entry.
            //
            // unary * is the opposite of unary &
            //
            // &a = address of a
            // *a = value at the address given by a.
            std::cout << *it << '\n';
        }
        // C++ uses half-open ranges: my_vector.begin() points to the first entry
        // my_vector.end() points to ONE PAST the last entry. This is useful 
        // because it allows us to easily represent empty containers (i.e., 
        // begin() == end()).
    }
}
