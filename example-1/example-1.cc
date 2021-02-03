// example program 1
//
// This example does the following:
// 1. Creates a quadrilateral discretization of the disk
// 2. Refines (by splitting each quadrilateral into 4) each cell five times
// 3. Saves the output to the disk in SVG format (which can be viewed with any web browser)
//
// New ideas:
// 1. Namespaces: In C++, libraries put all their functions and classes into namespaces.
//    For example: everything in the standard library (like std::ofstream, used below)
//    is in the "std" (short for "standard") namespace. Similarly, everything in
//    deal.II is available in the "dealii" namespace (like dealii::GridOut).
// 2. Headers: in C++, in order to use a class or function, its declaration must be available.
//    For example, here we include the "deal.II/grid/grid_generator.h" header so
//    that the functions in dealii::GridGenerator are available for us to use.

// Need the declaration of dealii::Triangulation
#include <deal.II/grid/tria.h>

// Need declarations of various GridGenerator functions
#include <deal.II/grid/grid_generator.h>

// deal.II header for writing grids
#include <deal.II/grid/grid_out.h>

// allows me to read and write files: defines std::ofstream
#include <fstream>

int main()
{
    // Create an empty triangulation
    dealii::Triangulation<2> tria;
    // Replace it with a grid of a circle
    dealii::GridGenerator::hyper_ball(tria);
    // split each cell into 4 cells, five times. As a result, each cell in the
    // initial grid will ultimately have 4^5 = 1024 active child cells by the
    // time we get to the end.

    for (unsigned int i = 0; i < 5; ++i)
    {
        // Create a GridOut object, and then write an SVG file.
        dealii::GridOut go;
        std::ofstream out("ball-" + std::to_string(i) + ".svg");
        go.write_svg(tria, out);
        tria.refine_global(1);
    }

    // save the last grid
    std::ofstream out("ball-5.svg");
    dealii::GridOut go;
    go.write_svg(tria, out);
    tria.refine_global(1);
}
