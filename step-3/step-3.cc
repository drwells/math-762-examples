/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2020 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------

 *
 * Authors: Wolfgang Bangerth, 1999,
 *          Guido Kanschat, 2011
 */


// @sect3{Many new include files}

// These include files are already known to you. They declare the classes
// which handle triangulations and enumeration of degrees of freedom:
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
// And this is the file in which the functions are declared that create grids:
#include <deal.II/grid/grid_generator.h>

// The next three files contain classes which are needed for loops over all
// cells and to get the information from the cell objects. The first two have
// been used before to get geometric information from cells; the last one is
// new and provides information about the degrees of freedom local to a cell:
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>

// This file contains the description of the Lagrange interpolation finite
// element:
#include <deal.II/fe/fe_q.h>

// And this file is needed for the creation of sparsity patterns of sparse
// matrices, as shown in previous examples:
#include <deal.II/dofs/dof_tools.h>

// The next two files are needed for assembling the matrix using quadrature on
// each cell. The classes declared in them will be explained below:
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>

// The following three include files we need for the treatment of boundary
// values:
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

// We're now almost to the end. The second to last group of include files is
// for the linear algebra which we employ to solve the system of equations
// arising from the finite element discretization of the Laplace equation. We
// will use vectors and full matrices for assembling the system of equations
// locally on each cell, and transfer the results into a sparse matrix. We
// will then use a Conjugate Gradient solver to solve the problem, for which
// we need a preconditioner (in this program, we use the identity
// preconditioner which does nothing, but we need to include the file anyway):
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

// Finally, this is for output to a file and to the console:
#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iostream>

// ...and this is to import the deal.II namespace into the global scope:
using namespace dealii;

// @sect3{The <code>Step3</code> class}

// Instead of the procedural programming of previous examples, we encapsulate
// everything into a class for this program. The class consists of functions
// which each perform certain aspects of a finite element program, a `main`
// function which controls what is done first and what is done next, and a
// list of member variables.

// The public part of the class is rather short: it has a constructor and a
// function `run` that is called from the outside and acts as something like
// the `main` function: it coordinates which operations of this class shall be
// run in which order. Everything else in the class, i.e. all the functions
// that actually do anything, are in the private section of the class:
class Step3
{
public:
  Step3();

  void run();

  // Then there are the member functions that mostly do what their names
  // suggest and whose have been discussed in the introduction already. Since
  // they do not need to be called from outside, they are made private to this
  // class.

private:
  void make_grid();
  void setup_system();
  void assemble_system();
  void solve();
  void output_results() const;

  // And finally we have some member variables. There are variables describing
  // the triangulation and the global numbering of the degrees of freedom (we
  // will specify the exact polynomial degree of the finite element in the
  // constructor of this class)...
  Triangulation<2> triangulation;
  FE_Q<2>          fe;
  DoFHandler<2>    dof_handler;

  // ...variables for the sparsity pattern and values of the system matrix
  // resulting from the discretization of the Laplace equation...
  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

  // ...and variables which will hold the right hand side and solution
  // vectors.
  Vector<double> solution;
  Vector<double> system_rhs;
};

// @sect4{Step3::Step3}

// Here comes the constructor. It does not much more than first to specify
// that we want bi-linear elements (denoted by the parameter to the finite
// element object, which indicates the polynomial degree), and to associate
// the dof_handler variable to the triangulation we use. (Note that the
// triangulation isn't set up with a mesh at all at the present time, but the
// DoFHandler doesn't care: it only wants to know which triangulation it will
// be associated with, and it only starts to care about an actual mesh once
// you try to distribute degree of freedom on the mesh using the
// distribute_dofs() function.) All the other member variables of the Step3
// class have a default constructor which does all we want.
Step3::Step3()
  : fe(1)
  , dof_handler(triangulation)
{}


// @sect4{Step3::make_grid}

// Now, the first thing we've got to do is to generate the triangulation on
// which we would like to do our computation and number each vertex with a
// degree of freedom. We have seen these two steps in step-1 and step-2
// before, respectively.
//
// This function does the first part, creating the mesh.  We create the grid
// and refine all cells five times. Since the initial grid (which is the
// square $[-1,1] \times [-1,1]$) consists of only one cell, the final grid
// has 32 times 32 cells, for a total of 1024.
//
// Unsure that 1024 is the correct number? We can check that by outputting the
// number of cells using the <code>n_active_cells()</code> function on the
// triangulation.
void Step3::make_grid()
{
  GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(5);

  std::cout << "Number of active cells: " << triangulation.n_active_cells()
            << std::endl;
}

// @note We call the Triangulation::n_active_cells() function, rather than
// Triangulation::n_cells(). Here, <i>active</i> means the cells that aren't
// refined any further. We stress the adjective "active" since there are more
// cells, namely the parent cells of the finest cells, their parents, etc, up
// to the one cell which made up the initial grid. Of course, on the next
// coarser level, the number of cells is one quarter that of the cells on the
// finest level, i.e. 256, then 64, 16, 4, and 1. If you called
// <code>triangulation.n_cells()</code> instead in the code above, you would
// consequently get a value of 1365 instead. On the other hand, the number of
// cells (as opposed to the number of active cells) is not typically of much
// interest, so there is no good reason to print it.


// @sect4{Step3::setup_system}

// Next we enumerate all the degrees of freedom and set up matrix and vector
// objects to hold the system data. Enumerating is done by using
// DoFHandler::distribute_dofs(), as we have seen in the step-2 example. Since
// we use the FE_Q class and have set the polynomial degree to 1 in the
// constructor, i.e. bilinear elements, this associates one degree of freedom
// with each vertex. While we're at generating output, let us also take a look
// at how many degrees of freedom are generated:
void Step3::setup_system()
{
  dof_handler.distribute_dofs(fe);
  std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;
  // There should be one DoF for each vertex. Since we have a 32 times 32
  // grid, the number of DoFs should be 33 times 33, or 1089.

  // As we have seen in the previous example, we set up a sparsity pattern by
  // first creating a temporary structure, tagging those entries that might be
  // nonzero, and then copying the data over to the SparsityPattern object
  // that can then be used by the system matrix.
  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  // for general unstructured grids - the sparsity pattern might be complicated.
  DoFTools::make_sparsity_pattern(dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);

  // Note that the SparsityPattern object does not hold the values of the
  // matrix, it only stores the places where entries are. The entries
  // themselves are stored in objects of type SparseMatrix, of which our
  // variable system_matrix is one.
  //
  // The distinction between sparsity pattern and matrix was made to allow
  // several matrices to use the same sparsity pattern. This may not seem
  // relevant here, but when you consider the size which matrices can have,
  // and that it may take some time to build the sparsity pattern, this
  // becomes important in large-scale problems if you have to store several
  // matrices in your program.
  system_matrix.reinit(sparsity_pattern);

  // The last thing to do in this function is to set the sizes of the right
  // hand side vector and the solution vector to the right values:
  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
}

// @sect4{Step3::assemble_system}


// The next step is to compute the entries of the matrix and right hand side
// that form the linear system from which we compute the solution. This is the
// central function of each finite element program and we have discussed the
// primary steps in the introduction already.
//
// The general approach to assemble matrices and vectors is to loop over all
// cells, and on each cell compute the contribution of that cell to the global
// matrix and right hand side by quadrature. The point to realize now is that
// we need the values of the shape functions at the locations of quadrature
// points on the real cell. However, both the finite element shape functions
// as well as the quadrature points are only defined on the reference
// cell. They are therefore of little help to us, and we will in fact hardly
// ever query information about finite element shape functions or quadrature
// points from these objects directly.
//
// Rather, what is required is a way to map this data from the reference cell
// to the real cell. Classes that can do that are derived from the Mapping
// class, though one again often does not have to deal with them directly:
// many functions in the library can take a mapping object as argument, but
// when it is omitted they simply resort to the standard bilinear Q1
// mapping. We will go this route, and not bother with it for the moment (we
// come back to this in step-10, step-11, and step-12).
//
// So what we now have is a collection of three classes to deal with: finite
// element, quadrature, and mapping objects. That's too much, so there is one
// type of class that orchestrates information exchange between these three:
// the FEValues class. If given one instance of each three of these objects
// (or two, and an implicit linear mapping), it will be able to provide you
// with information about values and gradients of shape functions at
// quadrature points on a real cell.
//
// Using all this, we will assemble the linear system for this problem in the
// following function:
void Step3::assemble_system()
{
  QGauss<2> quadrature_formula(fe.degree + 1); // 2D Gauss Rule, computed by tensor product, that uses fe.degree + 1 points in each coordinate direction
  // HW2: wrote QUAD9 as a tensor product of two 1D elements
  // QGauss<2>: tensor product of two 1D quadrature formulas
  FEValues<2> fe_values(fe,
                        quadrature_formula,
                        update_values | update_gradients | update_JxW_values);
  // FEValues : values of a finite element (finite element values)
  //
  // update_values | update_gradients | update_JxW_values
  //
  // bitwise OR
  //
  //
  // update_values     = 00000001 in binary
  // update_gradients  = 00000010 in binary
  // update_JxW_values = 01000000 in binary
  //
  // bitwise or: we get  01000011 as the value of update_values | update_gradients | update_JxW_values
  //
  // tells FEValues which values to compute: we don't need Hessians (third value
  // from right is zero), but we do need values (first value is 1) and gradients
  // (second value is 1).
  //
  // | is bitwise OR, & is bitwise AND

  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);

      cell_matrix = 0;
      cell_rhs    = 0;

      for (unsigned int q_index = 0; q_index < quadrature_formula.size(); ++q_index)
        {
          for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
            for (unsigned int j = 0; j < fe.dofs_per_cell; ++j)
              cell_matrix(i, j) +=
                (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q) (if phi_i is bilinear, its derivatives in x and y are also linear)
                 fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
                 fe_values.JxW(q_index));           // dx (the mapping is linear, so if we take its jacobian we will still get something linear)

          for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
            cell_rhs(i) += (fe_values.shape_value(i, q_index) * // phi_i(x_q) linear
                            1. *                                // f(x_q)     constant
                            fe_values.JxW(q_index));            // dx         linear
        }
      cell->get_dof_indices(local_dof_indices);

      for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
        for (unsigned int j = 0; j < fe.dofs_per_cell; ++j)
          system_matrix.add(local_dof_indices[i],
                            local_dof_indices[j],
                            cell_matrix(i, j));

      for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
        system_rhs(local_dof_indices[i]) += cell_rhs(i);
    }


  // Now almost everything is set up for the solution of the discrete
  // system. However, we have not yet taken care of boundary values (in fact,
  // Laplace's equation without Dirichlet boundary values is not even uniquely
  // solvable, since you can add an arbitrary constant to the discrete
  // solution). We therefore have to do something about the situation.
  //
  // For this, we first obtain a list of the degrees of freedom on the
  // boundary and the value the shape function shall have there. For
  // simplicity, we only interpolate the boundary value function, rather than
  // projecting it onto the boundary. There is a function in the library which
  // does exactly this: VectorTools::interpolate_boundary_values(). Its
  // parameters are (omitting parameters for which default values exist and
  // that we don't care about): the DoFHandler object to get the global
  // numbers of the degrees of freedom on the boundary; the component of the
  // boundary where the boundary values shall be interpolated; the boundary
  // value function itself; and the output object.
  //
  // The component of the boundary is meant as follows: in many cases, you may
  // want to impose certain boundary values only on parts of the boundary. For
  // example, you may have inflow and outflow boundaries in fluid dynamics, or
  // clamped and free parts of bodies in deformation computations of
  // bodies. Then you will want to denote these different parts of the
  // boundary by indicators, and tell the interpolate_boundary_values
  // function to only compute the boundary values on a certain part of the
  // boundary (e.g. the clamped part, or the inflow boundary). By default,
  // all boundaries have a 0 boundary indicator, unless otherwise specified. If
  // sections of the boundary have different boundary conditions, you have to
  // number those parts with different boundary indicators. The function call
  // below will then only determine boundary values for those parts of the
  // boundary for which the boundary indicator is in fact the zero specified as
  // the second argument.
  //
  // The function describing the boundary values is an object of type Function
  // or of a derived class. One of the derived classes is
  // Functions::ZeroFunction, which describes (not unexpectedly) a function
  // which is zero everywhere. We create such an object in-place and pass it to
  // the VectorTools::interpolate_boundary_values() function.
  //
  // Finally, the output object is a list of pairs of global degree of freedom
  // numbers (i.e. the number of the degrees of freedom on the boundary) and
  // their boundary values (which are zero here for all entries). This mapping
  // of DoF numbers to boundary values is done by the <code>std::map</code>
  // class.
  std::map<types::global_dof_index, double> boundary_values;
  VectorTools::interpolate_boundary_values(dof_handler,
                                           0,
                                           Functions::ZeroFunction<2>(),
                                           boundary_values);
  // Now that we got the list of boundary DoFs and their respective boundary
  // values, let's use them to modify the system of equations
  // accordingly. This is done by the following function call:
  MatrixTools::apply_boundary_values(boundary_values,
                                     system_matrix,
                                     solution,
                                     system_rhs);
}


// @sect4{Step3::solve}

// The following function simply solves the discretized equation. As the
// system is quite a large one for direct solvers such as Gauss elimination or
// LU decomposition, we use a Conjugate Gradient algorithm. You should
// remember that the number of variables here (only 1089) is a very small
// number for finite element computations, where 100.000 is a more usual
// number.  For this number of variables, direct methods are no longer usable
// and you are forced to use methods like CG.
void Step3::solve()
{
  // First, we need to have an object that knows how to tell the CG algorithm
  // when to stop. This is done by using a SolverControl object, and as
  // stopping criterion we say: stop after a maximum of 1000 iterations (which
  // is far more than is needed for 1089 variables; see the results section to
  // find out how many were really used), and stop if the norm of the residual
  // is below $10^{-12}$. In practice, the latter criterion will be the one
  // which stops the iteration:
  SolverControl solver_control(1000, 1e-12);
  // Then we need the solver itself. The template parameter to the SolverCG
  // class is the type of the vectors, and leaving the empty angle brackets
  // would indicate that we are taking the default argument (which is
  // <code>Vector@<double@></code>). However, we explicitly mention the template
  // argument:
  SolverCG<Vector<double>> solver(solver_control);

  // Now solve the system of equations. The CG solver takes a preconditioner
  // as its fourth argument. We don't feel ready to delve into this yet, so we
  // tell it to use the identity operation as preconditioner:
  solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
  // Now that the solver has done its job, the solution variable contains the
  // nodal values of the solution function.
}


// @sect4{Step3::output_results}

// The last part of a typical finite element program is to output the results
// and maybe do some postprocessing (for example compute the maximal stress
// values at the boundary, or the average flux across the outflow, etc). We
// have no such postprocessing here, but we would like to write the solution
// to a file.
void Step3::output_results() const
{
  // To write the output to a file, we need an object which knows about output
  // formats and the like. This is the DataOut class, and we need an object of
  // that type:
  DataOut<2> data_out;
  // Now we have to tell it where to take the values from which it shall
  // write. We tell it which DoFHandler object to use, and the solution vector
  // (and the name by which the solution variable shall appear in the output
  // file). If we had more than one vector which we would like to look at in
  // the output (for example right hand sides, errors per cell, etc) we would
  // add them as well:
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");
  // After the DataOut object knows which data it is to work on, we have to
  // tell it to process them into something the back ends can handle. The
  // reason is that we have separated the frontend (which knows about how to
  // treat DoFHandler objects and data vectors) from the back end (which knows
  // many different output formats) and use an intermediate data format to
  // transfer data from the front- to the backend. The data is transformed
  // into this intermediate format by the following function:
  data_out.build_patches();

  // Now we have everything in place for the actual output. Just open a file
  // and write the data into it, using VTK format (there are many other
  // functions in the DataOut class we are using here that can write the
  // data in postscript, AVS, GMV, Gnuplot, or some other file
  // formats):
  std::ofstream output("solution.vtk");
  data_out.write_vtk(output);
}


// @sect4{Step3::run}

// Finally, the last function of this class is the main function which calls
// all the other functions of the <code>Step3</code> class. The order in which
// this is done resembles the order in which most finite element programs
// work. Since the names are mostly self-explanatory, there is not much to
// comment about:
void Step3::run()
{
  make_grid();
  setup_system(); // 1. distribute (number) DoFs 2. setup sparsity pattern 3. use that to initialize our matrix and vectors
  assemble_system(); //
  solve();
  output_results();
}


// @sect3{The <code>main</code> function}

// This is the main function of the program. Since the concept of a
// main function is mostly a remnant from the pre-object oriented era
// before C++ programming, it often does not do much more than
// creating an object of the top-level class and calling its principle
// function.
//
// Finally, the first line of the function is used to enable output of
// some diagnostics that deal.II can generate.  The @p deallog
// variable (which stands for deal-log, not de-allog) represents a
// stream to which some parts of the library write output. For
// example, iterative solvers will generate diagnostics (starting
// residual, number of solver steps, final residual) as can be seen
// when running this tutorial program.
//
// The output of @p deallog can be written to the console, to a file,
// or both. Both are disabled by default since over the years we have
// learned that a program should only generate output when a user
// explicitly asks for it. But this can be changed, and to explain how
// this can be done, we need to explain how @p deallog works: When
// individual parts of the library want to log output, they open a
// "context" or "section" into which this output will be placed. At
// the end of the part that wants to write output, one exits this
// section again. Since a function may call another one from within
// the scope where this output section is open, output may in fact be
// nested hierarchically into these sections. The LogStream class of
// which @p deallog is a variable calls each of these sections a
// "prefix" because all output is printed with this prefix at the left
// end of the line, with prefixes separated by colons. There is always
// a default prefix called "DEAL" (a hint at deal.II's history as the
// successor of a previous library called "DEAL" and from which the
// LogStream class is one of the few pieces of code that were taken
// into deal.II).
//
// By default, @p logstream only outputs lines with zero prefixes --
// i.e., all output is disabled because the default "DEAL" prefix is
// always there. But one can set a different maximal number of
// prefixes for lines that should be output to something larger, and
// indeed here we set it to two by calling
// LogStream::depth_console(). This means that for all screen output,
// a context that has pushed one additional prefix beyond the default
// "DEAL" is allowed to print its output to the screen ("console"),
// whereas all further nested sections that would have three or more
// prefixes active would write to @p deallog, but @p deallog does not
// forward this output to the screen. Thus, running this example (or
// looking at the "Results" section), you will see the solver
// statistics prefixed with "DEAL:CG", which is two prefixes. This is
// sufficient for the context of the current program, but you will see
// examples later on (e.g., in step-22) where solvers are nested more
// deeply and where you may get useful information by setting the
// depth even higher.
int main()
{
  deallog.depth_console(2);

  Step3 laplace_problem;
  laplace_problem.run();

  return 0;
}
