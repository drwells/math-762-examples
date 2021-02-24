#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <iostream>
#include <fstream>

// goal today: write a Laplace solver from scratch

template <int dim>
class Laplace
{
  public:
        Laplace(const unsigned int n_refinements);
        void run();

    protected:

        void setup_triangulation();
        void setup_dofs(); // remember to BCs here

        void build_system();
        void solve_system();

        void make_pretty_graphics();

        // needed for setup_triangulation
        unsigned int n_refinements;

        // needed for setup_triangulation
        dealii::Triangulation<dim> triangulation;

        // needed for setup_dofs
        dealii::AffineConstraints<double> constraints;

        // needed for setup_dofs
        dealii::FE_Q<dim>       finite_element; // FE_Q - tensor product lagrange interpolatory polynomials
        // needed for setup_dofs
        dealii::DoFHandler<dim> dof_handler;

        // needed for build_system
        dealii::SparsityPattern      sparsity_pattern;
        // needed for build_system
        dealii::SparseMatrix<double> system_matrix;

        // needed for build_system
        dealii::Vector<double>       system_rhs;

        // needed for solve_system
        dealii::Vector<double>       solution;

        // The order in which we declare and set things up is important: they
        // are created in the order they are listed in and destroyed in the
        // reverse order.
        //
        // This is called RAII - Resource Allocation Is Initialization
};



template <int dim>
Laplace<dim>::Laplace(const unsigned int n)
    : n_refinements(n)
    , finite_element(1)
{}

template <int dim>
void
Laplace<dim>::setup_triangulation()
{
    std::vector<unsigned int> holes;
    holes.push_back(3);
    holes.push_back(2);
    dealii::GridGenerator::cheese(triangulation, holes);

    triangulation.refine_global(n_refinements);
}


template <int dim>
void
Laplace<dim>::setup_dofs()
{
    // question: how do I tell the DoFHandler about the Triangulation and the
    // finite element?
    dof_handler.initialize(triangulation, finite_element);

    // TODO finish this next class
}


template <int dim>
void
Laplace<dim>::build_system()
{}


template <int dim>
void
Laplace<dim>::solve_system()
{}


template <int dim>
void
Laplace<dim>::make_pretty_graphics()
{
    // "Always have a working code"
    std::ofstream out
        ("grid-" + std::to_string(n_refinements) + ".svg");
    dealii::GridOut go;
    go.write_svg(triangulation, out);
}


template <int dim>
void
Laplace<dim>::run()
{
    setup_triangulation();
    setup_dofs();
    build_system();
    solve_system();
    make_pretty_graphics();
}

int main()
{
    for (unsigned int i = 0; i < 5; ++i)
    {
        Laplace<2> laplace_solver(i);
        laplace_solver.run();
    }
}
