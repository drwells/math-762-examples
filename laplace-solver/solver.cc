#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_update_flags.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

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
    // tell the DoFHandler about the finite element and triangulation
    dof_handler.initialize(triangulation, finite_element);

    // set up the constraints:
    dealii::VectorTools::interpolate_boundary_values(
        dof_handler,
        0, // default boundary id
        dealii::Functions::ZeroFunction<dim>(),
        constraints);
    constraints.close();

    dealii::DynamicSparsityPattern dynamic_sparsity_pattern(dof_handler.n_dofs());
    dealii::DoFTools::make_sparsity_pattern(dof_handler,
                                            dynamic_sparsity_pattern,
                                            constraints,
                                            /* keep_constrained_dofs */ false);

    /*
       replace whatever I have with (keep_constrained_dofs = true)

       1 0 (not stored) (not stored) (not stored)
       0
       (not stored)
       (not stored)

       etc.

       instead do (keep_constrained_dofs = false)

       1 (not stored) (not stored) (not stored) (not stored)
       (not stored)
       (not stored)
       (not stored)
    */

    sparsity_pattern.copy_from(dynamic_sparsity_pattern);

    system_matrix.reinit(sparsity_pattern);

    // up to this point: solution and system_rhs are empty vectors. Now that we
    // know how many degrees of freedom there are we can set them to the correct
    // length:
    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
}


template <int dim>
void
Laplace<dim>::build_system()
{
    // setup quadrature
    // values (and gradients) of our basis functions on the physical cells (and Jacobian)
    // a mapping from the reference to physical cell
    //
    //
    // actually put things in the sparse matrix (and system_rhs)

    const dealii::QGauss<dim> quadrature(finite_element.degree + 1);
    dealii::MappingQGeneric<dim> mapping(1); // weird name, but this is a bilinear/trilinear mapping

    dealii::FEValues<dim> fe_values(mapping,
                                    finite_element,
                                    quadrature,
                                    dealii::update_values | dealii::update_gradients |
                                    dealii::update_JxW_values);

    std::vector<dealii::types::global_dof_index> cell_dof_indices(finite_element.dofs_per_cell);
    // int : 32 bits (about -2 billion to to 2 billion)
    // unsigned int : 32 bits (up to about 4 billion)
    //
    // long int (or just long): 64 bits
    // unsigned long int (or just long): 64 bits
    //
    // int could be 16 on very small machines


    // loop over: elements, quadrature points, test functions, trial functions
    for (const auto &cell : dof_handler.active_cell_iterators())
    {
        fe_values.reinit(cell);
        // cell is:
        // 1. a level number (for us: number of global refinements)
        // 2. an index inside that level
        // 3. a pointer to the DoFHandler
        cell_matrix = 0.0;
        cell_rhs = 0.0;

        for (unsigned int quad_point_n = 0; quad_point_n < quadrature.size();
             ++quad_point_n)
        {
            for (unsigned int i = 0; i < finite_element.dofs_per_cell; ++i)
            {
                for (unsigned int j = 0; j < finite_element.dofs_per_cell; ++j)
                {
                    cell_matrix(i, j) +=
                        fe_values.shape_grad(i, quad_point_n) *
                        fe_values.shape_grad(j, quad_point_n) *
                        fe_values.JxW(quad_point_n);
                }

                cell_rhs(i) += 1.0 // TODO: fix this with the actual forcing function
                    * fe_values.shape_value(i, quad_point_n)
                    * fe_values.JxW(quad_point_n);
            }
        }

        cell->get_dof_indices(cell_dof_indices);
        constraints.distribute_local_to_global(cell_matrix, cell_rhs, cell_dof_indices,
                                               system_matrix, system_rhs);
    }
}


template <int dim>
void
Laplace<dim>::solve_system()
{
    dealii::SolverControl solver_control(1000, 1e-12*system_rhs.l2_norm());
    dealii::SolverCG<dealii::Vector<double>>    solver(solver_control); // Same as SolverCG<>

    dealii::PreconditionJacobi<> preconditioner; // Same as PreconditionJacobi<SparseMatrix<double>>
    preconditioner.initialize(system_matrix);

    // We use CG (a Krylov method) here because we want our code to work
    // reasonably well in 3D. The alternative is to use SparseDirectUMFPACK.
    solver.solve(system_matrix, solution, system_rhs, preconditioner);

    constraints.distribute(solution);
}


template <int dim>
void
Laplace<dim>::make_pretty_graphics()
{
    // "Always have a working code"
    std::ofstream out
        ("grid-" + std::to_string(n_refinements) + ".svg");
    dealii::GridOut go;
    go.write_svg(triangulation, out);

    // write the solution
    {
      dealii::DataOut<dim> data_out;
      data_out.attach_dof_handler(dof_handler);
      data_out.add_data_vector(solution, "solution");
      data_out.add_data_vector(system_rhs, "system_rhs");
      // All finite element spaces are subsets of a discontinuous polynomial space
      data_out.build_patches();
      std::ofstream output("solution.vtu");
      data_out.write_vtu(output);
    }
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
