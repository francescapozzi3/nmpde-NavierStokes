#ifndef NAVIER_STOKES_2D_HPP
#define NAVIER_STOKES_2D_HPP

#include <limits>
#include <utility>
#include <cmath>

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

using namespace dealii;

// Class implementing a solver for the Stokes problem.
class NavierStokes2D
{
public:
  // Physical dimension (1D, 2D, 3D)
  static constexpr unsigned int dim = 2;

    struct ProblemData
  {
    double nu  = 0.001;
    double rho = 1.0;

    double Um = 0.0;

    unsigned int choice = 1;

    // Forcing term
    std::function<Tensor<1, dim>(const Point<dim> &, const double &)> f;

    // Neumann BC on outlet
    std::function<Tensor<1, dim>(const Point<dim> &, const double &)> h;

    // Inlet profile
    std::function<Tensor<1, dim>(const Point<dim> &, const double &)> u_in;

    // Initial condition
    std::function<Tensor<1, dim>(const Point<dim> &)> u0;
  };


  // Function for inlet velocity.
   class InletVelocity : public Function<dim>
  {
  public:
    InletVelocity(const std::function<Tensor<1, dim>(const Point<dim> &, const double &)> &fun)
      : Function<dim>(dim + 1)
      , fun(fun)
    {}

    void set_time(const double t) override
    { 
      time = t; 
    }

    void vector_value(const Point<dim> &p, Vector<double> &values) const override
    {
      const Tensor<1, dim> u = fun(p, time);

      values = 0.0;
      for (unsigned int d = 0; d < dim; ++d)
        values[d] = u[d];
    }

  private:
    std::function<Tensor<1, dim>(const Point<dim> &, const double &)> fun;
    double time = 0.0;
  };

  // Identity preconditioner.
  class PreconditionIdentity
  {
  public:
    // Application of the preconditioner
    void
    vmult(TrilinosWrappers::MPI::BlockVector       &dst,
          const TrilinosWrappers::MPI::BlockVector &src) const
    {
      dst = src;
    }

  protected:
  };

  // Block-diagonal preconditioner.
  class PreconditionBlockDiagonal
  {
  public:
    // Initialize the preconditioner.
    void
    initialize(const TrilinosWrappers::SparseMatrix &velocity_stiffness_,
               const TrilinosWrappers::SparseMatrix &pressure_mass_)
    {
      velocity_stiffness = &velocity_stiffness_;
      pressure_mass      = &pressure_mass_;

      preconditioner_velocity.initialize(velocity_stiffness_);
      preconditioner_pressure.initialize(pressure_mass_);
    }

    // Application of the preconditioner.
    void
    vmult(TrilinosWrappers::MPI::BlockVector       &dst,
          const TrilinosWrappers::MPI::BlockVector &src) const
    {
      SolverControl solver_control_velocity(2000, 1e-2 * src.block(0).l2_norm());
      SolverGMRES<TrilinosWrappers::MPI::Vector> solver_gmres_velocity(
      solver_control_velocity);
      solver_gmres_velocity.solve(*velocity_stiffness,
                                dst.block(0),
                                src.block(0),
                                preconditioner_velocity);

      SolverControl solver_control_pressure(2000, 1e-2 * src.block(1).l2_norm());
      SolverCG<TrilinosWrappers::MPI::Vector> solver_cg_pressure(
        solver_control_pressure);
      solver_cg_pressure.solve(*pressure_mass,
                               dst.block(1),
                               src.block(1),
                               preconditioner_pressure);

    }

  protected:
    // Velocity stiffness matrix.
    const TrilinosWrappers::SparseMatrix *velocity_stiffness;

    // Preconditioner used for the velocity block.
    TrilinosWrappers::PreconditionILU preconditioner_velocity;

    // Pressure mass matrix.
    const TrilinosWrappers::SparseMatrix *pressure_mass;

    // Preconditioner used for the pressure block.
    TrilinosWrappers::PreconditionILU preconditioner_pressure;
  };

  // Block-triangular preconditioner.
  class PreconditionBlockTriangular
  {
  public:
    // Initialize the preconditioner.
    void
    initialize(const TrilinosWrappers::SparseMatrix &velocity_stiffness_,
               const TrilinosWrappers::SparseMatrix &pressure_mass_,
               const TrilinosWrappers::SparseMatrix &B_)
    {
      velocity_stiffness = &velocity_stiffness_;
      pressure_mass      = &pressure_mass_;
      B                  = &B_;

      preconditioner_velocity.initialize(velocity_stiffness_);
      preconditioner_pressure.initialize(pressure_mass_);
    }

    // Application of the preconditioner.
    void
    vmult(TrilinosWrappers::MPI::BlockVector       &dst,
          const TrilinosWrappers::MPI::BlockVector &src) const
    {
      SolverControl solver_control_velocity(2000, 1e-2 * src.block(0).l2_norm());

    SolverGMRES<TrilinosWrappers::MPI::Vector> solver_gmres_velocity(
    solver_control_velocity);

    solver_gmres_velocity.solve(*velocity_stiffness,
                                dst.block(0),
                                src.block(0),
                                preconditioner_velocity);

      tmp.reinit(src.block(1));
      B->vmult(tmp, dst.block(0));
      tmp.sadd(-1.0, src.block(1));

      SolverControl solver_control_pressure(2000, 1e-2 * src.block(1).l2_norm());
      SolverCG<TrilinosWrappers::MPI::Vector> solver_cg_pressure(solver_control_pressure);
      solver_cg_pressure.solve(*pressure_mass,
                               dst.block(1),
                               tmp,
                               preconditioner_pressure);
    }


  protected:

    // Velocity stiffness matrix.
    const TrilinosWrappers::SparseMatrix *velocity_stiffness;

    // Preconditioner used for the velocity block.
    TrilinosWrappers::PreconditionILU preconditioner_velocity;

    // Pressure mass matrix.
    const TrilinosWrappers::SparseMatrix *pressure_mass;

    // Preconditioner used for the pressure block.
    TrilinosWrappers::PreconditionILU preconditioner_pressure;

    // B matrix.
    const TrilinosWrappers::SparseMatrix *B;

    // Temporary vector.
    mutable TrilinosWrappers::MPI::Vector tmp;
  };


 // Constructor.
  NavierStokes2D(const std::string                               &mesh_file_name_,
                 const unsigned int                              &degree_velocity_,
                 const unsigned int                              &degree_pressure_,
                 const double                                    &T_,
                 const double                                    &theta_,
                 const double                                    &delta_t_,
                 const ProblemData                               &problem_data_)
    : mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
    , mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
    , pcout(std::cout, mpi_rank == 0)
    , mesh_file_name(mesh_file_name_)
    , degree_velocity(degree_velocity_)
    , degree_pressure(degree_pressure_)
    , T(T_)
    , theta(theta_)
    , delta_t(delta_t_)
    , data(problem_data_)
    , mesh(MPI_COMM_WORLD)
  {}

  // Run the time-dependent simulation.
  void
  run();

  // Setup system.
  void
  setup();

  // Assemble system.
  void
  assemble();

  // Solve system.
  void
  solve();

  // Output results.
  void
  output();

  // Benchmark quantities.
  double compute_reynolds_number() const;
  std::pair<double, double> compute_drag_lift_forces() const;
  double compute_pressure_difference() const;
  double compute_recirculation_length() const;
  double compute_strouhal_number() const;
  double compute_delta_p_at_half_period() const;

  struct BenchmarkResult
  {
    double Re     = 0.0;
    double cD     = 0.0;
    double cL     = 0.0;
    double dP     = 0.0;
    double La     = std::numeric_limits<double>::quiet_NaN();
    double cD_max = 0.0;
    double cL_max = 0.0;
    double St     = std::numeric_limits<double>::quiet_NaN();
    double dP_half = std::numeric_limits<double>::quiet_NaN();
  };

  BenchmarkResult compute_benchmark_result() const;
  void print_benchmark_quantities(const BenchmarkResult &res) const;

protected:

  // MPI parallel. /////////////////////////////////////////////////////////////

  // Number of MPI processes.
  const unsigned int mpi_size;

  // This MPI process.
  const unsigned int mpi_rank;

  // Parallel output stream.
  ConditionalOStream pcout;

  double cD_max = 0.0;
  double cL_max = 0.0;

  // Time series of drag and lift coefficients and pressure difference.
  std::vector<double> time_history;
  std::vector<double> cL_history;
  std::vector<double> cD_history;
  std::vector<double> dP_history;


  // Discretization. ///////////////////////////////////////////////////////////

  // Mesh file name.
  const std::string mesh_file_name;

  // Polynomial degree used for velocity.
  const unsigned int degree_velocity;

  // Polynomial degree used for pressure.
  const unsigned int degree_pressure;
  // Final time.
  const double T;

  // Theta parameter for the theta method.
  const double theta;

  // Time step.
  const double delta_t;

  // Data for the problem.
  ProblemData data;

  // Current time.
  double time = 0.0;

  // Current timestep number.
  unsigned int timestep_number = 0;

  // Mesh.
  parallel::fullydistributed::Triangulation<dim> mesh;

  // Finite element space.
  std::unique_ptr<FiniteElement<dim>> fe;

  // Quadrature formula.
  std::unique_ptr<Quadrature<dim>> quadrature;

  // Quadrature formula for face integrals.
  std::unique_ptr<Quadrature<dim - 1>> quadrature_face;

  // DoF handler.
  DoFHandler<dim> dof_handler;

  // DoFs owned by current process.
  IndexSet locally_owned_dofs;

  // DoFs owned by current process in the velocity and pressure blocks.
  std::vector<IndexSet> block_owned_dofs;

  // DoFs relevant to the current process (including ghost DoFs).
  IndexSet locally_relevant_dofs;

  // DoFs relevant to current process in the velocity and pressure blocks.
  std::vector<IndexSet> block_relevant_dofs;

  // System matrix.
  TrilinosWrappers::BlockSparseMatrix system_matrix;

  // Pressure mass matrix, needed for preconditioning. We use a block matrix for
  // convenience, but in practice we only look at the pressure-pressure block.
  TrilinosWrappers::BlockSparseMatrix pressure_mass;

  // Right-hand side vector in the linear system.
  TrilinosWrappers::MPI::BlockVector system_rhs;

  // System solution (without ghost elements).
  TrilinosWrappers::MPI::BlockVector solution_owned;

  // System solution (including ghost elements).
  TrilinosWrappers::MPI::BlockVector solution;

  // Physical constants for the benchmark.
  static constexpr double D_cylinder = 0.1;

  // 2D benchmark points.
  static constexpr double x_front_cylinder = 0.15;
  static constexpr double y_probe          = 0.20;
  static constexpr double x_back_cylinder  = 0.25;

  // Search interval for the recirculation zone.
  static constexpr double x_wake_start = x_back_cylinder;
  static constexpr double x_wake_end   = 2.20;

  double reference_velocity() const
  {
    // From the benchmark: U = 2 U(0,H/2,t) / 3.
    return 2.0 * data.Um / 3.0;
  }

};

#endif