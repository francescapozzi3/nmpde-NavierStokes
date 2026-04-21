#include "Stokes.hpp"

// Main function.
int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  constexpr unsigned int dim = Stokes::dim;


  const std::string  mesh_file_name  = "../mesh/mesh-2D-cylinder-circular-cs-0.0205.msh";
  const unsigned int degree_velocity = 2;
  const unsigned int degree_pressure = 1;
  const double T = 1;
  const double theta = 1.0;
  const double delta_t = 0.01;


 // TEST 2D - 1 (steady): f = 0,  h = 0, nu = 0.001, u0 = 0 
    const auto f  = [](const Point<dim>  &/*p*/, const double  &/*t*/) {
    Tensor<1, dim> F;
    F = 0.0;
    return F;
  };

  Stokes problem(mesh_file_name, degree_velocity, degree_pressure, T, theta, delta_t, f);

  problem.run();

  return 0;
}