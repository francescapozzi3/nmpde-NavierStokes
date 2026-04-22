#include "Stokes.hpp"

void
Stokes::setup()
{
  // Create the mesh.
  {
    pcout << "Initializing the mesh" << std::endl;

    Triangulation<dim> mesh_serial;

    GridIn<dim> grid_in;
    grid_in.attach_triangulation(mesh_serial);

    std::ifstream grid_in_file(mesh_file_name);
    grid_in.read_msh(grid_in_file);

    GridTools::partition_triangulation(mpi_size, mesh_serial);
    const auto construction_data = TriangulationDescription::Utilities::
      create_description_from_triangulation(mesh_serial, MPI_COMM_WORLD);
    mesh.create_triangulation(construction_data);

    pcout << "  Number of elements = " << mesh.n_global_active_cells()
          << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the finite element space.
  {
    pcout << "Initializing the finite element space" << std::endl;

    const FE_SimplexP<dim> fe_scalar_velocity(degree_velocity);
    const FE_SimplexP<dim> fe_scalar_pressure(degree_pressure);
    fe = std::make_unique<FESystem<dim>>(fe_scalar_velocity,
                                         dim,
                                         fe_scalar_pressure,
                                         1);

    pcout << "  Velocity degree:           = " << fe_scalar_velocity.degree
          << std::endl;
    pcout << "  Pressure degree:           = " << fe_scalar_pressure.degree
          << std::endl;
    pcout << "  DoFs per cell              = " << fe->dofs_per_cell
          << std::endl;

    quadrature = std::make_unique<QGaussSimplex<dim>>(fe->degree + 1);

    pcout << "  Quadrature points per cell = " << quadrature->size()
          << std::endl;

    quadrature_face = std::make_unique<QGaussSimplex<dim - 1>>(fe->degree + 1);

    pcout << "  Quadrature points per face = " << quadrature_face->size()
          << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the DoF handler.
  {
    pcout << "Initializing the DoF handler" << std::endl;

    dof_handler.reinit(mesh);
    dof_handler.distribute_dofs(*fe);

    // We want to reorder DoFs so that all velocity DoFs come first, and then
    // all pressure DoFs.
    std::vector<unsigned int> block_component(dim + 1, 0);
    block_component[dim] = 1;
    DoFRenumbering::component_wise(dof_handler, block_component);

    locally_owned_dofs = dof_handler.locally_owned_dofs();
    locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dof_handler);

    // Besides the locally owned and locally relevant indices for the whole
    // system (velocity and pressure), we will also need those for the
    // individual velocity and pressure blocks.
    std::vector<types::global_dof_index> dofs_per_block =
      DoFTools::count_dofs_per_fe_block(dof_handler, block_component);
    const unsigned int n_u = dofs_per_block[0];
    const unsigned int n_p = dofs_per_block[1];

    block_owned_dofs.resize(2);
    block_relevant_dofs.resize(2);
    block_owned_dofs[0]    = locally_owned_dofs.get_view(0, n_u);
    block_owned_dofs[1]    = locally_owned_dofs.get_view(n_u, n_u + n_p);
    block_relevant_dofs[0] = locally_relevant_dofs.get_view(0, n_u);
    block_relevant_dofs[1] = locally_relevant_dofs.get_view(n_u, n_u + n_p);

    pcout << "  Number of DoFs: " << std::endl;
    pcout << "    velocity = " << n_u << std::endl;
    pcout << "    pressure = " << n_p << std::endl;
    pcout << "    total    = " << n_u + n_p << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the linear system.
  {
    pcout << "Initializing the linear system" << std::endl;

    pcout << "  Initializing the sparsity pattern" << std::endl;

    // Velocity DoFs interact with other velocity DoFs (the weak formulation has
    // terms involving u times v), and pressure DoFs interact with velocity DoFs
    // (there are terms involving p times v or u times q). However, pressure
    // DoFs do not interact with other pressure DoFs (there are no terms
    // involving p times q). We build a table to store this information, so that
    // the sparsity pattern can be built accordingly.
    Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1);
    for (unsigned int c = 0; c < dim + 1; ++c)
      {
        for (unsigned int d = 0; d < dim + 1; ++d)
          {
            if (c == dim && d == dim) // pressure-pressure term
              coupling[c][d] = DoFTools::none;
            else // other combinations
              coupling[c][d] = DoFTools::always;
          }
      }

    TrilinosWrappers::BlockSparsityPattern sparsity(block_owned_dofs,
                                                    MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern(dof_handler, coupling, sparsity);
    sparsity.compress();

    // We also build a sparsity pattern for the pressure mass matrix.
    for (unsigned int c = 0; c < dim + 1; ++c)
      {
        for (unsigned int d = 0; d < dim + 1; ++d)
          {
            if (c == dim && d == dim) // pressure-pressure term
              coupling[c][d] = DoFTools::always;
            else // other combinations
              coupling[c][d] = DoFTools::none;
          }
      }
    TrilinosWrappers::BlockSparsityPattern sparsity_pressure_mass(
      block_owned_dofs, MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern(dof_handler,
                                    coupling,
                                    sparsity_pressure_mass);
    sparsity_pressure_mass.compress();

    pcout << "  Initializing the matrices" << std::endl;
    system_matrix.reinit(sparsity);
    pressure_mass.reinit(sparsity_pressure_mass);

    pcout << "  Initializing the system right-hand side" << std::endl;
    system_rhs.reinit(block_owned_dofs, MPI_COMM_WORLD);
    pcout << "  Initializing the solution vector" << std::endl;
    solution_owned.reinit(block_owned_dofs, MPI_COMM_WORLD);
    solution.reinit(block_owned_dofs, block_relevant_dofs, MPI_COMM_WORLD);
  }
}

void
Stokes::assemble()
{
  pcout << "===============================================" << std::endl;
  pcout << "Assembling the system" << std::endl;

  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q           = quadrature->size();
  const unsigned int n_q_face      = quadrature_face->size();

  FEValues<dim>     fe_values(*fe,
                          *quadrature,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);
  FEFaceValues<dim> fe_face_values(*fe,
                                   *quadrature_face,
                                   update_values | update_normal_vectors |
                                     update_JxW_values);

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> cell_pressure_mass_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  system_matrix = 0.0;
  system_rhs    = 0.0;
  pressure_mass = 0.0;

  FEValuesExtractors::Vector velocity(0);
  FEValuesExtractors::Scalar pressure(dim);
  // Evaluation of the old solution on quadrature nodes of current cell.
  std::vector<Tensor<1, dim>> solution_old_values(n_q);

  // Evaluation of the gradient of the old solution on quadrature nodes of
  // current cell.
  std::vector<Tensor<dim, dim>> solution_old_grads(n_q);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      fe_values.reinit(cell);

      cell_matrix               = 0.0;
      cell_rhs                  = 0.0;
      cell_pressure_mass_matrix = 0.0;

      // Evaluate the old solution and its gradient on quadrature nodes.
      fe_values[velocity].get_function_values(solution, solution_old_values);
      fe_values[velocity].get_function_gradients(solution, solution_old_grads);


      for (unsigned int q = 0; q < n_q; ++q)
        {

          
           const Tensor<1, dim> f_old_loc = f(fe_values.quadrature_point(q), time - delta_t);
           const Tensor<1, dim> f_new_loc = f(fe_values.quadrature_point(q), time);


          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {

                 // Time derivative: (1/dt) (u^{n+1}, v)
                  cell_matrix(i, j) += (1.0 / delta_t) *             //
                                       fe_values[velocity].value(i, q) * //
                                       fe_values[velocity].value(j, q) * //
                                       fe_values.JxW(q);

                  // Viscosity term: nu * theta * (grad u^{n+1}, grad v)
                  cell_matrix(i, j) +=
                    nu * theta *
                    scalar_product(fe_values[velocity].gradient(i, q),
                                   fe_values[velocity].gradient(j, q)) *
                    fe_values.JxW(q);


                  // Convection term: ((u^n · grad) u^{n+1}, v)
                  cell_matrix(i, j) += (fe_values[velocity].gradient(j, q) * 
                  solution_old_values[q]) * fe_values[velocity].value(i, q) * 
                  fe_values.JxW(q);

            


                 // Pressure term in the momentum equation: -(p^{n+1}, div v)
                  cell_matrix(i, j) -= fe_values[velocity].divergence(i, q) *
                                       fe_values[pressure].value(j, q) *
                                       fe_values.JxW(q);

                  // Pressure term in the continuity equation: -(div u^{n+1}, q)
                  cell_matrix(i, j) -= fe_values[velocity].divergence(j, q) *
                                       fe_values[pressure].value(i, q) *
                                       fe_values.JxW(q);


                  // Pressure mass matrix
                   cell_pressure_mass_matrix(i, j) += (nu + (1.0 / delta_t)) *
                    fe_values[pressure].value(i, q) *
                    fe_values[pressure].value(j, q) *
                    fe_values.JxW(q);
                 
               
                  }

                 // Time derivative: (1/dt)(u^n, v)
                  cell_rhs(i) += (1.0 / delta_t) *             //
                             fe_values[velocity].value(i, q) * //
                             solution_old_values[q] *      //
                             fe_values.JxW(q);

                  // - nu (1-theta) (grad u^n, grad v)
                  cell_rhs(i) -= (1.0 - theta) * nu *                   //
                             scalar_product(fe_values[velocity].gradient(i, q), //
                                            solution_old_grads[q]) *    //
                             fe_values.JxW(q);
                  
              // Forcing term: (f^{n+theta}, v)
              cell_rhs(i) +=  (theta * f_new_loc + (1.0 - theta) * f_old_loc) *                      //
                             fe_values[velocity].value(i, q) * //
                             fe_values.JxW(q);


            }
        }

      // Boundary integral for Neumann BCs.
      if (cell->at_boundary())
        {
          for (unsigned int f = 0; f < cell->n_faces(); ++f)
            {
              
              if (cell->face(f)->at_boundary() &&
                  (cell->face(f)->boundary_id() == 2 ))
                {
                  fe_face_values.reinit(cell, f);

                  for (unsigned int q = 0; q < n_q_face; ++q)
                    {
                    //const Tensor<1, dim> h_val = h(fe_face_values.quadrature_point(q), time);
                     
                    for (unsigned int i = 0; i < dofs_per_cell; ++i)
                        {
                          // mettere h_val al posto di h e quest0
                          //cell_rhs(i) += (h_val * fe_face_values[velocity].value(i, q)) * fe_face_values.JxW(q);

                          cell_rhs(i) +=
                            -h *
                            scalar_product(fe_face_values.normal_vector(q),
                                           fe_face_values[velocity].value(i, q)) *
                            fe_face_values.JxW(q);
                        }
                    }
                }
            }
        }

      cell->get_dof_indices(dof_indices);

      system_matrix.add(dof_indices, cell_matrix);
      system_rhs.add(dof_indices, cell_rhs);
      pressure_mass.add(dof_indices, cell_pressure_mass_matrix);
    }

  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);
  pressure_mass.compress(VectorOperation::add);

  // Dirichlet boundary conditions.
  {
    std::map<types::global_dof_index, double>           boundary_values;
    std::map<types::boundary_id, const Function<dim> *> boundary_functions;
    
    InletVelocity inlet;
    Functions::ZeroFunction<dim> zero_function(dim + 1);
 

    boundary_functions[1] = &inlet; // inlet
    boundary_functions[3] = &zero_function; // top/bottom walls
    boundary_functions[4] = &zero_function; // cylinder

   ComponentMask mask_velocity(dim + 1, true);
   mask_velocity.set(dim, false);

    VectorTools::interpolate_boundary_values(dof_handler,
                                             boundary_functions,
                                             boundary_values,
                                             mask_velocity);
    

    MatrixTools::apply_boundary_values(
      boundary_values, system_matrix, solution_owned, system_rhs, false);
  }
}

void
Stokes::solve()
{
  pcout << "===============================================" << std::endl;

  SolverControl solver_control(10000, 1e-4 * system_rhs.l2_norm());
  
  SolverGMRES<TrilinosWrappers::MPI::BlockVector> solver(solver_control);

 /* PreconditionBlockDiagonal preconditioner;
  preconditioner.initialize(system_matrix.block(0, 0),
                             pressure_mass.block(1, 1));
 */
 
 PreconditionBlockTriangular preconditioner;
  preconditioner.initialize(system_matrix.block(0, 0),
                            pressure_mass.block(1, 1),
                            system_matrix.block(1, 0));
  

  pcout << "Solving the linear system" << std::endl;
  solver.solve(system_matrix, solution_owned, system_rhs, preconditioner);
  pcout << "  " << solver_control.last_step() << " GMRES iterations"
        << std::endl;

  solution = solution_owned;
}

void
Stokes::output()
{
  pcout << "===============================================" << std::endl;

  DataOut<dim> data_out;

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    interpretation(dim,
                   DataComponentInterpretation::component_is_part_of_vector);
  interpretation.push_back(DataComponentInterpretation::component_is_scalar);

  std::vector<std::string> names(dim, "velocity");
  names.push_back("pressure");

  data_out.add_data_vector(dof_handler, solution, names, interpretation);

  std::vector<unsigned int> partition_int(mesh.n_active_cells());
  GridTools::get_subdomain_association(mesh, partition_int);
  const Vector<double> partitioning(partition_int.begin(), partition_int.end());
  data_out.add_data_vector(partitioning, "partitioning");

  data_out.build_patches();

  const std::string output_file_name = "output-stokes";
  data_out.write_vtu_with_pvtu_record(/* folder = */ "./",
                                      /* basename = */ output_file_name,
                                      /* index = */ timestep_number,
                                      MPI_COMM_WORLD);

  pcout << "Output written to " << output_file_name << std::endl;
  pcout << "===============================================" << std::endl;
}

void
Stokes::run()
{
  // Setup initial conditions.
  {
    setup();

    Functions::ZeroFunction<dim> zero_initial(dim + 1);
   VectorTools::interpolate(dof_handler, zero_initial, solution_owned);
   solution = solution_owned;


    time            = 0.0;
    timestep_number = 0;

    // Output initial condition.
    output();
  }

  pcout << "===============================================" << std::endl;

  // Time-stepping loop.
  while (time < T - 0.5 * delta_t)
    {
      time += delta_t;
      ++timestep_number;

      pcout << "Timestep " << std::setw(3) << timestep_number
            << ", time = " << std::setw(4) << std::fixed << std::setprecision(2)
            << time << " : ";

      assemble();
      solve();

      // Perform parallel communication to update the ghost values of the
      // solution vector.
      solution = solution_owned;

      output();
    }

// Compute benchmark quantities and write them to file.
      const double cD = compute_drag_coefficient();
      const double cL = compute_lift_coefficient();
      const double dP = compute_pressure_difference();
      const double La = compute_recirculation_length();

      std::ofstream out("benchmark_quantities.csv", std::ios::app);
      out << timestep_number  << ","
          << time << ","
          << cD << ","
          << cL << ","
          << La << ","
          << dP << '\n';

  print_benchmark_quantities();

}


std::pair<double, double>
Stokes::compute_drag_lift_forces() const
{
  FEFaceValues<dim> fe_face_values(*fe,
                                   *quadrature_face,
                                   update_values | update_gradients |
                                     update_normal_vectors |
                                     update_JxW_values);

  FEValuesExtractors::Vector velocity(0);
  FEValuesExtractors::Scalar pressure(dim);

  const unsigned int n_q_face = quadrature_face->size();

  std::vector<Tensor<1, dim>> velocity_values(n_q_face);
  std::vector<Tensor<2, dim>> velocity_grads(n_q_face);
  std::vector<double> pressure_values(n_q_face);

  double drag = 0.0;
  double lift = 0.0;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      if (!cell->at_boundary())
        continue;

      for (unsigned int f = 0; f < cell->n_faces(); ++f)
        {
          if (!(cell->face(f)->at_boundary() && cell->face(f)->boundary_id() == 4))
            continue;

          fe_face_values.reinit(cell, f);

          fe_face_values[velocity].get_function_values(solution, velocity_values);
          fe_face_values[velocity].get_function_gradients(solution, velocity_grads);
          fe_face_values[pressure].get_function_values(solution, pressure_values);

               
        for (unsigned int q = 0; q < n_q_face; ++q)
        {
          const Tensor<1, dim> n = fe_face_values.normal_vector(q);
          const double p = pressure_values[q];
          const Tensor<2, dim> grad_u = velocity_grads[q];

          // Force exerted by fluid on the cylinder
          // F = \int (p * n - rho * nu * grad_u * n) dS
          const double dF_D = (p * n[0] - rho * nu * (grad_u[0] * n)) * fe_face_values.JxW(q);
          const double dF_L = (p * n[1] - rho * nu * (grad_u[1] * n)) * fe_face_values.JxW(q);

          drag += dF_D;
          lift += dF_L;
        }

        }


    }

  // Sum contributions across MPI tasks.
  const double drag_global = Utilities::MPI::sum(drag, MPI_COMM_WORLD);
  const double lift_global = Utilities::MPI::sum(lift, MPI_COMM_WORLD);

  return {drag_global, lift_global};
}

double
Stokes::compute_drag_force() const
{
  return compute_drag_lift_forces().first;
}

double
Stokes::compute_lift_force() const
{
  return compute_drag_lift_forces().second;
}

double
Stokes::compute_drag_coefficient() const
{
  const double U = reference_velocity();
  return 2.0 * compute_drag_force() / (rho * U * U * D_cylinder);
}

double
Stokes::compute_lift_coefficient() const
{
  const double U = reference_velocity();
  return 2.0 * compute_lift_force() / (rho * U * U * D_cylinder);
}

double
Stokes::compute_pressure_difference() const
{
  // point_value() is the simplest option; use solution (ghosted vector).

  Vector<double> values(dim + 1);


  const Point<dim> p_front(x_front_cylinder, y_probe);
  const Point<dim> p_back(x_back_cylinder, y_probe);

  VectorTools::point_value(dof_handler, solution, p_front, values);
  const double p_a = values[dim];

  VectorTools::point_value(dof_handler, solution, p_back, values);
  const double p_e = values[dim];

  return p_a - p_e;
}


  double Stokes::compute_recirculation_length() const 
{ 
  // Standard numerical interpretation of x_r: 
  // first zero crossing of u_x along y = 0.2 downstream of the cylinder. 
  const unsigned int n_samples = 400; 
  const double y = y_probe; 
  auto u_x_at = [&](const double x) -> double { 
      Vector<double> values(dim + 1); 
      const Point<dim> p(x, y); 
      VectorTools::point_value(dof_handler, solution, p, values); 
      return values[0]; 
  }; 

  double x_prev = x_wake_start + 1e-6; 
  double u_prev = u_x_at(x_prev); 
  for (unsigned int k = 1; k <= n_samples; ++k) { 
      const double x = x_wake_start + (x_wake_end - x_wake_start) * static_cast<double>(k) / n_samples; 
      const double u = u_x_at(x); // Look for a sign change from negative to non-negative. 
      if (u_prev < 0.0 && u >= 0.0) 
      { 
          const double xr = x_prev - u_prev * (x - x_prev) / (u - u_prev); 
          return xr - x_back_cylinder; 
      } 
      x_prev = x; u_prev = u; 
  } 
  return std::numeric_limits<double>::quiet_NaN(); 
}

double
Stokes::compute_reynolds_number() const
{
  const double U = reference_velocity();

  return U * D_cylinder / nu;
}

void
Stokes::print_benchmark_quantities() const
{
  pcout << std::fixed << std::setprecision(4);
  
  pcout << "Benchmark quantities at t = " << time << '\n';
  pcout << "  Re = " << compute_reynolds_number() << '\n';
  pcout << "  cD = " << compute_drag_coefficient() << '\n';
  pcout << "  cL = " << compute_lift_coefficient() << '\n';
  pcout << "  ΔP = " << compute_pressure_difference() << '\n';
  pcout << "  La = " << compute_recirculation_length() << '\n';
}