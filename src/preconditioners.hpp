#ifndef NAVIER_STOKES_PRECONDITIONERS_HPP
#define NAVIER_STOKES_PRECONDITIONERS_HPP

#include <deal.II/base/subscriptor.h>   
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_vector.h> 
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_cg.h>

#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>

using namespace dealii;

namespace NSPreconditioners {

    // Common abstract base class to enable runtime polymorphism
    class NavierStokesPreconditionerBase : public Subscriptor {
        public:
            virtual ~NavierStokesPreconditionerBase() = default;

            virtual void vmult(TrilinosWrappers::MPI::BlockVector       &dst, 
                               const TrilinosWrappers::MPI::BlockVector &src) const = 0;
    };

    // Implementation of the Identity preconditioner
    class PreconditionIdentity : public NavierStokesPreconditionerBase
    {
    public:
        PreconditionIdentity() = default;

        virtual void vmult(TrilinosWrappers::MPI::BlockVector       &dst,
                           const TrilinosWrappers::MPI::BlockVector &src) const override
        {
            dst = src;
        }
    };


    // Implementation of the block-diagonal preconditioner
    class PreconditionBlockDiagonal : public NavierStokesPreconditionerBase
    {
    public:
        PreconditionBlockDiagonal() = default;

        // Initialize the preconditioner
        void initialize(const TrilinosWrappers::SparseMatrix &velocity_stiffness_,
                        const TrilinosWrappers::SparseMatrix &pressure_mass_)
        {
            velocity_stiffness = &velocity_stiffness_;
            pressure_mass      = &pressure_mass_;

            preconditioner_velocity.initialize(velocity_stiffness_);
            preconditioner_pressure.initialize(pressure_mass_);
        }

        // Application of the preconditioner
        virtual void vmult(TrilinosWrappers::MPI::BlockVector       &dst,
                           const TrilinosWrappers::MPI::BlockVector &src) const override
        {
            dst.block(0) = 0.0;
            dst.block(1) = 0.0;

            SolverControl solver_control_velocity(2000, 1e-2 * src.block(0).l2_norm());
            SolverGMRES<TrilinosWrappers::MPI::Vector> solver_gmres_velocity(solver_control_velocity);
            
            solver_gmres_velocity.solve(*velocity_stiffness,
                                        dst.block(0),
                                        src.block(0),
                                        preconditioner_velocity);

            SolverControl solver_control_pressure(2000, 1e-2 * src.block(1).l2_norm());
            SolverCG<TrilinosWrappers::MPI::Vector> solver_cg_pressure(solver_control_pressure);
            
            solver_cg_pressure.solve(*pressure_mass,
                                     dst.block(1),
                                     src.block(1),
                                     preconditioner_pressure);
        }

    protected:
        // Velocity stiffness matrix
        const TrilinosWrappers::SparseMatrix *velocity_stiffness;
        // Preconditioner used for the velocity block
        TrilinosWrappers::PreconditionILU preconditioner_velocity;
        // Pressure mass matrix
        const TrilinosWrappers::SparseMatrix *pressure_mass;
        // Preconditioner used for the pressure block
        TrilinosWrappers::PreconditionILU preconditioner_pressure;
    };


    // Implementation of the block-triangular preconditioner
    class PreconditionBlockTriangular : public NavierStokesPreconditionerBase
    {
    public:
        PreconditionBlockTriangular() = default;

        // Initialize the preconditioner
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

            tmp.reinit(pressure_mass_.locally_owned_domain_indices(), pressure_mass_.get_mpi_communicator());
        }

        // Application of the preconditioner
        virtual void
        vmult(TrilinosWrappers::MPI::BlockVector       &dst,
              const TrilinosWrappers::MPI::BlockVector &src) const override
        {
            dst.block(0) = 0.0;
            dst.block(1) = 0.0;

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
        // Velocity stiffness matrix
        const TrilinosWrappers::SparseMatrix *velocity_stiffness;

        // Preconditioner used for the velocity block
        TrilinosWrappers::PreconditionILU preconditioner_velocity;

        // Pressure mass matrix
        const TrilinosWrappers::SparseMatrix *pressure_mass;

        // Preconditioner used for the pressure block
        TrilinosWrappers::PreconditionILU preconditioner_pressure;

        // B matrix
        const TrilinosWrappers::SparseMatrix *B;

        // Temporary vector
        mutable TrilinosWrappers::MPI::Vector tmp;
    };

    
    // Implementation of the SIMPLE preconditioner
    class PreconditionSIMPLE : public NavierStokesPreconditionerBase {
        public:
            PreconditionSIMPLE() {}

            // Initialize the preconditioner:
            void initialize(const TrilinosWrappers::BlockSparseMatrix &matrix,
                            const TrilinosWrappers::PreconditionBase  &F_inv_prec,
                            const TrilinosWrappers::PreconditionBase  &S_inv_prec) {
                system_matrix = &matrix;
                prec_F_inv = &F_inv_prec;
                prec_S_inv = &S_inv_prec;

                // Temporary vectors for block operations
                tmp_u.reinit(matrix.block(0,0).locally_owned_domain_indices(), matrix.block(0,0).get_mpi_communicator());
                tmp_u_2.reinit(matrix.block(0,0).locally_owned_domain_indices(), matrix.block(0,0).get_mpi_communicator());
                tmp_p.reinit(matrix.block(1,1).locally_owned_domain_indices(), matrix.block(1,1).get_mpi_communicator());
            }

            // Application of the preconditioner
            virtual void vmult(TrilinosWrappers::MPI::BlockVector       &dst, 
                               const TrilinosWrappers::MPI::BlockVector &src) const override {
                dst.block(0) = 0.0;
                dst.block(1) = 0.0;

                // Solve F * u_tmp = src_u
                prec_F_inv->vmult(dst.block(0), src.block(0));

                // Compute rhs_p = src_p - B * u_tmp
                system_matrix->block(1, 0).vmult(tmp_p, dst.block(0));
                tmp_p.sadd(-1.0, 1.0, src.block(1));

                // Solve S * dst_p = tmp_p
                prec_S_inv->vmult(dst.block(1), tmp_p);

                // Update velocity: dst_u = dst_u - F^-1 * B^T * dst_p
                system_matrix->block(0, 1).vmult(tmp_u, dst.block(1));
                prec_F_inv->vmult(tmp_u_2, tmp_u);
                dst.block(0).add(-1.0, tmp_u_2);
            }
        
        private: 
            const TrilinosWrappers::BlockSparseMatrix *system_matrix;
            const TrilinosWrappers::PreconditionBase  *prec_F_inv;
            const TrilinosWrappers::PreconditionBase  *prec_S_inv;

            mutable TrilinosWrappers::MPI::Vector tmp_u;
            mutable TrilinosWrappers::MPI::Vector tmp_u_2;
            mutable TrilinosWrappers::MPI::Vector tmp_p;
    };


    // Implementation of the Yosida preconditioner
    class PreconditionYosida : public NavierStokesPreconditionerBase {
        public:
            PreconditionYosida() {}

            // Initialize the preconditioner:
            void initialize(const TrilinosWrappers::BlockSparseMatrix &matrix,
                            const TrilinosWrappers::PreconditionBase  &F_inv_prec,
                            const TrilinosWrappers::PreconditionBase  &S_inv_prec) {
                system_matrix = &matrix;
                prec_F_inv = &F_inv_prec;
                prec_S_inv = &S_inv_prec;

                // Temporary vectors for block operations
                tmp_u.reinit(matrix.block(0,0).locally_owned_domain_indices(), matrix.block(0,0).get_mpi_communicator());
                tmp_u_2.reinit(matrix.block(0,0).locally_owned_domain_indices(), matrix.block(0,0).get_mpi_communicator());
                tmp_p.reinit(matrix.block(1,1).locally_owned_domain_indices(), matrix.block(1,1).get_mpi_communicator());
            }

            // Application of the preconditioner
            virtual void vmult(TrilinosWrappers::MPI::BlockVector       &dst, 
                               const TrilinosWrappers::MPI::BlockVector &src) const override {
                // Solve F * dst_u = src_u
                prec_F_inv->vmult(dst.block(0), src.block(0));

                // Compute rhs_p = src_p - B * dst_u
                system_matrix->block(1, 0).vmult(tmp_p, dst.block(0));
                tmp_p.sadd(-1.0, 1.0, src.block(1));

                // Solve S * dst_p = tmp_p
                prec_S_inv->vmult(dst.block(1), tmp_p);

                // Update velocity: dst_u = dst_u - F^-1 * B^T * dst_p
                system_matrix->block(0, 1).vmult(tmp_u, dst.block(1));
                prec_F_inv->vmult(tmp_u_2, tmp_u);
                dst.block(0).add(-1.0, tmp_u_2);
            }
        
        private: 
            const TrilinosWrappers::BlockSparseMatrix *system_matrix;
            const TrilinosWrappers::PreconditionBase  *prec_F_inv;
            const TrilinosWrappers::PreconditionBase  *prec_S_inv;

            mutable TrilinosWrappers::MPI::Vector tmp_u;
            mutable TrilinosWrappers::MPI::Vector tmp_u_2;
            mutable TrilinosWrappers::MPI::Vector tmp_p;
    };


    // Implementation of the PCD preconditioner
    class PreconditionPCD : public NavierStokesPreconditionerBase {
        public:
            PreconditionPCD() {}

            // Initialize the preconditioner:
            void initialize(const TrilinosWrappers::BlockSparseMatrix &matrix,       // full Navier-Stokes system matrix
                            const TrilinosWrappers::PreconditionBase  &F_inv_prec,   // preconditioner for the momentum block F
                            const TrilinosWrappers::PreconditionBase  &Mp_inv_prec,  // preconditioner for the pressure mass matrix Mp
                            const TrilinosWrappers::PreconditionBase  &Ap_inv_prec,  // preconditioner for the pressure Laplacian Ap
                            const TrilinosWrappers::SparseMatrix      &Fp_matrix) {  // pressure-space convection-diffusion operator Fp
                system_matrix = &matrix;
                prec_F_inv = &F_inv_prec;
                prec_Mp_inv = &Mp_inv_prec;
                prec_Ap_inv = &Ap_inv_prec;
                matrix_Fp = &Fp_matrix;

                // Temporary vectors for block operations
                tmp_u.reinit(matrix.block(0,0).locally_owned_domain_indices(), matrix.block(0,0).get_mpi_communicator());
                tmp_u2.reinit(matrix.block(0,0).locally_owned_domain_indices(), matrix.block(0,0).get_mpi_communicator());
                // tmp_p1.reinit(matrix.block(1,1).locally_owned_domain_indices(), matrix.block(1,1).get_mpi_communicator());
                // tmp_p2.reinit(matrix.block(1,1).locally_owned_domain_indices(), matrix.block(1,1).get_mpi_communicator());
                tmp_p1.reinit(Fp_matrix.locally_owned_range_indices(), Fp_matrix.get_mpi_communicator());
                tmp_p2.reinit(Fp_matrix.locally_owned_range_indices(), Fp_matrix.get_mpi_communicator());
            }

            // Application of the preconditioner
            virtual void vmult(TrilinosWrappers::MPI::BlockVector       &dst, 
                               const TrilinosWrappers::MPI::BlockVector &src) const override {
                // Solve F * dst_u = src_u
                prec_F_inv->vmult(dst.block(0), src.block(0));

                // Compute tmp_p1 = src_p - B * dst_u
                system_matrix->block(1, 0).vmult(tmp_p1, dst.block(0));
                tmp_p1.sadd(-1.0, 1.0, src.block(1));

                // Solve S * dst_p = tmp_p1
                // a) Apply Mp^-1
                prec_Mp_inv->vmult(tmp_p2, tmp_p1);
                // b) Apply Fp
                matrix_Fp->vmult(tmp_p1, tmp_p2);
                // c) Apply Ap^-1
                prec_Ap_inv->vmult(dst.block(1), tmp_p1);
                // d) Sign convention (Schur complement is negative)
                dst.block(1) *= -1.0;

                // Update velocity: dst_u = dst_u - F^-1 * B^T * dst_p
                system_matrix->block(0, 1).vmult(tmp_u, dst.block(1));
                prec_F_inv->vmult(tmp_u2, tmp_u);
                dst.block(0).add(-1.0, tmp_u2);
            }
        
        private: 
            const TrilinosWrappers::BlockSparseMatrix *system_matrix;
            const TrilinosWrappers::PreconditionBase  *prec_F_inv;
            const TrilinosWrappers::PreconditionBase  *prec_Mp_inv;
            const TrilinosWrappers::PreconditionBase  *prec_Ap_inv;
            const TrilinosWrappers::SparseMatrix      *matrix_Fp;

            mutable TrilinosWrappers::MPI::Vector tmp_u;
            mutable TrilinosWrappers::MPI::Vector tmp_u2;
            mutable TrilinosWrappers::MPI::Vector tmp_p1;
            mutable TrilinosWrappers::MPI::Vector tmp_p2;
    };


    // Implementation of the Additive Schwarz preconditioner
    class PreconditionAS : public NavierStokesPreconditionerBase {
        public:
            PreconditionAS() = default;

            // Accept both the system matrix and the specialized pressure matrix
            void initialize(const TrilinosWrappers::BlockSparseMatrix &matrix,
                            const TrilinosWrappers::SparseMatrix      &pressure_matrix,
                            const TrilinosWrappers::PreconditionILU::AdditionalData &data = 
                                  TrilinosWrappers::PreconditionILU::AdditionalData()) {
                
                as_preconditioners.resize(2);
                
                // Block 0 (Velocity): Uses the top-left fluid dynamics block (F)
                as_preconditioners[0].initialize(matrix.block(0, 0), data);
                
                // Block 1 (Pressure): Uses the valid pressure matrix instead of the zero block(1,1)
                as_preconditioners[1].initialize(pressure_matrix, data);
            }

            virtual void vmult(TrilinosWrappers::MPI::BlockVector       &dst, 
                               const TrilinosWrappers::MPI::BlockVector &src) const override {
                for (unsigned int i = 0; i < as_preconditioners.size(); ++i)
                {
                    as_preconditioners[i].vmult(dst.block(i), src.block(i));
                }
            }
        
        private: 
            std::vector<dealii::TrilinosWrappers::PreconditionILU> as_preconditioners;
    };
    
}

#endif
