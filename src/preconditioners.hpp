#ifndef NAVIER_STOKES_PRECONDITIONERS_HPP
#define NAVIER_STOKES_PRECONDITIONERS_HPP

#include <deal.II/base/subscriptor.h>   
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_vector.h> 
#include <deal.II/lac/trilinos_solver.h>
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


    // Implementation of the SIMPLE preconditioner
    //   P_SIMPLE = [ F   0  ] [ I   D^-1 B^T ]
    //              [ B  -S~ ] [ 0   alpha I  ]
    // with D = diag(F), S~ = B D^-1 B^T
    class PreconditionSIMPLE : public NavierStokesPreconditionerBase {
        public:
            PreconditionSIMPLE() = default;

            // F_       : system_matrix.block(0,0)
            // Bt_      : system_matrix.block(0,1)   (this is B^T)
            // B_       : system_matrix.block(1,0) 
            // alpha_   : pressure damping parameter in (0,1], default 0.6      

            // Initialize the preconditioner:
            void
            initialize(const TrilinosWrappers::SparseMatrix &F_,
                       const TrilinosWrappers::SparseMatrix &Bt_,
                       const TrilinosWrappers::SparseMatrix &B_,
                       const TrilinosWrappers::MPI::BlockVector &vec_owned,
                       const double                          alpha_ = 0.6) 
            {
            F     = &F_;
            Bt    = &Bt_;
            B     = &B_;
            alpha = alpha_;
            
            // Approximate inverse of F via ILU (used inside an inner GMRES).
            preconditioner_F.initialize(F_);

            // D^{-1}: inverse of the diagonal of F.
            D_inv.reinit(vec_owned.block(0));
            for (const auto i : D_inv.locally_owned_elements())
                D_inv[i] = 1.0 / F_.diag_element(i);
            D_inv.compress(VectorOperation::insert);

            // Explicitly build the approximate (pseudo-Laplacian) Schur complement
            B_.mmult(S, *Bt, D_inv);

            preconditioner_S.initialize(S);
        }

            // Application of the preconditioner
            void
            vmult(TrilinosWrappers::MPI::BlockVector       &dst,
                  const TrilinosWrappers::MPI::BlockVector &src) const override
            {
                // Lower block-triangular solve: L z = src 
                tmp_u.reinit(src.block(0));
                SolverControl solver_control_u(2000, 1e-2 * src.block(0).l2_norm());
                SolverGMRES<TrilinosWrappers::MPI::Vector> solver_u(solver_control_u);
                solver_u.solve(*F, tmp_u, src.block(0), preconditioner_F);

                tmp_p.reinit(src.block(1));
                B->vmult(tmp_p, tmp_u);      // tmp_p = B * z_u
                tmp_p -= src.block(1);       // tmp_p = B z_u - r_p

                SolverControl solver_control_p(2000, 1e-2 * tmp_p.l2_norm());
                SolverCG<TrilinosWrappers::MPI::Vector> solver_p(solver_control_p);
                solver_p.solve(S, dst.block(1), tmp_p, preconditioner_S);

                // Upper block-triangular solve: U dst = z 
                dst.block(1) /= alpha;   // x_p = z_p / alpha

                tmp_u2.reinit(src.block(0));
                Bt->vmult(tmp_u2, dst.block(1));  // B^T x_p
                tmp_u2.scale(D_inv);              // D^{-1} B^T x_p
                dst.block(0)  = tmp_u;
                dst.block(0) -= tmp_u2;            // x_u = z_u - D^{-1} B^T x_p
                            
            }
        
        private: 
            const TrilinosWrappers::SparseMatrix *F;
            const TrilinosWrappers::SparseMatrix *B;
            const TrilinosWrappers::SparseMatrix *Bt;

            double alpha;

            TrilinosWrappers::MPI::Vector D_inv;
            TrilinosWrappers::SparseMatrix S;

            TrilinosWrappers::PreconditionILU preconditioner_F;
            TrilinosWrappers::PreconditionILU preconditioner_S;

            mutable TrilinosWrappers::MPI::Vector tmp_u;
            mutable TrilinosWrappers::MPI::Vector tmp_u2;
            mutable TrilinosWrappers::MPI::Vector tmp_p;
    };


    // Implementation of the Yosida preconditioner
    //   P_Y = [ F   0  ] [ I   F^-1 B^T ]
    //         [ B  -S~ ] [ 0     I      ]
    // with  S~ = Delta_t * B * diag(M_u)^-1 * B^T
    class PreconditionYosida : public NavierStokesPreconditionerBase {
        public:
            PreconditionYosida() = default;

            // F_            : system_matrix.block(0,0)
            // Bt_           : system_matrix.block(0,1)
            // B_            : system_matrix.block(1,0)
            // velocity_mass_: velocity mass matrix M_u (block(0,0) of a
            //                 separately assembled mass matrix)
            // dt_           : current timestep Delta_t

            // Initialize the preconditioner:
            void
            initialize(const TrilinosWrappers::SparseMatrix &F_,
                       const TrilinosWrappers::SparseMatrix &Bt_,
                       const TrilinosWrappers::SparseMatrix &B_,
                       const TrilinosWrappers::SparseMatrix &velocity_mass_,
                       const TrilinosWrappers::MPI::BlockVector &vec_owned,
                       const double                          dt_)
            {
                F     = &F_;
                Bt    = &Bt_;
                B     = &B_;
                dt    = dt_;

                // Preconditioner for velocity inversion
                preconditioner_F.initialize(F_);

                // Allocate and compute lumped inverse mass matrix
                Dmu_inv.reinit(vec_owned.block(0));
                for (const auto i : Dmu_inv.locally_owned_elements())
                {
                    double lumped = 0.0;
                    for (auto entry = velocity_mass_.begin(i); entry != velocity_mass_.end(i); ++entry)
                        lumped += std::abs(entry->value());
                    Dmu_inv[i] = 1.0 / lumped;
                }
                Dmu_inv.compress(VectorOperation::insert);

                // Build positive-definite Schur complement S = dt * B * D_Mu^-1 * B^T
                B_.mmult(schur_matrix, *Bt, Dmu_inv);
                schur_matrix *= dt;

                // S behaves like a standard elliptic Laplacian operator
                preconditioner_S.initialize(schur_matrix);
            }

            // Application of the preconditioner
            void
            vmult(TrilinosWrappers::MPI::BlockVector       &dst,
                  const TrilinosWrappers::MPI::BlockVector &src) const override
            {
                // Predictor velocity z_u = F^-1 r_u via GMRES
                tmp_u.reinit(src.block(0));
                SolverControl solver_control_F1(2000, 1e-2 * src.block(0).l2_norm());
                SolverGMRES<TrilinosWrappers::MPI::Vector> solver_gmres_F1(solver_control_F1);
                solver_gmres_F1.solve(*F, tmp_u, src.block(0), preconditioner_F);

                // Pressure RHS = B * z_u - r_p
                tmp_p.reinit(src.block(1));
                B->vmult(tmp_p, tmp_u);
                tmp_p -= src.block(1);

                // Solve for pressure dst_p = S^-1 RHS via Conjugate Gradient
                SolverControl solver_control_S(2000, 1e-2 * tmp_p.l2_norm());
                SolverCG<TrilinosWrappers::MPI::Vector> solver_cg_S(solver_control_S);
                solver_cg_S.solve(schur_matrix, dst.block(1), tmp_p, preconditioner_S);

                // Velocity Projection Correction: dst_u = z_u - F^-1 (B^T dst_p)
                tmp_u2.reinit(src.block(0));
                Bt->vmult(tmp_u2, dst.block(1)); 
                
                correction.reinit(src.block(0));
                SolverControl solver_control_F2(2000, 1e-2 * tmp_u2.l2_norm());
                SolverGMRES<TrilinosWrappers::MPI::Vector> solver_gmres_F2(solver_control_F2);
                solver_gmres_F2.solve(*F, correction, tmp_u2, preconditioner_F);

                dst.block(0) = tmp_u;
                dst.block(0) -= correction;

            }
        
        private: 
            const TrilinosWrappers::SparseMatrix *F;
            const TrilinosWrappers::SparseMatrix *B;
            const TrilinosWrappers::SparseMatrix *Bt;
            double dt;

            TrilinosWrappers::SparseMatrix schur_matrix;
            TrilinosWrappers::MPI::Vector  Dmu_inv;

            TrilinosWrappers::PreconditionILU preconditioner_F;
            TrilinosWrappers::PreconditionILU preconditioner_S;

            // Persistent storage to eliminate dynamic parallel vectors on the heap
            mutable TrilinosWrappers::MPI::Vector tmp_u;
            mutable TrilinosWrappers::MPI::Vector tmp_u2;
            mutable TrilinosWrappers::MPI::Vector correction;
            mutable TrilinosWrappers::MPI::Vector tmp_p;
    };


    // Implementation of the PCD preconditioner
    //   P_PCD = [ F   B^T             ]
    //           [ 0   -M_p F_p^-1 A_p ]
    class PreconditionPCD : public NavierStokesPreconditionerBase {
        public:
            PreconditionPCD() = default;

            // F_        : system_matrix.block(0,0)
            // B_        : system_matrix.block(1,0)
            // Bt_       : system_matrix.block(0,1)
            // Ap_       : pressure_laplacian
            // Cp_       : pressure_convection (pure convection term only)
            // Mp_       : pressure_mass.block(1,1)
            // nu_       : kinematic viscosity
            // dt_       : current timestep

            // Initialize the preconditioner:
            void
            initialize(const TrilinosWrappers::SparseMatrix &F_,
                       const TrilinosWrappers::SparseMatrix &B_,
                       const TrilinosWrappers::SparseMatrix &Bt_,
                       const TrilinosWrappers::SparseMatrix &Ap_,
                       const TrilinosWrappers::SparseMatrix &Cp_,
                       const TrilinosWrappers::SparseMatrix &Mp_,
                       const double                          nu_,
                       const double                          dt_)
            {
                F  = &F_;
                B  = &B_;
                Bt = &Bt_;
                Ap = &Ap_;
                Mp = &Mp_;

                // Assemble Fp = (1/dt) Mp + nu * Ap + Cp
                Fp_matrix.copy_from(Mp_);
                Fp_matrix *= 1.0 / dt_;
                Fp_matrix.add(nu_, Ap_);
                Fp_matrix.add(1.0, Cp_);
                Fp = &Fp_matrix;

                preconditioner_F.initialize(F_);
                preconditioner_Mp.initialize(Mp_);
                preconditioner_Ap.initialize(Ap_);
                }

            // Application of the preconditioner
            void
            vmult(TrilinosWrappers::MPI::BlockVector       &dst,
                  const TrilinosWrappers::MPI::BlockVector &src) const override
            {
                // w1 = Mp^-1 r_p via CG
                w1.reinit(src.block(1));
                SolverControl solver_control_Mp(2000, 1e-2 * src.block(1).l2_norm());
                SolverCG<TrilinosWrappers::MPI::Vector> solver_cg_Mp(solver_control_Mp);
                solver_cg_Mp.solve(*Mp, w1, src.block(1), preconditioner_Mp);

                // w2 = Fp * w1
                w2.reinit(src.block(1));
                Fp->vmult(w2, w1);

                // w3 = Ap^-1 w2 via CG
                w3.reinit(src.block(1));
                SolverControl solver_control_Ap(2000, 1e-2 * w2.l2_norm());
                SolverCG<TrilinosWrappers::MPI::Vector> solver_cg_Ap(solver_control_Ap);
                solver_cg_Ap.solve(*Ap, w3, w2, preconditioner_Ap);

                // y_p = -w3 
                dst.block(1) = w3;
                dst.block(1) *= -1.0;

                // y_u = F^-1 (r_u - B^T y_p)
                rhs_u.reinit(src.block(0));
                Bt->vmult(rhs_u, dst.block(1));
                
                // rhs_u = r_u - B^T y_p
                rhs_u.sadd(-1.0, src.block(0)); 

                SolverControl solver_control_F(2000, 1e-2 * rhs_u.l2_norm());
                SolverGMRES<TrilinosWrappers::MPI::Vector> solver_gmres_F(solver_control_F);
                solver_gmres_F.solve(*F, dst.block(0), rhs_u, preconditioner_F);
            }
        
        private: 
            const TrilinosWrappers::SparseMatrix *F;
            const TrilinosWrappers::SparseMatrix *B;
            const TrilinosWrappers::SparseMatrix *Bt;
            const TrilinosWrappers::SparseMatrix *Ap;
            const TrilinosWrappers::SparseMatrix *Fp;
            const TrilinosWrappers::SparseMatrix *Mp;

            TrilinosWrappers::PreconditionILU preconditioner_F;
            TrilinosWrappers::PreconditionILU preconditioner_Mp;
            TrilinosWrappers::PreconditionAMG preconditioner_Ap;
            TrilinosWrappers::SparseMatrix Fp_matrix;

            mutable TrilinosWrappers::MPI::Vector w1;
            mutable TrilinosWrappers::MPI::Vector w2;
            mutable TrilinosWrappers::MPI::Vector w3;
            mutable TrilinosWrappers::MPI::Vector rhs_u;
    };

}

#endif
