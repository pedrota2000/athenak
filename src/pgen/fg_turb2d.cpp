//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file fg_turb2d.cpp
//! \brief Problem generator for 2D driven turbulence in a (nearly) conformal special
//! relativistic fluid.  Written to produce the data set used to build the fluid/gravity
//! metric coefficients, which requires the four-velocity u^\mu and the fluid-frame
//! energy density e (the eigenvalue of T^{\mu\nu}: T = diag(-e,p,p,p) in the rest frame).
//!
//! CONFORMAL EOS.  AthenaK implements an ideal gas, p = (\Gamma-1)*e_int, where the
//! fluid-frame energy density is e = \rho + e_int (\rho = rest-mass density).  Hence
//!            T^\mu_\mu = -e + (d-1) p = -\rho + [(d-1)(\Gamma-1) - 1] e_int .
//! Choosing \Gamma = d/(d-1) (i.e. \Gamma = 3/2 for d = 3, 4/3 for d = 4) kills the
//! e_int term identically, leaving
//!            T^\mu_\mu = -\rho   exactly,   and   p/e = 1/[(d-1)(1 + 1/\epsilon)],
//! so the conformal (traceless) limit is approached as \epsilon = e_int/\rho -> \infty.
//! The run is therefore initialized "hot": T0 = p/\rho >> 1, which makes the residual
//! trace anomaly \rho/e = 1/(1 + T0*(d-1)) small and *measurable* (\rho is output, so
//! the anomaly is known pointwise and can be subtracted or monitored).
//!
//! The rest-mass density plays no other dynamical role in this limit; it is advected as
//! \partial_\mu (\rho u^\mu) = 0 and is used only as the conformal-anomaly diagnostic.
//!
//! Requires <coord>/special_rel = true, a 2D mesh (nx3 = 1), and (for the inverse
//! cascade) a <turb_driving> block with driving_type = 1, which is the only driving
//! type in turb_driver.cpp that produces a strictly in-plane (F_z = 0), exactly
//! divergence-free force on a 2D mesh.  See inputs/fluid_gravity_runs/README.md.

#include <iostream>
#include <string>
#include <cmath>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "coordinates/cell_locations.hpp"
#include "eos/eos.hpp"
#include "eos/ideal_c2p_hyd.hpp"
#include "hydro/hydro.hpp"
#include "pgen.hpp"

// user-defined history function
void FGTurb2DHistory(HistoryData *pdata, Mesh *pm);

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::UserProblem()
//! \brief uniform, hot, static conformal fluid; the turbulence driver does the rest

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;

  // enroll user history function (also on restarts)
  user_hist_func = FGTurb2DHistory;
  if (restart) return;

  // sanity checks -----------------------------------------------------------------
  if (pmbp->phydro == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "fg_turb2d requires a <hydro> block" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (!(pmbp->pcoord->is_special_relativistic)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "fg_turb2d requires <coord>/special_rel = true" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (!(pmy_mesh_->two_d)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "fg_turb2d is a 2D problem: set <mesh>/nx3 = 1 and nx2 > 1" << std::endl;
    exit(EXIT_FAILURE);
  }

  // parameters --------------------------------------------------------------------
  // d = number of space-time dimensions the conformal EOS refers to.  d=3 (2+1) gives
  // p = e/2 and Gamma = 3/2; d=4 (3+1) gives p = e/3 and Gamma = 4/3.
  Real dst  = pin->GetOrAddReal("problem", "d_spacetime", 3.0);
  Real dens0 = pin->GetOrAddReal("problem", "dens0", 1.0);
  Real temp0 = pin->GetOrAddReal("problem", "temp0", 100.0);   // T = p/rho

  EOS_Data &eos = pmbp->phydro->peos->eos_data;
  Real gm1 = eos.gamma - 1.0;
  Real gamma_conformal = dst/(dst - 1.0);
  if (std::fabs(eos.gamma - gamma_conformal) > 1.0e-8) {
    std::cout << "### WARNING in " << __FILE__ << ": <hydro>/gamma = " << eos.gamma
              << " is NOT the conformal value d/(d-1) = " << gamma_conformal
              << " for <problem>/d_spacetime = " << dst << ".  The EOS is not conformal."
              << std::endl;
  }
  Real p0 = dens0*temp0;
  Real eint0 = p0/gm1;
  Real anomaly = dens0/(dens0 + eint0);   // |T^mu_mu|/e at t=0
  if (global_variable::my_rank == 0) {
    std::cout << "fg_turb2d: rho0 = " << dens0 << ", p0 = " << p0 << ", e0 = "
              << (dens0 + eint0) << ", p/e = " << p0/(dens0 + eint0)
              << " (conformal value " << 1.0/(dst - 1.0) << "), trace anomaly rho/e = "
              << anomaly << std::endl;
  }

  // initial conditions ------------------------------------------------------------
  auto &indcs = pmy_mesh_->mb_indcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  auto &u0 = pmbp->phydro->u0;
  Real gam = eos.gamma;

  par_for("pgen_fg_turb2d", DevExeSpace(), 0, (pmbp->nmb_thispack-1), ks, ke, js, je,
  is, ie, KOKKOS_LAMBDA(int m, int k, int j, int i) {
    // static, uniform, hot fluid.  In SR the primitives are (rho, u^x, u^y, u^z, e_int)
    HydPrim1D w;
    w.d  = dens0;
    w.vx = 0.0;
    w.vy = 0.0;
    w.vz = 0.0;
    w.e  = eint0;
    // conserved vars are (D, m^i, E-D); use the EOS routine so the SR conversion (and
    // in particular the "evolve E-D" convention) cannot be gotten wrong here
    HydCons1D u;
    SingleP2C_IdealSRHyd(w, gam, u);
    u0(m,IDN,k,j,i) = u.d;
    u0(m,IM1,k,j,i) = u.mx;
    u0(m,IM2,k,j,i) = u.my;
    u0(m,IM3,k,j,i) = u.mz;
    u0(m,IEN,k,j,i) = u.e;
  });

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FGTurb2DHistory()
//! \brief volume-averaged diagnostics for the conformal 2D turbulence runs.  All entries
//! are true volume averages (the sums are divided by the box volume).

void FGTurb2DHistory(HistoryData *pdata, Mesh *pm) {
  pdata->nhist = 10;
  pdata->label[0] = "rho";      // <rho>              rest-mass density
  pdata->label[1] = "e";        // <e> = <rho+eint>   fluid-frame energy density
  pdata->label[2] = "pres";     // <p>
  pdata->label[3] = "p/e";      // <p>/<e>            = 1/(d-1) for a conformal fluid
  pdata->label[4] = "anom";     // <rho>/<e>          |T^mu_mu|/e, the trace anomaly
  pdata->label[5] = "usq";      // <u_i u^i>          4-velocity^2 (relativistic-ness)
  pdata->label[6] = "W";        // <u^t>              Lorentz factor
  pdata->label[7] = "Ekin";     // <(e+p) u_i u^i>    energy in the flow
  pdata->label[8] = "Etot";     // <E> = <(E-D) + D>  conserved lab-frame energy
  pdata->label[9] = "omsq";     // <omega_z^2>        enstrophy of the 3-velocity

  auto &w0_ = pm->pmb_pack->phydro->w0;
  auto &u0_ = pm->pmb_pack->phydro->u0;
  auto &size = pm->pmb_pack->pmb->mb_size;
  Real gm1 = pm->pmb_pack->phydro->peos->eos_data.gamma - 1.0;
  int &nhist_ = pdata->nhist;

  auto &indcs = pm->pmb_pack->pmesh->mb_indcs;
  int is = indcs.is; int nx1 = indcs.nx1;
  int js = indcs.js; int nx2 = indcs.nx2;
  int ks = indcs.ks; int nx3 = indcs.nx3;
  const int nmkji = (pm->pmb_pack->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji  = nx2*nx1;

  array_sum::GlobalSum sum_this_mb;
  Kokkos::parallel_reduce("FGHistSums", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, array_sum::GlobalSum &mb_sum) {
    int m = (idx)/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;
    Real vol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;

    Real rho  = w0_(m,IDN,k,j,i);
    Real eint = w0_(m,IEN,k,j,i);
    Real ux   = w0_(m,IVX,k,j,i);
    Real uy   = w0_(m,IVY,k,j,i);
    Real uz   = w0_(m,IVZ,k,j,i);
    Real usq  = ux*ux + uy*uy + uz*uz;
    Real ut   = sqrt(1.0 + usq);
    Real pres = gm1*eint;
    Real ener = rho + eint;

    // in-plane vorticity of the 3-velocity v^i = u^i/u^t
    Real utp = sqrt(1.0 + SQR(w0_(m,IVX,k,j,i+1)) + SQR(w0_(m,IVY,k,j,i+1))
                        + SQR(w0_(m,IVZ,k,j,i+1)));
    Real utm = sqrt(1.0 + SQR(w0_(m,IVX,k,j,i-1)) + SQR(w0_(m,IVY,k,j,i-1))
                        + SQR(w0_(m,IVZ,k,j,i-1)));
    Real dvydx = 0.5*(w0_(m,IVY,k,j,i+1)/utp - w0_(m,IVY,k,j,i-1)/utm)
                 /size.d_view(m).dx1;
    utp = sqrt(1.0 + SQR(w0_(m,IVX,k,j+1,i)) + SQR(w0_(m,IVY,k,j+1,i))
                   + SQR(w0_(m,IVZ,k,j+1,i)));
    utm = sqrt(1.0 + SQR(w0_(m,IVX,k,j-1,i)) + SQR(w0_(m,IVY,k,j-1,i))
                   + SQR(w0_(m,IVZ,k,j-1,i)));
    Real dvxdy = 0.5*(w0_(m,IVX,k,j+1,i)/utp - w0_(m,IVX,k,j-1,i)/utm)
                 /size.d_view(m).dx2;

    array_sum::GlobalSum hvars;
    hvars.the_array[0] = vol*rho;
    hvars.the_array[1] = vol*ener;
    hvars.the_array[2] = vol*pres;
    hvars.the_array[3] = vol*pres/ener;
    hvars.the_array[4] = vol*rho/ener;
    hvars.the_array[5] = vol*usq;
    hvars.the_array[6] = vol*ut;
    hvars.the_array[7] = vol*(ener + pres)*usq;
    hvars.the_array[8] = vol*(u0_(m,IEN,k,j,i) + u0_(m,IDN,k,j,i));
    hvars.the_array[9] = vol*SQR(dvydx - dvxdy);

    for (int n=nhist_; n<NHISTORY_VARIABLES; ++n) {
      hvars.the_array[n] = 0.0;
    }
    mb_sum += hvars;
  }, Kokkos::Sum<array_sum::GlobalSum>(sum_this_mb));
  Kokkos::fence();

  // divide by box volume to turn the integrals into averages
  Real tvol = (pm->mesh_size.x1max - pm->mesh_size.x1min)
             *(pm->mesh_size.x2max - pm->mesh_size.x2min)
             *(pm->mesh_size.x3max - pm->mesh_size.x3min);
  for (int n=0; n<pdata->nhist; ++n) {
    pdata->hdata[n] = sum_this_mb.the_array[n]/tvol;
  }
  return;
}
