#include "MFP_bc.H"

using namespace amrex;
//
// Components are:
//  Interior,  Inflow,  Outflow,   Symmetry, SlipWall, NoSlipWall, Invert (only for scalar)
//
// Note: The component for inflow is given here as a outflow, this is because we handle
// the inflow condition ourselves (see State::update_face_prim)
//
int scalar_bc[] = {BCType::int_dir,      BCType::foextrap,
                          BCType::foextrap,     BCType::reflect_even,
                          BCType::reflect_even, BCType::reflect_even, BCType::reflect_odd};

int norm_vel_bc[] = {BCType::int_dir,     BCType::foextrap,
                            BCType::foextrap,    BCType::reflect_odd,
                            BCType::reflect_odd, BCType::reflect_odd};

int tang_vel_bc[] = {BCType::int_dir,      BCType::foextrap,
                            BCType::foextrap,     BCType::reflect_even,
                            BCType::reflect_even, BCType::reflect_odd};

void set_scalar_bc(BCRec& bc, const BCRec& phys_bc) {
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  for (int i = 0; i < AMREX_SPACEDIM; i++) {
    bc.setLo(i, scalar_bc[lo_bc[i]]);
    bc.setHi(i, scalar_bc[hi_bc[i]]);
  }
}

void set_x_vel_bc(BCRec& bc, const BCRec& phys_bc) {
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  bc.setLo(0, norm_vel_bc[lo_bc[0]]);
  bc.setHi(0, norm_vel_bc[hi_bc[0]]);
#if (AMREX_SPACEDIM >= 2)
  bc.setLo(1, tang_vel_bc[lo_bc[1]]);
  bc.setHi(1, tang_vel_bc[hi_bc[1]]);
#endif
#if (AMREX_SPACEDIM == 3)
  bc.setLo(2, tang_vel_bc[lo_bc[2]]);
  bc.setHi(2, tang_vel_bc[hi_bc[2]]);
#endif
}

void set_y_vel_bc(BCRec& bc, const BCRec& phys_bc) {
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  bc.setLo(0, tang_vel_bc[lo_bc[0]]);
  bc.setHi(0, tang_vel_bc[hi_bc[0]]);
#if (AMREX_SPACEDIM >= 2)
  bc.setLo(1, norm_vel_bc[lo_bc[1]]);
  bc.setHi(1, norm_vel_bc[hi_bc[1]]);
#endif
#if (AMREX_SPACEDIM == 3)
  bc.setLo(2, tang_vel_bc[lo_bc[2]]);
  bc.setHi(2, tang_vel_bc[hi_bc[2]]);
#endif
}

void set_z_vel_bc(BCRec& bc, const BCRec& phys_bc) {
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  bc.setLo(0, tang_vel_bc[lo_bc[0]]);
  bc.setHi(0, tang_vel_bc[hi_bc[0]]);
#if (AMREX_SPACEDIM >= 2)
  bc.setLo(1, tang_vel_bc[lo_bc[1]]);
  bc.setHi(1, tang_vel_bc[hi_bc[1]]);
#endif
#if (AMREX_SPACEDIM == 3)
  bc.setLo(2, norm_vel_bc[lo_bc[2]]);
  bc.setHi(2, norm_vel_bc[hi_bc[2]]);
#endif
}

// Field boundary types
// Components are:
//  Interior, Inflow, Edge, Symmetry, Anti-Symmetry
//
// Note: The component for inflow is given here as a outflow, this is because we handle
// the inflow condition ourselves (see State::update_face_prim)
//

int norm_E_field_bc[] = {BCType::int_dir, BCType::foextrap, BCType::foextrap,
                                BCType::reflect_odd, BCType::reflect_even};

int tang_E_field_bc[] = {BCType::int_dir, BCType::foextrap, BCType::foextrap,
                                BCType::reflect_even, BCType::reflect_odd};

int norm_B_field_bc[] = {BCType::int_dir, BCType::foextrap, BCType::foextrap,
                                BCType::reflect_even, BCType::reflect_odd};

int tang_B_field_bc[] = {BCType::int_dir, BCType::foextrap, BCType::foextrap,
                                BCType::reflect_odd, BCType::reflect_even};

void set_x_D_bc(BCRec& bc, const BCRec& phys_bc) {
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  bc.setLo(0, norm_E_field_bc[lo_bc[0]]);
  bc.setHi(0, norm_E_field_bc[hi_bc[0]]);
#if (AMREX_SPACEDIM >= 2)
  bc.setLo(1, tang_E_field_bc[lo_bc[1]]);
  bc.setHi(1, tang_E_field_bc[hi_bc[1]]);
#endif
#if (AMREX_SPACEDIM == 3)
  bc.setLo(2, tang_E_field_bc[lo_bc[2]]);
  bc.setHi(2, tang_E_field_bc[hi_bc[2]]);
#endif
}

void set_x_B_bc(BCRec& bc, const BCRec& phys_bc) {
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  bc.setLo(0, norm_B_field_bc[lo_bc[0]]);
  bc.setHi(0, norm_B_field_bc[hi_bc[0]]);
#if (AMREX_SPACEDIM >= 2)
  bc.setLo(1, tang_B_field_bc[lo_bc[1]]);
  bc.setHi(1, tang_B_field_bc[hi_bc[1]]);
#endif
#if (AMREX_SPACEDIM == 3)
  bc.setLo(2, tang_B_field_bc[lo_bc[2]]);
  bc.setHi(2, tang_B_field_bc[hi_bc[2]]);
#endif
}

void set_y_D_bc(BCRec& bc, const BCRec& phys_bc) {
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  bc.setLo(0, tang_E_field_bc[lo_bc[0]]);
  bc.setHi(0, tang_E_field_bc[hi_bc[0]]);
#if (AMREX_SPACEDIM >= 2)
  bc.setLo(1, norm_E_field_bc[lo_bc[1]]);
  bc.setHi(1, norm_E_field_bc[hi_bc[1]]);
#endif
#if (AMREX_SPACEDIM == 3)
  bc.setLo(2, tang_E_field_bc[lo_bc[2]]);
  bc.setHi(2, tang_E_field_bc[hi_bc[2]]);
#endif
}

void set_y_B_bc(BCRec& bc, const BCRec& phys_bc) {
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  bc.setLo(0, tang_B_field_bc[lo_bc[0]]);
  bc.setHi(0, tang_B_field_bc[hi_bc[0]]);
#if (AMREX_SPACEDIM >= 2)
  bc.setLo(1, norm_B_field_bc[lo_bc[1]]);
  bc.setHi(1, norm_B_field_bc[hi_bc[1]]);
#endif
#if (AMREX_SPACEDIM == 3)
  bc.setLo(2, tang_B_field_bc[lo_bc[2]]);
  bc.setHi(2, tang_B_field_bc[hi_bc[2]]);
#endif
}

void set_z_D_bc(BCRec& bc, const BCRec& phys_bc) {
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  bc.setLo(0, tang_E_field_bc[lo_bc[0]]);
  bc.setHi(0, tang_E_field_bc[hi_bc[0]]);
#if (AMREX_SPACEDIM >= 2)
  bc.setLo(1, tang_E_field_bc[lo_bc[1]]);
  bc.setHi(1, tang_E_field_bc[hi_bc[1]]);
#endif
#if (AMREX_SPACEDIM == 3)
  bc.setLo(2, norm_E_field_bc[lo_bc[2]]);
  bc.setHi(2, norm_E_field_bc[hi_bc[2]]);
#endif
}

void set_z_B_bc(BCRec& bc, const BCRec& phys_bc) {
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  bc.setLo(0, tang_B_field_bc[lo_bc[0]]);
  bc.setHi(0, tang_B_field_bc[hi_bc[0]]);
#if (AMREX_SPACEDIM >= 2)
  bc.setLo(1, tang_B_field_bc[lo_bc[1]]);
  bc.setHi(1, tang_B_field_bc[hi_bc[1]]);
#endif
#if (AMREX_SPACEDIM == 3)
  bc.setLo(2, norm_B_field_bc[lo_bc[2]]);
  bc.setHi(2, norm_B_field_bc[hi_bc[2]]);
#endif
}
