#ifndef MFP_BC_H
#define MFP_BC_H

#include <AMReX.H>
#include <AMReX_BCRec.H>


using namespace amrex;

typedef void (*set_bc)(BCRec& bc, const BCRec& phys_bc);


// Hydro boundary types

void set_scalar_bc(BCRec& bc, const BCRec& phys_bc);

void set_x_vel_bc(BCRec& bc, const BCRec& phys_bc);

void set_y_vel_bc(BCRec& bc, const BCRec& phys_bc);

void set_z_vel_bc(BCRec& bc, const BCRec& phys_bc);

// Field boundary types


void set_x_D_bc(BCRec& bc, const BCRec& phys_bc);

void set_x_B_bc(BCRec& bc, const BCRec& phys_bc);

void set_y_D_bc(BCRec& bc, const BCRec& phys_bc);

void set_y_B_bc(BCRec& bc, const BCRec& phys_bc);

void set_z_D_bc(BCRec& bc, const BCRec& phys_bc);

void set_z_B_bc(BCRec& bc, const BCRec& phys_bc);

#endif // MFP_BC_H
