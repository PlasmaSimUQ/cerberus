#ifndef MFP_INTERPOLATION_H
#define MFP_INTERPOLATION_H
#ifdef SYMPLECTIC

    #include <AMReX_GpuQualifiers.H>
    #include <AMReX_REAL.H>

AMREX_GPU_HOST_DEVICE amrex::Real W1(amrex::Real x);
AMREX_GPU_HOST_DEVICE amrex::Real Wp(amrex::Real x);
AMREX_GPU_HOST_DEVICE amrex::Real I_W1(amrex::Real a, amrex::Real b);
AMREX_GPU_HOST_DEVICE amrex::Real I_Wp(amrex::Real a, amrex::Real b);
extern const int interpolation_range;

#endif
#endif  // MFP_INTERPOLATION_H
