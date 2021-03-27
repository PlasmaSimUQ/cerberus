
#include <AMReX_LevelBld.H>
#include <MFP.H>

using namespace amrex;

class MFPBld : public LevelBld {
  virtual void variableSetUp() override;
  virtual void variableCleanUp() override;
  virtual AmrLevel* operator()() override;
  virtual AmrLevel* operator()(Amr& papa, int lev, const Geometry& level_geom,
                               const BoxArray& ba,
                               const DistributionMapping& dm,
                               Real time) override;
};

MFPBld MFP_bld;

LevelBld* getLevelBld() { return &MFP_bld; }

void MFPBld::variableSetUp() { MFP::variableSetUp(); }

void MFPBld::variableCleanUp() { MFP::variableCleanUp(); }

AmrLevel* MFPBld::operator()() { return new MFP; }

AmrLevel* MFPBld::operator()(Amr& papa, int lev, const Geometry& level_geom,
                             const BoxArray& ba, const DistributionMapping& dm,
                             Real time) {
  return new MFP(papa, lev, level_geom, ba, dm, time);
}
