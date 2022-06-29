#ifndef COLLISIONS_RAMBO_H
#define COLLISIONS_RAMBO_H

#include "MFP_collisions.H"

class CollisionsRambo : public Collisions {
public:
  CollisionsRambo();
  CollisionsRambo(const int idx, const sol::table &def);
  ~CollisionsRambo();

  virtual std::string get_tag() const override { return tag; }
  virtual ActionType get_type() const override {
    return ActionType::Collisions;
  }

  virtual void calc_update() override;

  Vector<Vector<Real>> ccs;

  static bool registered;
  static std::string tag;

  Real n0, Debye, c1, c2, c3, lnC;
};

#endif // COLLISIONS_RAMBO_H
