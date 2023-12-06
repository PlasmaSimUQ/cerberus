#ifndef SYMPLECTIC_H
#define SYMPLECTIC_H

#ifdef SYMPLECTIC

    #include "MFP_action.H"
    #include "MFP_chargedparticle.H"
    #include "MFP_field.H"

    #define X 0
    #define Y 1
    #define Z 2

typedef void (*F_INTERIOR)(
  const Geometry, Box const&, Array4<Real const> const&, Array4<Real> const&, double);

template <F_INTERIOR InteriorF>
void push_ff(const Geometry geom,
             Box const& bx,
             Array4<Real const> const& Source,
             Array4<Real> const& Target,
             double dt)
{
    InteriorF(geom, bx, Source, Target, dt);
}

class Symplectic : public Action
{
  public:
    Symplectic();
    Symplectic(const int idx, const sol::table& def);
    ~Symplectic();

    void Theta_B(MFP* mfp, Real time, Real dt) const;
    void Theta_E(MFP* mfp, Real time, Real dt) const;
    void Theta_R(MFP* mfp, Real time, Real dt, int comp) const;

    virtual void apply_change(MFP* mfp, const Real time, const Real dt) override;

    virtual std::string get_tag() const override { return tag; }
    virtual ActionType get_type() const override { return ActionType::Symplectic; }

    FieldState* field;
    Vector<ChargedParticle*> species;

    int time_order;
    int allow_field_generation;  // allow generation of fields from the motion of charged particles

    static bool registered;
    static std::string tag;
};

#endif
#endif  // SYMPLECTIC_H
