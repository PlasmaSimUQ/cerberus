#include "MFP_source.H"

#include "MFP_global.H"
#include "MFP_ode_solver.h"
#include "MFP_ode_system.h"
#include "MFP_lua.H"
#include "MFP_state.H"
#include "Core"

using GD = GlobalData;
using dual = autodiff::dual;

//---------------------------------------------------------------------------------------------

SourceTerm::SourceTerm(){}

SourceTerm::SourceTerm(Vector<int> apply, const std::string label)
{
    name = label;

    offsets.resize(apply.size());
    for (int idx=0; idx<apply.size(); ++idx) {
        offsets[idx].local = idx;
        offsets[idx].global = apply[idx];
    }
}

SourceTerm::~SourceTerm(){}

void SourceTerm::get_includes(const sol::table& def, state_valid valid, Vector<std::string>& includes, Vector<int>& index)
{

    // get the type
    std::string source_type = def["type"].get_or<std::string>("");

    if (source_type.empty())
        Abort("Source '"+name+"' failed as source type is not defined. See 'default.lua' for valid options.");

    // do we have an 'include' list?
    get_numbered(def, includes);
    if (includes.empty())
        Abort("No input states defined for source '"+name+"'");

    // make sure our includes are available
    std::pair<bool, int> find;
    for (const auto& state_name : includes) {
        find = findInVector(GD::state_names, state_name);
        if (!find.first) {
            Abort("Source '"+name+"' failed as state '"+state_name+"' is not available. Available states are:"+vec2str(GD::state_names));
        }
    }

    // iterate over all states and grab the index of those that are listed and applicable
    // or, if the list is empty, all those that are applicable

    // if the list of states is empty then assume we will apply to all


    for (const auto& name : includes) {
        if ( GD::state_index.find(name) == GD::state_index.end() ) {
            Abort("Attempting to reference a state that doesn't exist");
        }

        const int idx = GD::state_index[name];

        // check if this state matches the source term
        if (valid) {
            if (!(*valid)(idx)) {
                Abort("Input '"+name+"' is not a valid state for source '"+name+"'");
            }
        }

        index.push_back(GD::state_index[name]);


    }

    if (index.empty())
        Abort("Source '"+name+"' failed as it is not applied to any states");


    // update the states so they know what they are attached to
    for (const int& i1 : index) {
        State& state1 = GD::get_state(i1);
        for (const int& i2 : index) {
            if (i1 == i2) continue;
            State& state2 = GD::get_state(i2);

            Vector<int>& alist = state1.associated_state[state2.get_association_type()];

            if (std::find(alist.begin(), alist.end(), i2) == alist.end()) {
              // i2 not in 'alist', add it
              alist.push_back(i2);
            }
        }
    }

    return;
}

void SourceTerm::set_reconstruction(const sol::table& def)
{
    PhysicsFactory<Reconstruction> rfact = GetReconstructionFactory();

    sol::optional<std::string> rec = def["reconstruction"];

    if (rec) {
        reconstruction = rfact.Build(rec.value(), def);

        if (!reconstruction)
            Abort("Invalid reconstruction option '"+rec.value()+"'. Options are "+vec2str(rfact.getKeys()));
    }
}

void SourceTerm::set_local_idx_validity(const int local_idx, const bool valid)
{
    int global_idx = parent->local2global_index[local_idx];
    set_global_idx_validity(global_idx, valid);
}

void SourceTerm::set_global_idx_validity(const int global_idx, const bool valid)
{
    // need to search through list for item corresponding to the passed in global index
    for (auto& offset : offsets) {
        if (offset.global == global_idx) {
            offset.valid = valid;
        }
    }
}

int SourceTerm::fun_rhs(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &ydot, Real dt) const
{
    return 0;
}

int SourceTerm::fun_jac(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &J) const
{
    return 0;
}

int SourceTerm::num_jac(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &ydot, Vector<Real> &J) const
{
    int n_terms = y0.size();

    Vector<Real> f_p(n_terms);
    Vector<Real> y_p(y0.begin(), y0.end());
    ydot.resize(n_terms);
    J.resize(n_terms*n_terms);

    fun_rhs(x, y, z, t, y0, ydot);

    const Real eps = 1e-12;

    for (int i=0; i<n_terms; ++i) {
        y_p[i] += eps;
        fun_rhs(x, y, z, t, y_p, f_p);
        for (int j=0; j<n_terms; ++j) {
            J[j + i*n_terms] = (f_p[j] - ydot[j])/eps;
        }
        y_p[i] = y0[i];
    }

    return 0;

}

bool SourceTerm::get_hydro_prim(Vector<Real> &y, Vector<Real> &hydro_prim, const int offset)
{

    int y_cnt = 0; // offset within the y vector
    int h_cnt = 0; // offset within the hydro_prim vector

    Vector<Real> Q;
    bool valid = true;

    for (const auto &idx : offsets) {
        State &istate = *GD::states[idx.global];
        const int t = istate.get_type();
        const int num = istate.n_cons();
        if ((t == +StateType::isHydro) || (t == +StateType::isHydro2P)) {
            Q.resize(num);

            // get a copy of the species state
            std::copy(&y.at(y_cnt + offset), &y.at(y_cnt + offset) + num, Q.begin());

            // convert to primitive
            valid &= istate.cons2prim(Q);

            // load it into the overall vector
            std::copy(Q.begin(), Q.end(), &hydro_prim.at(h_cnt));

            h_cnt += num;
        }
        y_cnt += num;
    }

    return valid;

}

Real SourceTerm::get_max_freq(Vector<Real> &y) const
{
    return std::numeric_limits<Real>::min(); // a very small number
}

int SourceTerm::num_slopes() const
{
    return 0;
}

PhysicsFactory<SourceTerm>& GetSourceTermFactory()
{
    static PhysicsFactory<SourceTerm> F;
    return F;
}
