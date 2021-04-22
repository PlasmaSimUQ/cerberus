#include "MFP_udf.H"
#include "MFP_global.H"

using GD = GlobalData;

//---------------------------------------------------------------------------------------------

std::string UserDefinedSource::tag = "user";
bool UserDefinedSource::registered = GetSourceTermFactory().Register(UserDefinedSource::tag, SourceTermBuilder<UserDefinedSource>);

UserDefinedSource::UserDefinedSource(){}

UserDefinedSource::UserDefinedSource(const sol::table &def)
{

    name = def.get<std::string>("name");

    if (!UserDefinedSource::valid_solver(def["solver"])) {
        Abort("Error: Source '"+name+"' needs a different solver");
    }

    Vector<int> index;
    Vector<std::string> includes;

    get_includes(def, &UserDefinedSource::valid_state, includes, index);

    offsets.resize(index.size());
    for (int idx=0; idx<index.size(); ++idx) {
        offsets[idx].local = idx;
        offsets[idx].global = index[idx];
    }

    // get all of our manufactured solutions
    terms.resize(index.size());

    for (int i = 0; i<index.size(); ++i) {
        State &istate = GD::get_state(index[i]);
        Vector<std::string> cons_names = istate.get_cons_names();
        terms[i].resize(istate.n_cons());
        for (int j = 0; j<istate.n_cons(); ++j) {
            get_udf(def["value"][cons_names[j]], terms[i][j], 0.0);
        }
    }
    return;
}

UserDefinedSource::~UserDefinedSource()
{
    // do nothing
}

int UserDefinedSource::fun_rhs(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &ydot, Real dt) const
{
    BL_PROFILE("UserDefinedSource::fun_rhs");
    // cycle through the states and apply source terms

    std::map<std::string, Real> Q{{"x",x}, {"y",y}, {"z",z}, {"t",t}};

    for (const auto &idx : offsets) {

        State &istate = GD::get_state(idx.global);

        const Vector<std::string> &names = istate.get_cons_names();

        // load data into a dictionary for the lua function call
        for (int i=0; i<istate.n_cons(); ++i) {
            Q[names[i]] = y0[idx.solver + i];
        }

        // call the function
        for (int i=0; i<istate.n_cons(); ++i) {
            const Optional3D1VFunction &f = terms[idx.local][i];
            ydot[idx.solver + i] = f(Q);
        }

    }

    return 0;
}


bool UserDefinedSource::valid_state(const int global_idx)
{
    return true;
}

bool UserDefinedSource::valid_solver(const int solve_idx)
{
    if (solve_idx != +SolveType::Explicit) {
        return false;
    }
    return true;
}

std::string UserDefinedSource::print() const
{
    std::stringstream msg;

    msg << tag << " : " << name;

    for (auto& offset : offsets) {
        State &istate = GD::get_state(offset.global);
        msg << "\n        " << istate.name << " : ";
        Vector<std::string> cons_names = istate.get_cons_names();
        for (int j = 0; j<istate.n_cons(); ++j) {
            msg << "\n          " << cons_names[j] << " = " << terms[offset.local][j];
        }
    }

    return msg.str();
}
