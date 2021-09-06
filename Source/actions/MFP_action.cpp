#include "MFP_action.H"
#include "MFP.H"
#include "MFP_eulerian.H"

Action::Action(){}
Action::~Action(){}

ClassFactory<Action>& GetActionFactory() {
    static ClassFactory<Action> F;
    return F;
}
