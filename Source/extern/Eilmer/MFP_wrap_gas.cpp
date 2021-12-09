#ifdef EILMER_GAS
#include "MFP_wrap_gas.H"

bool EilmerGas::initialised = false;

void EilmerGas::initialise()
{
    if (!initialised) {
        cwrap_gas_init();
        initialised = true;
    }
}

#endif
