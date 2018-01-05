#ifndef STRUCTS_USED_HPP
#define STRUCTS_USED_HPP

// structure to decide integrator type
enum sdeIntegratorType{
    EU_MAR,
    RK4VEC_TV_STEP
};


// structure to decide the flow field being used
enum flowFieldType{
    DG_FLOW,
    TANK_FLOW,
    DG_FLOW_CONTROL
};
#endif
