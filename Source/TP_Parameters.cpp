#include "TwoPunctures.hpp"

namespace TP
{
Parameters::Parameters()
{
    // set all the parameters
    verbose = false;
    keep_u_around = false;
    give_bare_mass = false;
    adm_tol = 1e-10;
    grid_setup_method = "Taylor expansion";
    initial_lapse = "twopunctures-averaged";
    npoints_A = 30;
    npoints_B = 30;
    npoints_phi = 16; // has to be multiples of 4
    Newton_tol = 1e-10;
    Newton_maxit = 5;
    TP_epsilon = 0;
    TP_Tiny = 0;
    TP_Extend_Radius = 0;
    par_b = 1; // Position
    par_bv = 1;
    par_m_plus = 1;
    par_m_minus = 1;
    // target ADM masses
    target_M_plus = 0.5;
    target_M_minus = 0.5;
    for (int i = 0; i < 3; i++)
    {
        // Linear momenta
        par_P_plus[i] = 0.5;
        par_P_minus[i] = 0.5;
        // Spins
        par_S_plus[i] = 0;
        par_S_minus[i] = 0;
        // Positions
        center_offset[i] = 0;
    }
    initial_lapse_psi_exponent = -2.0;
    swap_xz = false;
    use_sources = false;
    rescale_sources = true;
    use_external_initial_guess = false;
    do_residuum_debug_output = false;
    do_initial_debug_output = false;
    multiply_old_lapse = false;
    schedule_in_ADMBase_InitialData = true;
    solve_momentum_constraint = false;

    // conformal_state must be zero!
    metric_type = "something else";
    conformal_storage = "not conformal at all";

    // These are presumably pointers due limitations of cactus interface
    // description
    conformal_state = 0;
    mp = 0;
    mm = 0;
    mp_adm = 0;
    mm_adm = 0;

    // Move this print output to SimulationParameters
    /*
    TP_INFO("The two black holes have target ADM masses of %f and %f",
            target_M_minus, target_M_plus);
    TP_INFO("The corresponding linear momenta are:");
    TP_INFO("P_x1:%f\tP_x2:%f", par_P_minus[0], par_P_plus[0]);
    TP_INFO("P_y1:%f\tP_y2:%f", par_P_minus[1], par_P_plus[1]);
    TP_INFO("P_z1:%f\tP_z2:%f", par_P_minus[2], par_P_plus[2]);
    TP_INFO("The corresponding spins are:");
    TP_INFO("S_x1:%f\tS_x2:%f", par_S_minus[0], par_S_plus[0]);
    TP_INFO("S_y1:%f\tS_y2:%f", par_S_minus[1], par_S_plus[1]);
    TP_INFO("S_z1:%f\tS_z2:%f", par_S_minus[2], par_S_plus[2]);
    */
} // constructor
} // namespace TP
