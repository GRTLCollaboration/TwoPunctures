#ifndef __TP_PARAMETERS__
#define __TP_PARAMETERS__

#include <string>

// The TP Parameter set from the Cactus
// as header. Yeah, it's not nice, but for the time being okay
// for a first compilation

namespace TP
{
struct Parameters
{
    Parameters();

    // Cactus inherits are still problematic:
    // ie. there's a conformal_state variable probably inherited from thorn
    // StaticConformal

    bool verbose;
    bool keep_u_around;
    bool give_bare_mass;
    double adm_tol;
    std::string grid_setup_method;
    std::string initial_lapse;
    int npoints_A;
    int npoints_B;
    int npoints_phi;
    double Newton_tol;
    int Newton_maxit;
    double TP_epsilon;
    double TP_Tiny;
    double TP_Extend_Radius;
    double par_b;
    double par_bv;
    double par_m_plus;
    double par_m_minus;
    double target_M_plus;
    double target_M_minus;
    double par_P_plus[3];
    double par_P_minus[3];
    double par_S_plus[3];
    double par_S_minus[3];
    double center_offset[3];
    double initial_lapse_psi_exponent;
    bool swap_xz;
    bool use_sources;
    bool rescale_sources;
    bool use_external_initial_guess;
    bool do_residuum_debug_output;
    bool do_initial_debug_output;
    bool multiply_old_lapse;
    bool schedule_in_ADMBase_InitialData;
    bool solve_momentum_constraint;

    // From thorn Einsteinbase/StaticConformal:
    std::string metric_type;
    std::string conformal_storage;
    int conformal_state;

    // Interface of TwoPunctures:
    double J1, J2, J3;
    double mp, mm, mp_adm, mm_adm;

}; // class Parameters
} // namespace TP

#endif /* __TP_PARAMETERS__ */
