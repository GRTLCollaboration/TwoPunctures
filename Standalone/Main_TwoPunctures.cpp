/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "CH_Timer.H"
#include "parstream.H" //Gives us pout()
#include <chrono>
#include <iostream>

#include "GRParmParse.hpp"
#include "SimulationParameters.hpp"

// Problem specific includes:
#include "SmallDataIO.hpp"
#include "TwoPunctures.hpp"

int runTwoPunctures(int argc, char *argv[])
{
    // Load the parameter file and construct the SimulationParameter class
    // To add more parameters edit the SimulationParameters file.
    char *in_file = argv[1];
    GRParmParse pp(argc - 2, argv + 2, NULL, in_file);
    SimulationParameters sim_params(pp);

    using namespace TP;
    TwoPunctures two_punctures;

    two_punctures.Parameters::operator=(sim_params.tp_params);
    two_punctures.Run();

    // For debugging
    /*
    {
        SmallDataIO tp_x_data_file("two_punctures_x", 1.0, 0.0, 0.0,
                                   SmallDataIO::APPEND, true);
        tp_x_data_file.write_header_line({"g11", "K12", "lapse"}, "x");
        const double finest_dx =
            sim_params.coarsest_dx / pow(2.0, sim_params.max_level);
        const int num_points = sim_params.L / finest_dx + 1;
        for (int ix = 0; ix < num_points; ++ix)
        {
            std::array<double, CH_SPACEDIM> coords;
            coords[0] = -sim_params.L / 2.0 + ix * finest_dx;
            coords[1] = 0.0;
            coords[2] = 0.0;

            using namespace TP::Z4VectorShortcuts;
            double TP_state[Qlen];
            tp_amr.m_two_punctures.Interpolate(coords.data(), TP_state);
            tp_x_data_file.write_data_line(
                {TP_state[g11], TP_state[K12], TP_state[lapse]}, coords[0]);
        }
    }
    */

    CH_TIMER_REPORT(); // Report results when running with Chombo timers.

    return 0;
}

int main(int argc, char *argv[])
{
    int status = runTwoPunctures(argc, argv);

    if (status == 0)
        pout() << "TwoPunctures finished." << std::endl;
    else
        pout() << "TwoPunctures failed with return code " << status
               << std::endl;

    return status;
}
