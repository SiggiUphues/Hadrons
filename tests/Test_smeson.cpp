/*
 * Test_smeson.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Tristan Ueding <trisueding@web.de>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution
 * directory.
 */

/*  END LEGAL */

/*If you use a sink from Module MSink the SMeson (spatial correlator) only works
 for the sink SPoint since there is a slicing in z direction */

#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;

int main(int argc, char *argv[])
{
    // initialization //////////////////////////////////////////////////////////
    Grid_init(&argc, &argv);
    HadronsLogError.Active(GridLogError.isActive());
    HadronsLogWarning.Active(GridLogWarning.isActive());
    HadronsLogMessage.Active(GridLogMessage.isActive());
    HadronsLogIterative.Active(GridLogIterative.isActive());
    HadronsLogDebug.Active(GridLogDebug.isActive());
    LOG(Message) << "Grid initialized" << std::endl;

    // run setup ///////////////////////////////////////////////////////////////
    Application              application;
    std::vector<std::string> flavour = {"z"/* ,"l", "s", "c1", "c2", "c3"*/};
    std::vector<double>      mass    = {0./* ,.01, .04, .2  , .25 , .3*/  };

    // global parameters
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start = 1500;
    globalPar.trajCounter.end   = 1520;
    globalPar.trajCounter.step  = 20;
    globalPar.runId             = "test";
    application.setPar(globalPar);
    // gauge field
    application.createModule<MGauge::Unit>("gauge");
    // sources
    MSource::Point::Par ptPar;
    ptPar.position = "0 0 0 0";
    application.createModule<MSource::Point>("pt", ptPar);
    // sink
    MSink::SPoint::Par ssinkPar;
    ssinkPar.mom = "0 0 0";
    application.createModule<MSink::ScalarSPoint>("ssink", ssinkPar);

    // set fermion boundary conditions to be periodic space, antiperiodic time.
    std::string boundary = "1 1 1 -1";
    std::string twist = "0. 0. 0. 0.";

    for (unsigned int i = 0; i < flavour.size(); ++i)
    {
        // actions
        MAction::MobiusDWF::Par actionPar;
        actionPar.gauge = "gauge";
        actionPar.Ls    = 12 ;
        actionPar.M5    = 1.8 ;
        actionPar.b     = 1 ;
        actionPar.c     =0.5 ;
        actionPar.mass  = mass[i];
        actionPar.boundary = boundary;
        actionPar.twist = twist;
        application.createModule<MAction::MobiusDWF>("MobiusDWF_" + flavour[i], actionPar);

        // solvers
        MSolver::RBPrecCG::Par solverPar;
        solverPar.action       = "MobiusDWF_" + flavour[i];
        solverPar.residual     = 1.0e-16;
        solverPar.maxIteration = 10000;
        application.createModule<MSolver::RBPrecCG>("CG_" + flavour[i],
                                                    solverPar);

        // propagators
        MFermion::GaugeProp::Par quarkPar;
        quarkPar.solver = "CG_" + flavour[i];
        quarkPar.source = "pt";
        application.createModule<MFermion::GaugeProp>("Qpt_" + flavour[i], quarkPar);
    }
    for (unsigned int i = 0; i < flavour.size(); ++i)
    for (unsigned int j = i; j < flavour.size(); ++j)
    {
        MContraction::SMeson::Par smesPar;

        smesPar.output  = "smesons/pt_" + flavour[i] + flavour[j];
        smesPar.q1      = "Qpt_" + flavour[i];
        smesPar.q2      = "Qpt_" + flavour[j];
        smesPar.gammas  = "all";
        smesPar.sink    = "ssink";
        application.createModule<MContraction::SMeson>("smeson_pt_"
                                                      + flavour[i] + flavour[j],
                                                      smesPar);

    }

    // execution
    application.saveParameterFile("spectrum.xml");
    application.run();

    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();

    return EXIT_SUCCESS;
}
