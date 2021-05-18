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
    std::vector<double>      mass    = {0.1/* ,.01, .04, .2  , .25 , .3*/  };

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
    MSink::Point::Par sinkPar;
    sinkPar.mom = "0 0 0";
    application.createModule<MSink::ScalarPoint>("sink", sinkPar);

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
        MContraction::Meson::Par mesPar;

        mesPar.output  = "tmesons/pt_" + flavour[i] + flavour[j];
        mesPar.q1      = "Qpt_" + flavour[i];
        mesPar.q2      = "Qpt_" + flavour[j];
        mesPar.gammas  = "all";
        mesPar.sink    = "sink";
        application.createModule<MContraction::Meson>("tmeson_pt_"
                                                      + flavour[i] + flavour[j],
                                                      mesPar);

    }

    // execution
    application.saveParameterFile("t_meson_spectrum.xml");
    application.run();

    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();

    return EXIT_SUCCESS;
}
