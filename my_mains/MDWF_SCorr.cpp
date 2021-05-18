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

    // Set additional input parameters
    // default values
    int Ls_in=8;
    double M5_in=1.8;
    double b5_in=1.5;
    double c5_in=0.5;
    double ml_in=0.001;
    double ms_in=0.01;
    std::string conf_name_in="test";
    double conf_number_in =  1000;
    for(int i=0;i<argc;i++){
          if(std::string(argv[i]) == "-Ls"){
            std::stringstream ss(argv[i+1]); ss >> Ls_in;
            }
          if(std::string(argv[i]) == "-M5"){
            std::stringstream ss(argv[i+1]); ss >> M5_in;
            }
          if(std::string(argv[i]) == "-b5"){
            std::stringstream ss(argv[i+1]); ss >> b5_in;
            }
          if(std::string(argv[i]) == "-c5"){
            std::stringstream ss(argv[i+1]); ss >> c5_in;
            }
          if(std::string(argv[i]) == "-conf_name"){
            std::stringstream ss(argv[i+1]); ss >> conf_name_in;
            }
          if(std::string(argv[i]) == "-conf_number"){
            std::stringstream ss(argv[i+1]); ss >> conf_number_in;
            }
          if(std::string(argv[i]) == "-ml"){
            std::stringstream ss(argv[i+1]); ss >> ml_in;
            }
          if(std::string(argv[i]) == "-ms"){
            std::stringstream ss(argv[i+1]); ss >> ms_in;
            }
        }
    // Show additional input parameters
    LOG(Message) << "Ls = " << Ls_in << std::endl;
    LOG(Message) << "M5 = " << M5_in << std::endl;
    LOG(Message) << "b5 = " << b5_in << std::endl;
    LOG(Message) << "c5 = " << c5_in << std::endl;
    LOG(Message) << "ml = " << ml_in << std::endl;
    LOG(Message) << "ms = " << ms_in << std::endl;
    LOG(Message) << "conf_name = " << conf_name_in << std::endl;
    LOG(Message) << "conf_number = " << conf_number_in << std::endl;



    // run setup ///////////////////////////////////////////////////////////////
    Application              application;
    std::vector<std::string> flavour = {"l", "s" /* ,"l", "s", "c1", "c2", "c3"*/};
    std::vector<double>      mass    = {ml_in, ms_in /* ,.01, .04, .2  , .25 , .3*/  };

    // global parameters
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start = conf_number_in;
    globalPar.trajCounter.end   = conf_number_in + 1;
    globalPar.trajCounter.step  = 1;
    globalPar.runId             = conf_name_in;
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
        actionPar.Ls    = Ls_in ;
        actionPar.M5    = M5_in ;
        actionPar.b     = b5_in ;
        actionPar.c     = c5_in ;
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

        smesPar.output  = "smesons/pt_" + flavour[i] + flavour[j] + "_" + conf_name_in ;
        smesPar.q1      = "Qpt_" + flavour[i];
        smesPar.q2      = "Qpt_" + flavour[j];
        smesPar.gammas  = "all";
        smesPar.sink    = "ssink";
        application.createModule<MContraction::SMeson>("smeson_pt_"
                                                      + flavour[i] + flavour[j],
                                                      smesPar);

    }

    // execution
    application.saveParameterFile("smeson_spectrum.xml");
    application.run();

    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();

    return EXIT_SUCCESS;
}
