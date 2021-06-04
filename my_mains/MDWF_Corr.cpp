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
    std::vector<std::string> flavour;
    std::vector<double>      mass;
    int Ls_in=8;
    double M5_in=1.8;
    double b5_in=1.5;
    double c5_in=0.5;
    double ml_in=0.001;
    double ms_in=0.01;
    std::string path_conf_in = "" ;
    std::string conf_name_in="unit";
    double conf_number_in =  1000;
    bool tdir_in = false ;
    bool sdir_in = false ;
    bool mres_in = false;
    std::string paramstring = "" ; // String to put all settings into the outputname
    for(int i=0;i<argc;i++){
          if(std::string(argv[i]) == "-Ls"){
            std::stringstream ss(argv[i+1]); ss >> Ls_in;
            paramstring += "Ls" + to_string(Ls_in);
            }
          if(std::string(argv[i]) == "-M5"){
            std::stringstream ss(argv[i+1]); ss >> M5_in;
            paramstring += "M5" + to_string(M5_in * 10);
            }
          if(std::string(argv[i]) == "-b5"){
            std::stringstream ss(argv[i+1]); ss >> b5_in;
            paramstring += "b5" + to_string(b5_in * 10);
            }
          if(std::string(argv[i]) == "-c5"){
            std::stringstream ss(argv[i+1]); ss >> c5_in;
            paramstring += "c5" + to_string(c5_in * 10);
            }
          if(std::string(argv[i]) == "-conf_name"){
            std::stringstream ss(argv[i+1]); ss >> conf_name_in;
            }
          if(std::string(argv[i]) == "-conf_number"){
            std::stringstream ss(argv[i+1]); ss >> conf_number_in;
            }
          if(std::string(argv[i]) == "-ml"){
            std::stringstream ss(argv[i+1]); ss >> ml_in;
            flavour.push_back("l") ;
            mass.push_back(ml_in) ;
            std:string tmp_ml = to_string(ml_in) ;
            paramstring += "ml" + tmp_ml.substr(2);
            }
          if(std::string(argv[i]) == "-ms"){
            std::stringstream ss(argv[i+1]); ss >> ms_in;
            flavour.push_back("s") ;
            mass.push_back(ms_in) ;
            std:string tmp_ms = to_string(ms_in) ;
            paramstring += "ms" + tmp_ms.substr(2);
            }
          if(std::string(argv[i]) == "-path_conf"){
            std::stringstream ss(argv[i+1]); ss >> path_conf_in;
            }
          if(std::string(argv[i]) == "-tdir"){
              tdir_in = true;
              }
          if(std::string(argv[i]) == "-sdir"){
                  sdir_in = true;
                  }
          if(std::string(argv[i]) == "-mres"){
                  mres_in = true;
                  }
        }
    if (!tdir_in && !sdir_in && !mres_in){
      LOG(Message) << "You have to use at least -sdir or -tdir \
      to set a direction for the contraction or -mres to calculate mres." << std::endl;
      exit(0);
    }
    // Show additional input parameters
    LOG(Message) << "Ls = " << Ls_in << std::endl;
    LOG(Message) << "M5 = " << M5_in << std::endl;
    LOG(Message) << "b5 = " << b5_in << std::endl;
    LOG(Message) << "c5 = " << c5_in << std::endl;
    LOG(Message) << "paramstring = " << paramstring << std::endl;
    //LOG(Message) << "ml = " << ml_in << std::endl;
    //LOG(Message) << "ms = " << ms_in << std::endl;
    for (int i = 0; i< flavour.size() ; i++){
      LOG(Message) << flavour[i] << " = " << mass[i] << std::endl;
    }
    LOG(Message) << "conf_path = " << path_conf_in << std::endl;
    LOG(Message) << "conf_name = " << conf_name_in << std::endl;
    LOG(Message) << "conf_number = " << conf_number_in << std::endl;
    LOG(Message) << "calculate temporal correlator = " << tdir_in << std::endl;
    LOG(Message) << "calculate spatial correlator = " << sdir_in << std::endl;
    LOG(Message) << "calculate residual mass = " << mres_in << std::endl;



    // run setup ///////////////////////////////////////////////////////////////
    Application              application;
    //std::vector<std::string> flavour = {"l", "s" };
    //std::vector<double>      mass    = {ml_in, ms_in };

    // global parameters
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start = conf_number_in;
    globalPar.trajCounter.end   = conf_number_in + 1;
    globalPar.trajCounter.step  = 1;
    globalPar.runId             = conf_name_in;
    application.setPar(globalPar);
    // gauge field

    if(conf_name_in != "unit"){
        LOG(Message) << "Use nersc configuration" << std::endl;
	MIO::LoadNersc::Par loadnerscpar;
        // the dot and confnumber are added in the LoadNersc module
        loadnerscpar.file =  path_conf_in + "/" + conf_name_in;
        //std::cout << loadnerscpar.file << std::endl;
        application.createModule<MIO::LoadNersc>("gauge",loadnerscpar);

    }
    else
    {
       LOG(Message) << "Use unit configuration" << std::endl;
       application.createModule<MGauge::Unit>("gauge");
    }

    // sources
    MSource::Point::Par ptPar;
    ptPar.position = "0 0 0 0";
    application.createModule<MSource::Point>("pt", ptPar);
    // sink
    //MSink::SPoint::Par ssinkPar;
    //ssinkPar.mom = "0 0 0";
    //application.createModule<MSink::ScalarSPoint>("ssink", ssinkPar);


    if (tdir_in){
         // sink for the temporal direction
         MSink::Point::Par tsinkPar;
         tsinkPar.mom = "0 0 0";
         application.createModule<MSink::ScalarPoint>("tsink", tsinkPar);
    }
    if (sdir_in){
         // sink for the spatial direction
         MSink::SPoint::Par ssinkPar;
         ssinkPar.mom = "0 0 0";
         application.createModule<MSink::ScalarSPoint>("ssink", ssinkPar);
    }


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
    {

        for (unsigned int j = i; j < flavour.size(); ++j)
        {

           if (tdir_in){
              // sink for the temporal direction
              // MSink::Point::Par tsinkPar;
              // tsinkPar.mom = "0 0 0";
              //application.createModule<MSink::ScalarPoint>("tsink", tsinkPar);

              //Contraction in the spatial direction
              MContraction::Meson::Par tmesPar;

              tmesPar.output  = "tmesons//pt_" + flavour[i] + flavour[j] + "_t_"
                                 + paramstring + "_" + conf_name_in ;
              tmesPar.q1      = "Qpt_" + flavour[i];
              tmesPar.q2      = "Qpt_" + flavour[j];
              tmesPar.gammas  = "all";
              tmesPar.sink    = "tsink";
              application.createModule<MContraction::Meson>("tmeson_pt_"
                                                            + flavour[i] + flavour[j],
                                                            tmesPar);

            }

            if (sdir_in){
                // sink for the spatial direction
                //MSink::SPoint::Par ssinkPar;
                //ssinkPar.mom = "0 0 0";
                //application.createModule<MSink::ScalarSPoint>("ssink", ssinkPar);

                //Contraction in the spatial direction
                MContraction::SMeson::Par smesPar;

                smesPar.output  = "smesons/pt_" + flavour[i] + flavour[j] + "_s_"
                                  + paramstring + "_" + conf_name_in ;
                smesPar.q1      = "Qpt_" + flavour[i];
                smesPar.q2      = "Qpt_" + flavour[j];
                smesPar.gammas  = "all";
                smesPar.sink    = "ssink";
                application.createModule<MContraction::SMeson>("smeson_pt_"
                                                              + flavour[i] + flavour[j],
                                                              smesPar);

             }
        }
        // Calculate mres
        if (mres_in){
            MContraction::WardIdentity::Par wardIdpar;
            wardIdpar.output = "wardidentity/ward_pt_" + flavour[i] + "_"
                               + paramstring + "_" + conf_name_in ;
            wardIdpar.prop = "Qpt_" + flavour[i] +  "_5d";
            wardIdpar.action = "MobiusDWF_" + flavour[i];
            wardIdpar.source = "pt";
            wardIdpar.mass = mass[i];
            application.createModule<MContraction::WardIdentity>("ward_pt_"
                                                                   + flavour[i],
                                                                     wardIdpar);

       }
    }
    // execution
    application.saveParameterFile("smeson_spectrum.xml");
    application.run();

    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();

    return EXIT_SUCCESS;
}
