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
    unsigned int             Ns      = GridDefaultLatt()[Xp];
    unsigned int             Nt      = GridDefaultLatt()[Tp];
    // Set additional input parameters
    // default values
    std::vector<std::string> flavour;
    std::vector<std::string> mass_str;
    std::vector<double>      mass;
    int Ls_in=8;
    double M5_in=1.8;
    double b5_in=1.50;
    double c5_in=0.50;
    // container to read in light, charm and strange quark masses
    double ml_in;
    double ms_in;
    double mc_in;

    double residual_in=1.0e-16;
    double max_it_in=100000;
    std::string path_conf_in = "" ;
    std::string pre_folder_in = "" ;
    bool pre_folder_set = false ;
    std::string conf_name_in="unit";
    double conf_number_in =  1000;
    bool tdir_in = false ;
    std::string tdir_mom_low = "0.0.0" ;
    std::string tdir_mom_up = "0.0.0" ;
    std::vector<unsigned int> tmom_low;
    std::vector<unsigned int> tmom_up;

    bool sdir_in = false ;
    std::string sdir_mom_low = "0.0.0" ;
    std::string sdir_mom_up = "0.0.0" ;
    std::vector<unsigned int> smom_low;
    std::vector<unsigned int> smom_up;

    bool mres_in = false;
    std::string paramstring = "" ; // String to put all settings into the outputname
    std::string tmp_str ; // to tempolary storing strings
    std::string str_M5 ; // strings which will be added to the paramstring
    std::string str_b5 ;
    std::string str_c5 ;

    if( GridCmdOptionExists(argv,argv+argc,"-tdir-lowmom") ){
        tdir_mom_low= GridCmdOptionPayload(argv,argv+argc,"-tdir-lowmom");
    }
    if( GridCmdOptionExists(argv,argv+argc,"-tdir-upmom") ){
        tdir_mom_up= GridCmdOptionPayload(argv,argv+argc,"-tdir-upmom");
    }
    if( GridCmdOptionExists(argv,argv+argc,"-sdir-lowmom") ){
        sdir_mom_low= GridCmdOptionPayload(argv,argv+argc,"-sdir-lowmom");
    }
    if( GridCmdOptionExists(argv,argv+argc,"-sdir-upmom") ){
        sdir_mom_up= GridCmdOptionPayload(argv,argv+argc,"-sdir-upmom");
    }

    LOG(Message) << "Before the for loop" << std::endl;
    for(int i=0;i<argc;i++){
          if(std::string(argv[i]) == "-Ls"){
            std::stringstream ss(argv[i+1]); ss >> Ls_in;
            //paramstring += "Ls" + std::to_string(Ls_in);
            }
          //if(std::string(argv[i]) == "-tdir-lowmom"){
          //  std::stringstream ss(argv[i+1]); tdir_mom_low = ss.str();
          //  }
          //if(std::string(argv[i]) == "-tdir-upmom"){
          //  std::stringstream ss(argv[i+1]); tdir_mom_up = ss.str();
          //  }

          //if(std::string(argv[i]) == "-sdir-lowmom"){
          //  std::stringstream ss(argv[i+1]); sdir_mom_low = ss.str();
          //  }
          //if(std::string(argv[i]) == "-sdir-upmom"){
          //  std::stringstream ss(argv[i+1]); sdir_mom_up = ss.str();

          //  }
          if(std::string(argv[i]) == "-M5"){
            std::stringstream ss(argv[i+1]); ss >> M5_in;
	    //tmp_str = std::to_string(M5_in);
           // str_M5 = "M" + tmp_str.substr(0,1) + tmp_str.substr(2,1);
            }
          if(std::string(argv[i]) == "-b5"){
            std::stringstream ss(argv[i+1]); ss >> b5_in;
	    //tmp_str = std::to_string(b5_in);
            //str_b5 = "b" + tmp_str.substr(0,1) + tmp_str.substr(2,2);
            }
          if(std::string(argv[i]) == "-c5"){
            std::stringstream ss(argv[i+1]); ss >> c5_in;
	    //tmp_str = std::to_string(c5_in);
            //str_c5 = "c" + tmp_str.substr(0,1) + tmp_str.substr(2,2);

	  }
          if(std::string(argv[i]) == "-conf_name"){
            std::stringstream ss(argv[i+1]); ss >> conf_name_in;
            }
          if(std::string(argv[i]) == "-pre_folder"){
            std::stringstream ss(argv[i+1]); ss >> pre_folder_in;
	    pre_folder_set=true;
            }
          if(std::string(argv[i]) == "-conf_number"){
            std::stringstream ss(argv[i+1]); ss >> conf_number_in;
            }
          if(std::string(argv[i]) == "-residual"){
            std::stringstream ss(argv[i+1]); ss >> residual_in;
            }
          if(std::string(argv[i]) == "-max-it"){
            std::stringstream ss(argv[i+1]); ss >> max_it_in;
            }
          if(std::string(argv[i]) == "-ml"){
            std::stringstream ss(argv[i+1]); ss >> ml_in;
            flavour.push_back("l") ;
            mass.push_back(ml_in) ;
	    mass_str.push_back(ss.str());
            //std::string tmp_ml = ss.str();
            //std::string tmp_ml = std::to_string(ml_in) ;
	    //LOG(Message) << "tmp_ml = " << tmp_ml << std::endl ;
            //LOG(Message) << "tmp_ml_sub = " << tmp_ml.substr(2) << std::endl ;
	    //paramstring += "ml" + tmp_ml.substr(2);
            }
          if(std::string(argv[i]) == "-ms"){
            std::stringstream ss(argv[i+1]); ss >> ms_in;
            flavour.push_back("s") ;
            mass.push_back(ms_in) ;
	    mass_str.push_back(ss.str());
	    //std::string tmp_ms = std::to_string(ms_in) ;
            //std::string tmp_ms = ss.str();
	    //paramstring += "ms" + tmp_ms.substr(2);
            }
          if(std::string(argv[i]) == "-mc"){
            std::stringstream ss(argv[i+1]); ss >> mc_in;
            flavour.push_back("c") ;
            mass.push_back(mc_in) ;
	    mass_str.push_back(ss.str());
	    //std::string tmp_ms = std::to_string(ms_in) ;
            //std::string tmp_ms = ss.str();
	    //paramstring += "ms" + tmp_ms.substr(2);
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
    if (!pre_folder_set){
        pre_folder_in = conf_name_in ;
    }

   LOG(Message) << "After the for loop" << std::endl;
   // set the one part of the outputname

   str_M5 = std::to_string(M5_in);
   str_M5 = "M" + str_M5.substr(0,1) + str_M5.substr(2,1);
   str_b5 = std::to_string(b5_in);
   str_b5 = "b" + str_b5.substr(0,1) + str_b5.substr(2,2);
   str_c5 = std::to_string(c5_in);
   str_c5 = "c" + str_c5.substr(0,1) + str_c5.substr(2,2);

   paramstring= "Ls" + std::to_string(Ls_in) + str_b5 + str_c5 + str_M5 ;

   LOG(Message) << "After paramstring is set" << std::endl;
   // if (!tdir_in && !sdir_in && !mres_in){
   //   LOG(Message) << "You have to use at least -sdir or -tdir \
   //   to set a direction for the contraction or -mres to calculate mres." << std::endl;
   //   exit(0);
   // }
    // Show additional input parameters
    LOG(Message) << "Ns = " << Ns << std::endl;
    LOG(Message) << "Nt = " << Nt << std::endl;
    LOG(Message) << "Ls = " << Ls_in << std::endl;
    LOG(Message) << "M5 = " << M5_in << std::endl;
    LOG(Message) << "b5 = " << b5_in << std::endl;
    LOG(Message) << "c5 = " << c5_in << std::endl;
    LOG(Message) << "paramstring = " << paramstring << std::endl;
    for (int i = 0; i< flavour.size() ; i++){
      LOG(Message) << flavour[i] << " = " << mass[i] << std::endl;
    }
    LOG(Message) << "conf_path = " << path_conf_in << std::endl;
    LOG(Message) << "conf_name = " << conf_name_in << std::endl;
    LOG(Message) << "conf_number = " << conf_number_in << std::endl;
    LOG(Message) << "pre_folder = " << pre_folder_in << std::endl;
    LOG(Message) << "calculate temporal correlator = " << tdir_in << std::endl;
    LOG(Message) << "momentum fof temporal correlator goes from " << tdir_mom_low << " up to " << tdir_mom_up << std::endl;
    LOG(Message) << "calculate spatial correlator = " << sdir_in << std::endl;
    LOG(Message) << "momentum of spatial correlator goes from " << sdir_mom_low << " up to " << sdir_mom_up << std::endl;
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
       conf_name_in = conf_name_in + std::to_string(Ns) + std::to_string(Nt);
    }

    // sources
    MSource::Point::Par ptPar;
    ptPar.position = "0 0 0 0";
    application.createModule<MSource::Point>("pt", ptPar);


    // set fermion boundary conditions to be periodic space, antiperiodic time.
    std::string boundary = "1 1 1 -1";
    std::string twist = "0. 0. 0. 0.";

    // convert momenta strings to vectors
    //tmom_low = strToVec<unsigned int>(tdir_mom_low);
    //tmom_up  = strToVec<unsigned int>(tdir_mom_up);
    //smom_low = strToVec<unsigned int>(sdir_mom_low);
    //smom_up  = strToVec<unsigned int>(sdir_mom_up);
    GridCmdOptionIntVector(tdir_mom_low,tmom_low);
    GridCmdOptionIntVector(tdir_mom_up,tmom_up);
    GridCmdOptionIntVector(sdir_mom_low,smom_low);
    GridCmdOptionIntVector(sdir_mom_up,smom_up);

    //check if vectors have the right size
    if( tmom_low.size() != 3 ){
        std::cout << "The size of --tdir-lowmom has to be 3 and not " << tmom_low.size() << std::endl;
        exit(0);
    }
    if( tmom_up.size() != 3 ){
        std::cout << "The size of --tdir-upmom has to be 3 and not " << tmom_up.size() << std::endl;
        exit(0);
    }
    if( smom_low.size() != 3 ){
        std::cout << "The size of --sdir-lowmom has to be 3 and not " << smom_low.size() << std::endl;
        exit(0);
    }
    if( smom_up.size() != 3 ){
        std::cout << "The size of --sdir-upmom has to be 3 and not " << smom_up.size() << std::endl;
        exit(0);
    }




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
        solverPar.residual     = residual_in;
        solverPar.maxIteration = max_it_in;
        application.createModule<MSolver::RBPrecCG>("CG_" + flavour[i],
                                                    solverPar);

        // propagators
        MFermion::GaugeProp::Par quarkPar;
        quarkPar.solver = "CG_" + flavour[i];
        quarkPar.source = "pt";
        application.createModule<MFermion::GaugeProp>("Qpt_" + flavour[i], quarkPar);
    }

    bool first_time = true;
    for (unsigned int i = 0; i < flavour.size(); ++i)
    {
        std::string tparamstring_fin = "";
        std::string sparamstring_fin = "";
        for (unsigned int j = i; j < flavour.size(); ++j)
        {
	    if( flavour[i] != flavour[j] ){
	        tparamstring_fin= "m" + flavour[i] + mass_str[i].substr(2) + "m" + flavour[j] + mass_str[j].substr(2) + paramstring;
	        sparamstring_fin= "m" + flavour[i] + mass_str[i].substr(2) + "m" + flavour[j] + mass_str[j].substr(2) + paramstring;
	    }
	    else{
	        tparamstring_fin= "m" + flavour[i] + mass_str[i].substr(2) + paramstring ;
	        sparamstring_fin= "m" + flavour[i] + mass_str[i].substr(2) + paramstring ;
	    }


	   if (tdir_in){
            //LOG(Message) << tmom_low[0] <<  tmom_low[1] << tmom_low[2] << std::endl;
            //LOG(Message) << tmom_up[0] <<  tmom_up[1] << tmom_up[2] << std::endl;
       for (unsigned int kt_x = tmom_low[0];kt_x < tmom_up[0] +1 ; kt_x++)
       {
       for (unsigned int kt_y = tmom_low[1];kt_y < tmom_up[1] +1 ; kt_y++)
       {
       for (unsigned int kt_z = tmom_low[2];kt_z < tmom_up[2] +1 ; kt_z++)
       {
              std::string tmom_string = "";
              std::string tdir_mom_in = std::to_string(kt_x) + " " + std::to_string(kt_y) + " " + std::to_string(kt_z);

              if(tdir_mom_in != "0 0 0"){
                  std::string tmp_tmom = tdir_mom_in ;
                  tmp_tmom.erase(std::remove(tmp_tmom.begin(),tmp_tmom.end(), ' '),tmp_tmom.end());
                  tmom_string = "kt" + tmp_tmom;
              }

              //if( tmom_string != ""){
              //    tparamstring_fin = tparamstring_fin + tmom_string ;
              //}

              // sink for the temporal direction
              if(first_time){
                LOG(Message) << "Create tsink with kt = " << tmom_string << std::endl;
                MSink::Point::Par tsinkPar;
                tsinkPar.mom = tdir_mom_in;
                application.createModule<MSink::ScalarPoint>("tsink" + tmom_string, tsinkPar);
              }
              //Contraction in the temporal direction
              MContraction::Meson::Par tmesPar;

              tmesPar.output  = "../data/tmesons/" + pre_folder_in  + "/pt_" + flavour[i] + flavour[j] + "_t_"
                                 + tparamstring_fin + tmom_string + "_" + conf_name_in ;
              tmesPar.q1      = "Qpt_" + flavour[i];
              tmesPar.q2      = "Qpt_" + flavour[j];
              tmesPar.gammas  = "all";
              tmesPar.sink    = "tsink" + tmom_string;
              application.createModule<MContraction::Meson>("tmeson_pt_"
                                                            + flavour[i] + flavour[j] + tmom_string,
                                                            tmesPar);

      }
      }
      }
    }

            if (sdir_in){
            //LOG(Message) << smom_low[0] <<  smom_low[1] << smom_low[2] << std::endl;
            //LOG(Message) << smom_up[0] <<  smom_up[1] << smom_up[2] << std::endl;
              for (unsigned int ks_x = smom_low[0];ks_x < smom_up[0] +1 ; ks_x++)
              {
              for (unsigned int ks_y = smom_low[1];ks_y < smom_up[1] +1 ; ks_y++)
              {
              for (unsigned int ks_t = smom_low[2];ks_t < smom_up[2] +1 ; ks_t++)
              {
                std::string smom_string = "";
                std::string sdir_mom_in = std::to_string(ks_x) + " " + std::to_string(ks_y) + " " + std::to_string(ks_t);

                if(sdir_mom_in != "0 0 0"){
                   std::string tmp_smom = sdir_mom_in ;
                   tmp_smom.erase(std::remove(tmp_smom.begin(),tmp_smom.end(), ' '),tmp_smom.end());
                   smom_string = "ks" + tmp_smom;
                }

                //if( smom_string != ""){
                //   sparamstring_fin = sparamstring_fin + smom_string ;
                //}

                // sink for the spatial direction
                if(first_time){
                  LOG(Message) << "Create ssink with ks = " << smom_string << std::endl;
                  MSink::SPoint::Par ssinkPar;
                  ssinkPar.mom = sdir_mom_in;
                  application.createModule<MSink::ScalarSPoint>("ssink" + smom_string, ssinkPar);
                }
                //Contraction in the spatial direction
                MContraction::SMeson::Par smesPar;


                smesPar.output  = "../data/smesons/" + pre_folder_in  + "/pt_" + flavour[i] + flavour[j] + "_s_"
                                  + sparamstring_fin + smom_string  + "_" + conf_name_in ;
                smesPar.q1      = "Qpt_" + flavour[i];
                smesPar.q2      = "Qpt_" + flavour[j];
                smesPar.gammas  = "all";
                smesPar.sink    = "ssink" + smom_string;
                application.createModule<MContraction::SMeson>("smeson_pt_"
                                                              + flavour[i] + flavour[j] + smom_string,
                                                              smesPar);
               }
               }
               }
             }
        first_time = false;
        }

	tparamstring_fin= "m" + flavour[i] + mass_str[i].substr(2) + paramstring ;
        // Calculate mres
        if (mres_in){
            MContraction::WardIdentity::Par wardIdpar;
            wardIdpar.output = "../data/wardidentity/" + pre_folder_in  +"/ward_pt_" + flavour[i] + "_"
                               + tparamstring_fin + "_" + conf_name_in ;
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
