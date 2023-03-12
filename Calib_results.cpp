// Root/C++ script to peforme count estimate calibration analysis
// 
//
//
//
//
// This function will take an config, seed, and istope, as an output a Tree
#include <sys/stat.h>
#include <unistd.h>
#include <iostream>
#include <string>
#include <glob.h>
//#include <format>



using namespace std;

class nexo_runinfo;
std::string nexo_g4_root_path = "/home/wrkshp/nEXO/slurm_data";
std::string co_search="-c ExternalGammas_light -n 5000 -i Co60";
std::string simple_sens_path="/home/wrkshp/nexo/slurm_data/sens";

TChain* GetnEXO_g4Simtree(string config="",string isotope="",int nevts=0, int seed=0 );
vector<string> parse_str(string str,string delim=" ");
TChain *Freind_all(vector<nexo_runinfo> Runlist);
void setcut();

std::string SensCut ;
std::string RecCut;
void setcut(){
  SensCut="passed_xy_thresh == 1 && passed_z_thresh== 1";
  SensCut= SensCut + " && n_x_ch_abovenoise > 0 && n_y_ch_abovenoise > 0";
  SensCut += " && (m_nOPCal <1.077*m_nQ+313)";
  SensCut += " && (m_nOPCal >0.597*m_nQ-216)";
  SensCut += " && standoff > 0 ";
  //SensCut += " && energy > 700 && energy < 3500";
  SensCut += "&& m_DNNvalue > 0.85";
  
  RecCut = "ReconFlag==2 && NumReconClusters==1";

}




class nexo_runinfo {
 
  public:
    int seed;
    int nevts;
    string config;
    string isotope;
    string filename="";
    string basename="";
    string abc="abcd";
    string option="";
    string filename_recon="";
    string basename_recon="";
    string basename_sensrecon ="";
    string filename_sensrecon="";
    string filename_simplesens="";
    string basename_simplesens="";


    void set_main (string c,string v,string iso,int n,int s) 
	{   config=c;   isotope=iso;   seed=s;   nevts=n; 
		basename= config +v, "-"+isotope;
		basename = basename + "-"+nevts+ "-seed"+seed+".root";
		filename= nexo_g4_root_path + "/" + isotope + "/" + basename;

		basename_recon = config + "-"+isotope + "-"+nevts+ "-seed"+seed+"_recon.root";
		filename_recon= nexo_g4_root_path + "/recon/" + basename_recon;
	  basename_sensrecon = config+v + "-"+isotope + "-"+nevts+ "-seed"+seed+"_sensrecon.root";
    filename_sensrecon= nexo_g4_root_path + "/recon/" + basename_sensrecon;

    basename_simplesens=isotope+v+"-seed"+seed+".root";
    filename_simplesens=simple_sens_path+"/"+basename_simplesens;



  }
  void set_all (string c,string iso,int n,int s, string o, string base, string file )
  {   config=c;   isotope=iso;   seed=s;   nevts=n; option=o;
    basename= base; filename=file;

    basename_recon = config + "-"+isotope + "-"+nevts+ "-seed"+seed+"_recon.root";
    filename_recon= nexo_g4_root_path + "/recon/" + basename_recon;
    
    basename_sensrecon = config + "-"+isotope + "-"+nevts+ "-seed"+seed+"_sensrecon.root";
    filename_sensrecon= nexo_g4_root_path + "/recon/" + basename_sensrecon;
    basename_simplesens=isotope+"-seed"+seed+".root";
    filename_simplesens=simple_sens_path+"/"+basename_simplesens;
  }

    void set_basename();

    TChain * get_g4tree();
    TChain * get_g4tree2();
};


void nexo_runinfo::set_basename()
{
	basename= config + "-"+isotope;
	basename = basename + "-"+nevts+ "-seed"+seed+".root";
	filename= nexo_g4_root_path + "/" + isotope + "/" + basename;
}

class nexo_args {
 public:
  vector<int> list_of_nevts;
  vector<string> list_of_tgts;
  vector<string> list_of_configs;  
  vector<string> list_of_strings;
  vector<int> list_of_seeds;
  string option;


};
vector<int> tovectint(vector<string> data)
{
  std::vector<int> Data;
  std::transform(data.begin(), data.end(), std::back_inserter(Data),
               [](const std::string& str) { return std::stoi(str); });
  return Data;
}

nexo_args Get_args(string argument=""){

  nexo_args Args;
  vector<string> list_of_nevts_str;
  vector<string> list_of_tgts;
  vector<string> list_of_configs;  
  vector<string> list_of_strings;
  vector<string> list_of_seeds_Str;
  string option ="";
  //parse arg
  vector<string> args_vector = parse_str(argument," ");
  

  for( std::vector<string>::iterator it=args_vector.begin(); it != args_vector.end(); ++it)
  {
    cout <<*it<<endl;
    if (*it == "-c" ||*it =="--config"){
       // Config
        ++it;
        list_of_configs.push_back(*it);
      }
    else if (*it == "-n" || *it == "-maxevts"){
        // code block
        ++it;
        list_of_nevts_str.push_back(*it);
      }
    else if (*it == "-o" || *it == "--option"){
        // code block
        ++it;
        option = *it;
      }
    else if (*it == "-i" || *it == "--isotope"){
        // code block
        ++it;
        list_of_tgts.push_back(*it);
      }
   else if (*it == "-s" || *it == "--seed"){
       // code block
       ++it;
       list_of_seeds_Str.push_back(*it);
      }
  }

  if( list_of_tgts.size()>0){ Args.list_of_tgts=list_of_tgts;}
  if( list_of_configs.size()>0){Args.list_of_configs=list_of_configs;}
  if( list_of_nevts_str.size()>0){
    vector<int> nevtsint =tovectint(list_of_nevts_str);
    Args.list_of_nevts=nevtsint;}
  if( option!=""){Args.option=option;}

  return Args;
}




inline bool exists_test3 (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}




nexo_runinfo GetRunInfoFromFile(string filename)
{

  string FullPath = filename;
  while(FullPath.find("/")!=std::string::npos)
  {
    auto pos = FullPath.find("/");
    FullPath = FullPath.substr(pos+1);
   // cout << FullPath<<endl;
  }

  string basename = FullPath;
  string Path = filename;
  Path.replace( filename.find(basename),basename.size(),"");

  string info_name= basename;
  auto cof_last= info_name.find("-");
  auto config = info_name.substr(0,cof_last);
  info_name = info_name.substr(cof_last+1);
  auto iso_last= info_name.find("-");
  auto isotope = info_name.substr(0,iso_last);
  info_name = info_name.substr(iso_last+1);
  auto n_last = info_name.find("-");
  auto nevts = info_name.substr(0,n_last);
  info_name = info_name.substr(n_last+1);
  auto seed_last=0;
  auto is_opt=0;
  string option="";

  auto seed = to_string( stoi(info_name.substr(4)) );
  info_name = info_name.substr(4+seed.size());

  if(info_name.size()>5){option=info_name.substr(0,info_name.find("."));}
  else{option="";}
 // cout <<filename<<endl;

  nexo_runinfo RI;
  RI.set_all(config,isotope,stoi(nevts),stoi(seed),option,basename,filename);

  return RI;
}


TChain* nexo_runinfo::get_g4tree(){
	return GetnEXO_g4Simtree(config,isotope,nevts,seed);
}

TChain* nexo_runinfo::get_g4tree2(){
	TChain *TC = new TChain("Event/Sim/SimEvent");

	TC->Add(filename.c_str());
	return TC;	

}	

template<typename T>
std::vector<T> arange(T start, T stop, T step = 1) {
    std::vector<T> values;
    for (T value = start; value < stop; value += step)
        values.push_back(value);
    return values;
}


vector<string> parse_str(string str,string delim=" "){
  vector <string> vec_str;
  size_t pos = 0;
  std::string token;
  while ((pos = str.find(delim)) != std::string::npos) {
    token = str.substr(0, pos);
    vec_str.push_back(token);
    str.erase(0, pos + delim.length());
  }
  if(str.size()>=1) vec_str.push_back(str);

return vec_str;
}

//OPens up the g4 event tree
TChain* GetnEXO_g4Simtree(string config="",string isotope="",int nevts=0, int seed=0 )
{

	TChain * TC = new TChain("Event/Sim/SimEvent");
	string file_name = nexo_g4_root_path + "/"+isotope+"/" + config +"-";
	file_name = file_name + isotope + "-"+nevts+ "-seed"+seed+".root";
	TFile *TF = new TFile(file_name.c_str());
	TDirectoryFile *Event = (TDirectoryFile*)TF->FindObjectAny("Event");
	TDirectoryFile *sim = (TDirectoryFile*)Event->FindObjectAny("Sim");
	TTree *T =(TTree*)sim->FindObjectAny("SimEvent");
	TBranchElement *SE = (TBranchElement*)T->GetBranch("SimEvent");
	TTree *newT = (TTree*)SE->GetTree();
	TC->Add(file_name.c_str());
	//TC->Add(file_name.c_str());
	return TC;
}

vector<nexo_runinfo> SearchForList(string arg="")
{
  vector<nexo_runinfo> Runlist;
  string config="ExternalGammas_light";
  string isotope="*";
  string nevts ="5000";
  string option="";                          
  string seed="*";
  string seedf_str="";
  string seedl_str="";
  int seedf=-1;
  int seedl=-1;
  string seed_last="";
  //parse arg
  vector<string> args = parse_str(arg," ");
  vector<string> list_of_nevts;
  vector<string> list_of_tgts;

  for( std::vector<string>::iterator it=args.begin(); it != args.end(); ++it)
  {
    cout <<*it<<endl;
    if (*it == "-c" ||*it =="--config"){
       // Config
        ++it;
        config = *it;
      }
    else if (*it == "-n" || *it == "-maxevts"){
        // code block
        ++it;
        nevts = *it;
        list_of_nevts.push_back(*it);
      }
    else if (*it == "-o" || *it == "--option"){
        // code block
        ++it;
        option = *it;
      }
    else if (*it == "-i" || *it == "--isotope"){
        // code block
        ++it;
        isotope = *it;
        list_of_tgts.push_back(*it);
      }
   else if (*it == "-sf" || *it == "--seedf"){
        // code block
        ++it;
        seedf_str = *it;
        //seeds.push_back(*it);
      }
   else if (*it == "-sl" || *it == "--seedl"){
        // code block
        ++it;
        seedl_str = *it;
        //seeds.push_back(*it);
      }

  }

  if(seedl_str!=""){seedl=stoi(seedl_str);}
  if(seedf_str!=""){seedf=stoi(seedf_str);}
  
//  if 




  vector<string> filenames;
  vector<string> num_opt= {"?","??","???","????"};
  for(auto mm=0;mm <=list_of_tgts.size()-1;mm++)
  {
    string m = list_of_tgts[mm];
    //cout<< list_of_nevts.size()<<endl;
    for( auto nn=0;nn <=list_of_nevts.size()-1;nn++)
    {
      string n = list_of_nevts[nn];
      printf("\t\t%s\n",n.c_str());
      for ( auto i =0; i <= num_opt.size()-1;i++)
      {

        string glob_key =nexo_g4_root_path+"/recon/"+config+"-";
        glob_key+=m+"-"+n+"-seed"+num_opt[i]+option+"sensrecon.root";
        printf("%s\n",glob_key.c_str());
        // glob struct resides on the stack
        glob_t glob_result;
        memset(&glob_result, 0, sizeof(glob_result));
        // do the glob operation
        int return_value = glob(glob_key.c_str(), 0, NULL, &glob_result);
        if(return_value == 0) 
        {
          // collect all the filenames into a std::list<std::string>
          for(size_t i = 0; i < glob_result.gl_pathc; ++i) 
          {

            


            filenames.push_back(string(glob_result.gl_pathv[i]));
            nexo_runinfo RL = GetRunInfoFromFile(filenames[i]);
            int checkseed = RL.seed;

            if(seedl>1)
              {if(checkseed > seedl){continue;}}
            if(seedf>1)
              {if(checkseed < seedf){continue;}}
            if(exists_test3(RL.filename_recon)&&exists_test3(RL.filename_sensrecon))
            {
              Runlist.push_back(RL);
              cout <<filenames[i] <<endl;
            }
          }//End of glob list
          filenames.clear();
          
          //cout <<filenames[0]<<endl;
          //cout <<filenames[filenames.size()-1]<<endl;
          // cleanup
          globfree(&glob_result);
        }
      // done
      }  
    }
  }  
  cout <<endl<<endl;
  //Argument searcher
  std::cout << Runlist.size() <<"\n";

  return Runlist;
}




std::vector<string> Target_list = {"Tl208","Ra224","Pb212","Bi212"};
std::vector<string> Cofig_list  = { "ExternalGammas"};
std::vector<int>    Nevts_list  = { 10000};
auto Seed_list      = arange<int>(1,251,1);


vector<nexo_runinfo> GetRunList(vector<string> list_configs = Cofig_list,
                                vector<string> list_tgts    = Target_list,
                                vector<int>    list_nevts   = Nevts_list,
                                vector<int> list_seeds      = Seed_list)
{
  vector<nexo_runinfo> Runlist;
  
  for(int i=0;i< list_configs.size(); i++)
  {
  	for(int j=0; j<list_tgts.size(); j++)
  	{
      for (int k=0; k<list_nevts.size(); k++)
  	  { 
  	    for (int m=0; m<list_seeds.size(); m++)
  	    {

  	    	nexo_runinfo RunInfo; 
  	    	cout <<list_configs[i]<<" "<<list_tgts[j]<<" "<<list_nevts[k]<<" "<<list_seeds[m]<<endl;
  	    	RunInfo.set_main(list_configs[i],"",list_tgts[j],list_nevts[k],list_seeds[m]);
  	    	Runlist.push_back(RunInfo);
        }	//End of Seed
  	  }//End of nevts
  	}//End of tgts
  }//End of configs


  return Runlist;
}

vector<nexo_runinfo> GetRunList(vector<string> list_ofconfigs,
                                vector<string> list_oftgts,
                                vector<string> list_ofnevts,      
                                vector<string> list_ofseeds )
{
  vector<nexo_runinfo> Runlist;
  
  for(int i=0;i< list_ofconfigs.size(); i++)
  {
    for(int j=0; j<list_oftgts.size(); j++)
    {
      for (int k=0; k<list_ofnevts.size(); k++)
      { 
        for (int m=0; m<list_ofseeds.size(); m++)
        {

          nexo_runinfo RunInfo; 
          //cout <<list_configs[i]<<" "<<list_tgts[j]<<" "<<list_nevts[k]<<" "<<list_seeds[m]<<endl;
          //RunInfo.set_main(list_ofconfigs[i],list_oftgts[j],list_ofnevts[k],list_ofseeds[m]);
          Runlist.push_back(RunInfo);
        } //End of Seed
      }//End of nevts
    }//End of tgts
  }//End of configs


  return Runlist;
}

vector<nexo_runinfo> GetRunList(string list_configs = "ExternalGammas_light-v2",
                                string list_tgts    = "Tl208",
                                int    list_nevts   = 5000,
                                string   seed_list_str      = "1,2,3,4")
{
  vector<nexo_runinfo> Runlist;
  vector<string> seed_list=parse_str(seed_list_str,".");

  
  for (int m=0; m<seed_list.size(); m++)
  {

    nexo_runinfo RunInfo; 
    cout <<list_configs<<" "<<list_tgts<<" "<<list_nevts<<" "<<seed_list[m]<<endl;
    RunInfo.set_main(list_configs,"",list_tgts,list_nevts,stoi(seed_list[m]));
    Runlist.push_back(RunInfo);
      } //End of Seed
  

  return Runlist;
}

vector<nexo_runinfo> GetRunList(string list_configs = "ExternalGammas_light",
                                string version = "-v2",
                                string list_tgts    = "Tl208",
                                int    list_nevts   = 5000,
                                int   first_seed    =  1,
                                int last_seed       =1)
{
  vector<nexo_runinfo> Runlist;
  //vector<string> seed_list=parse_str(seed_list_str,".");
  auto seed_list      = arange<int>(first_seed,last_seed+1,1);
  
  for (int m=0; m<seed_list.size(); m++)
  {

    nexo_runinfo RunInfo; 
    cout <<list_configs<<" "<<list_tgts<<" "<<list_nevts<<" "<<seed_list[m]<<endl;
    RunInfo.set_main(list_configs,version,list_tgts,list_nevts,seed_list[m]);
    Runlist.push_back(RunInfo);
      } //End of Seed
  

  return Runlist;
}




TChain * Load_Runlist(vector<nexo_runinfo> runlist, string treename = "Event/Sim/SimEvent")
{
  TChain *g4chain = new TChain(treename.c_str());

  for(int i=0; i<runlist.size();i++)
  {
  	string filename=runlist[i].filename;
    g4chain->Add(filename.c_str());

  }

  return g4chain;
}


TChain * Load_Recon(nexo_runinfo runinfo, string treename = "reconTree")
{
  TChain *g4chain = new TChain(treename.c_str());
  string filename=runinfo.filename_recon;
  g4chain->Add(filename.c_str());
  return g4chain;
}

TChain * Load_Sens(nexo_runinfo runinfo,string branch="", string treename = "Event/Recon")
{

  string branchname = treename + branch;
  TChain *g4chain = new TChain(branchname.c_str());
  string filename=runinfo.filename_sensrecon;

  g4chain->Add(filename.c_str());
  return g4chain;
}




TChain * Load_Runlist_Recon(vector<nexo_runinfo> runlist, string treename = "reconTree")
{
  TChain *g4chain = new TChain(treename.c_str());

  for(int i=0; i<runlist.size();i++)
  {
  	string filename=runlist[i].filename_recon;
    g4chain->Add(filename.c_str());

  }

  return g4chain;
}


/*Trees :: Event/Recon/
Energy/Energy



DNNTagroot 
*/

TChain * Load_Runlist_sens(vector<nexo_runinfo> runlist, string treename = "reconTree")
{
  TChain *g4chain = new TChain(treename.c_str());

  for(int i=0; i<runlist.size();i++)
  {
    string filename=runlist[i].filename_sensrecon;
    g4chain->Add(filename.c_str());

  }

  return g4chain;
}

//TChain * Load_all_sens(vector<nexo_runinfo> runlist)
//{
vector<nexo_runinfo> Co60rl(){
  return SearchForList("-c ExternalGammas_light -i Co60 -n 5000");
}

vector<nexo_runinfo> Cs137_Runlist(){
  return SearchForList("-c ExternalGammas_light -i Cs137 -n 5000");
}


//}



void writetofile(TChain *TC,string name="")
{

  if(name==""){ name = "Calib_analysis.root";}
  TFile *file = TFile::Open(name.c_str(),"RECREATE"); 
  TC->CloneTree(-1,"fast"); file->Write();
}



vector<vector<double>> hist_to_vector(TH1 *histo)
{
	auto Nbinsx= histo->GetNbinsX();	
	auto Nbinsy= histo->GetNbinsY();
	auto Nbinsz= histo->GetNbinsZ();

	int NumD=1;
	if(Nbinsz >1){ NumD++;}
	if(Nbinsy >1){ NumD++;}

	vector<double> bincenter;
	vector<double> value;
	for( auto zbin =0; zbin < Nbinsz; zbin++)
	{
		for( auto ybin =0; ybin < Nbinsy; ybin++)
		{
			for( auto xbin =0; xbin < Nbinsx; xbin++)
			{
				bincenter.push_back(histo->GetBinCenter(xbin));
				value.push_back(histo->GetBinContent(xbin));
			}//End of x
		}//End of y
	}//End of z

	vector<vector<double>> vect;
	vect.push_back(bincenter);
	vect.push_back(value);

	return vect;

}


void vect_to_file(vector<std::vector<double>> v, string filename="test.dat",string info="")
{
	auto N = v[0].size();
 	std::ofstream f(filename);

  if( info!=""){info = "#"+info+"\n";}
  f<<info;

	f << ",\"E\"" <<","<<"\"Count\""<<"\n";
	for(auto i =0; i<N;i++)
	{
		f<<i<<"," << v[0][i] <<","<<v[1][i]<<"\n";
	}
	f.close();

}

double  Bino_Err(double p, double n){return sqrt(p*(1-p)/(n*1.0));}

vector<double> propci_wilson_cc(double count, double nobs, double alpha=0.05)
{
    auto n = nobs;
    auto p = count;
    auto q = 1.-p;
    auto z = 2.8070337683438042;
    auto z2 = z*z;   
    auto denom = 2*(n+z2);
    auto num = 2.*n*p+z2-1.-z*sqrt(z2-2-1./n+4*p*(n*q+1));    
    auto ci_l = num/denom;
    num = 2.*n*p+z2+1.+z*sqrt(z2+2-1./n+4*p*(n*q-1));
    auto ci_u = num/denom;
    if (p == 0){ci_l = 0.;}
    else if (p == 1){ci_u = 1.;}
    vector<double> cc;
    cc.push_back(ci_l-p);
    cc.push_back(ci_u-p);

    return cc;
}

void print_calibinfo(TChain *g4tree, TChain *reconTree)
{
	double mod = (0.25+ 0.75*0.36);

	double TotalEvents= g4tree->GetEntries();

	TString Trig_cut ="fNumDeposits>0";

	double Triggered     = g4tree->GetEntries(Trig_cut);
	auto Triggered_p   = Triggered/TotalEvents/mod;
	auto Triggered_bino= Bino_Err(Triggered_p,TotalEvents);

	auto Num_recon = reconTree->GetEntries("ReconFlag==2");
    auto Num_recon_ptotal = Num_recon/TotalEvents/mod;
    auto Num_recon_ptotal_bino = Bino_Err(Num_recon_ptotal,TotalEvents);

	auto MS = reconTree->GetEntries("ReconFlag==2 & NumReconClusters>1");
    auto MS_ptotal = MS/TotalEvents/mod;
    auto MS_ptotal_bino=Bino_Err(MS_ptotal,TotalEvents);

    auto SS = reconTree->GetEntries("ReconFlag==2 & NumReconClusters==1");
    auto SS_ptotal = SS/TotalEvents/mod;
    auto SS_ptotal_bino=Bino_Err(SS_ptotal,TotalEvents);
    
    auto SS_peak = reconTree->GetEntries("ReconFlag==2 & NumReconClusters==1 & TotalReconEnergy*1000>=2610 & TotalReconEnergy*1000<2620");
    auto SS_peak_ptotal = SS_peak/TotalEvents/mod;
    auto SS_peak_ptotal_bino=Bino_Err(SS_peak_ptotal,TotalEvents);
    
	auto SS_peak_ws = propci_wilson_cc(SS_peak_ptotal,TotalEvents);


    auto SS_deep = reconTree->GetEntries("ReconFlag==2 & NumReconClusters==1 & StandoffDistance>=203");
    auto SS_deep_ptotal = SS_deep/TotalEvents/mod;
    auto SS_deep_ptotal_bino=Bino_Err(SS_deep_ptotal,TotalEvents);
    
    auto SS_deep_ws = propci_wilson_cc(SS_deep_ptotal,TotalEvents);

    cout <<"SS Deep "<< SS_deep <<endl;
    auto SS_dpeak = reconTree->GetEntries("ReconFlag==2 & NumReconClusters==1 & StandoffDistance>=203 & TotalReconEnergy*1000>=2610 & TotalReconEnergy*1000<2620");
    auto SS_dpeak_ptotal = SS_dpeak/TotalEvents/mod;
    auto SS_dpeak_ptotal_bino=Bino_Err(SS_dpeak_ptotal,TotalEvents);
    
    auto SS_dpeak_ws = propci_wilson_cc(SS_dpeak_ptotal,TotalEvents);


    auto  Deeps = reconTree->GetEntries("StandoffDistance>=203");

    cout <<" SS deep peak "<< SS_dpeak <<endl;

    printf( "Number of Disintegrations :: %f\n",TotalEvents);
    string trig=Form("Fraction of Triggered Events \t\t %0.4f +- %0.4f",Triggered_p, Triggered_bino);
    string rec =Form("Fraction of Reconstructed Events \t %0.4f +- %0.4f",Num_recon_ptotal,Num_recon_ptotal_bino);
    string MS_st =Form("Fraction of MS Events \t\t\t %0.4f +- %0.5f",MS_ptotal,MS_ptotal_bino);
    string SS_st =Form("Fraction of SS Events  \t\t\t %0.4f +- %0.5f",SS_ptotal,SS_ptotal_bino);
    string SS_peak_st =Form("Fraction of SS Peak E Events \t\t %0.2e +- %0.2e , wilson score %.2e,%.2e",SS_peak_ptotal,SS_peak_ptotal_bino,SS_peak_ws[0],SS_peak_ws[1]);
    string SS_deep_st =Form("Fraction of SS deep Events \t\t %0.2e +- %0.2e , wilson score %.2e,%.2e",SS_deep_ptotal,SS_deep_ptotal_bino,SS_deep_ws[0],SS_deep_ws[1]);
    string SS_dpeak_st =Form("Fraction of SS deep peek E Events \t %0.2e +- %0.2e, wilson score %.2e,%.2e",SS_dpeak_ptotal,SS_dpeak_ptotal_bino, SS_dpeak_ws[0], SS_dpeak_ws[1]);
    printf("%s\n",trig.c_str());
    printf("%s\n",rec.c_str());
    printf("%s\n",MS_st.c_str());
    printf("%s\n",SS_st.c_str());
    printf("%s\n",SS_peak_st.c_str());
    printf("%s\n",SS_deep_st.c_str());
    printf("%s\n",SS_dpeak_st.c_str());

}





void calib_results(vector<string> list_configs=Cofig_list, vector<string> list_tgts=Target_list, vector<int> list_nevts=Nevts_list,vector<int> list_seeds=Seed_list, std::string file_name="")
{
	string fname;
	if( file_name=="")
	{fname="Th228";}
	else{ fname=file_name;}


	vector<nexo_runinfo> Runlist;
	Runlist = GetRunList(list_configs,list_tgts,list_nevts,list_seeds); 

	print_calibinfo(Load_Runlist(Runlist),Load_Runlist_Recon(Runlist));

	TH1F *MS = new TH1F("MS","MultiSite",30000,0,3000);
	TH1F *SS = new TH1F("SS","SingleSite",30000,0,3000);
	TChain *recontree = Load_Runlist_Recon(Runlist);
	recontree->Draw("TotalReconEnergy*1000>>MS","ReconFlag==2 & NumReconClusters>1" ,"");
	recontree->Draw("TotalReconEnergy*1000>>SS","ReconFlag==2 & NumReconClusters==1","");
	vector<vector<double>> MS_vect =hist_to_vector(MS);  
	vect_to_file(MS_vect,Form("%s_MS.csv",fname.c_str()));
	vector<vector<double>> SS_vect =hist_to_vector(SS);  
	vect_to_file(SS_vect,Form("%s_SS.csv",fname.c_str()));

}

void calib_results(vector<nexo_runinfo> Runlist, std::string file_name="")
{
  string fname;
  if( file_name=="")
  {fname="Th228";}
  else{ fname=file_name;}
  double a=0.0;
  for( auto i =0; i<Runlist.size();i++){a+=Runlist[i].nevts;}

  //print_calibinfo(Load_Runlist(Runlist),Load_Runlist_Recon(Runlist));

  TH1F *MS = new TH1F("MS","MultiSite",30000,0,3000);
  TH1F *SS = new TH1F("SS","SingleSite",30000,0,3000);
  TChain *recontree = Load_Runlist_Recon(Runlist);
  int numofevents = recontree->GetEntries();

  recontree->Draw("TotalReconEnergy*1000>>MS",Form("(ReconFlag==2 & NumReconClusters>1) /%f",a) ,"");
  recontree->Draw("TotalReconEnergy*1000>>SS",Form("(ReconFlag==2 & NumReconClusters==1) /%f",a),"");
  vector<vector<double>> MS_vect =hist_to_vector(MS);  
  vect_to_file(MS_vect,Form("%s_MS.csv",fname.c_str()),Form("%.1f,%i",a,numofevents));
  vector<vector<double>> SS_vect =hist_to_vector(SS);  
  vect_to_file(SS_vect,Form("%s_SS.csv",fname.c_str()),Form("%.1f,%i",a,numofevents));



  TChain *SimTree = Freind_all(Runlist);
  setcut();
  TCut Sens_tcut = Form("%s",SensCut.c_str());

  TH1F *MS_sens = new TH1F("MS_sens","MultiSite",30000,0,3000);
  TH1F *SS_sens = new TH1F("SS_sens","SingleSite",30000,0,3000);
  SimTree->Draw("TotalReconEnergy*1000>>SS_sens",Sens_tcut,"");
  vector<std::vector<double>> SS_sens_vect=hist_to_vector(SS_sens);
  vect_to_file(SS_sens_vect,Form("%s_sense_SS.csv",fname.c_str()),Form("%f",a));


}


void Run_all_light()
{

vector<vector<string>> Target_list ={{"Tl208","Bi212","Ra224","Pb212"}};


vector<string> Config ;




}



/*
void nexo_runinfo::set_main(string config,string isotope,int nevts,int seed) 
{
	cout << "poop "<< isotope << "  " <<abc <<endl;

  config=config;
  isotope=isotope;
  seed=seed;
  nevts=nevts;

}*/

void save_tchain(TChain *T, string name="TChain.root")
{
  TFile *file = TFile::Open(name.c_str(),"RECREATE");
  T->CloneTree(-1,"fast"); 
  file->Write(); delete file;
}


TChain *Freind_all(vector<nexo_runinfo> Runlist)
{
  //'/Event/Recon/Standoff/Standoff'
  //'/Event/Recon/Energy/Energy'
  //'/Event/Recon/ChargeQuanta/ChargeQuanta'
  //'/Event/Recon/Photons/Photons'
  //Info for tree naming
  string TreePath="Event/Recon";
  vector<string> Treenames={"Standoff","Energy","ChargeQuanta","Photons","DNNTag"};


  TChain *G4   =Load_Runlist(Runlist);
  G4->AddFriend(Load_Runlist_Recon(Runlist));

  for( auto i =0;i<Treenames.size();i++)
  {
    G4->AddFriend(Load_Runlist_sens(Runlist,Form("%s/%s/%s",TreePath.c_str(),Treenames[i].c_str(),Treenames[i].c_str())) );

  }
  setcut();
  return G4;
}

TChain *Freind_sens(vector<nexo_runinfo> Runlist)
{
  //'/Event/Recon/Standoff/Standoff'
  //'/Event/Recon/Energy/Energy'
  //'/Event/Recon/ChargeQuanta/ChargeQuanta'
  //'/Event/Recon/Photons/Photons'
  //Info for tree naming
  string TreePath="Event/Recon";
  vector<string> Treenames={"Standoff","Energy","ChargeQuanta","Photons","DNNTag"};


  TChain *G4   =Load_Runlist_Recon(Runlist);

  for( auto i =0;i<Treenames.size();i++)
  {
    G4->AddFriend(Load_Runlist_sens(Runlist,Form("%s/%s/%s",TreePath.c_str(),Treenames[i].c_str(),Treenames[i].c_str())) );

  }
  setcut();
  return G4;
}




vector<nexo_runinfo> Ra226_Runlist(string args="")
{
  vector<nexo_runinfo> RList;
  nexo_args Args = Get_args(args);
  Args.list_of_tgts=parse_str("Ra226 Pb214 Bi214 Bi210 Pb206");
  if (Args.list_of_configs.size()==0)Args.list_of_configs={"ExternalGammas_light"};
  if(Args.list_of_seeds.size()==0)Args.list_of_seeds=arange(1,500,1);
  //if(Args.list_of_nevts.size()==0)Args.list_of_nevts.push_back(5000);

  cout << Args.list_of_configs[0] <<endl;
  vector<nexo_runinfo> templist;
  for(auto i=0;i<Args.list_of_tgts.size();i++)
  {
    string tgt=Args.list_of_tgts[i];
    templist = SearchForList(Form("-c ExternalGammas_light -i %s -n 5000",tgt.c_str()));
    RList.insert(end(RList),begin(templist),end(templist));
    templist.clear();
  }
  //RList = GetRunList(Args.list_of_configs, Args.list_of_tgts, Args.list_of_nevts, Args.list_of_seeds);


  return RList;

}

vector<nexo_runinfo> Th228_Runlist(string args="")
{
  vector<nexo_runinfo> RList;
  nexo_args Args = Get_args(args);
  Args.list_of_tgts=parse_str("Th228 Tl208 Ra224 Pb212 Bi212");
  if (Args.list_of_configs.size()==0)Args.list_of_configs={"ExternalGammas_light"};
  if(Args.list_of_seeds.size()==0)Args.list_of_seeds=arange(1,500,1);
  //if(Args.list_of_nevts.size()==0)Args.list_of_nevts.push_back(5000);

  cout << Args.list_of_configs[0] <<endl;
  vector<nexo_runinfo> templist;
  for(auto i=0;i<Args.list_of_tgts.size();i++)
  {
    string tgt=Args.list_of_tgts[i];
    templist = SearchForList(Form("-c ExternalGammas_light -i %s -n 5000",tgt.c_str()));
    RList.insert(end(RList),begin(templist),end(templist));
    templist.clear();
  }
  //RList = GetRunList(Args.list_of_configs, Args.list_of_tgts, Args.list_of_nevts, Args.list_of_seeds);

  //a.insert(std::end(a), std::begin(b), std::end(b));


  return RList;

}


int SimTotal(vector<nexo_runinfo> Runlist)
{
  int TotalEvents=0;

  for ( auto i =0; i<Runlist.size();i++)
  {
    TotalEvents+=Runlist[i].nevts;
  }


  return TotalEvents;
}



void weight_recon(string filename)
{
  float Weight_factor;


  TFile* file = new TFile (filename.c_str(), "update"); 
  TTree* t1 = (TTree*)file->Get("reconTree");
  if ( t1->FindBranch("reconWeightFactor")!=nullptr){return;}
  if( filename.find("Tl208")>0){Weight_factor=0.36;}
  else{ Weight_factor=1.0;}

  TBranch *newBranch = t1->Branch("reconWeightFactor", &Weight_factor, "reconWeightFactor/F");
  for (Int_t i = 0; i < t1->GetEntries(); i++)
    {newBranch->Fill();}
  
  t1->Write("", TObject::kOverwrite);

}

void addSeed(nexo_runinfo RI)
{
  int Seed = RI.seed;
  string filename = RI.filename_recon;
  cout << filename <<" "<< Seed <<endl;

  TFile* file = new TFile (filename.c_str(), "update"); 
  TTree* t1 = (TTree*)file->Get("reconTree");
  //return;
  //cout << t1->FindBranch("Seed") <<endl;
  if ( t1->FindBranch("Seed")!=nullptr){
    cout << "Seed exist" <<endl;
    //TBranch *b = t1->GetBranch("Seed");   
    //cout <<"Remove branch"<<endl;
    //t1->GetListOfBranches()->Remove(b);
    return;
  }
  TBranch *newBranch = t1->Branch("Seed", &Seed, "Seed/I");
  for (Int_t i = 0; i < t1->GetEntries(); i++)
    {newBranch->Fill();}
  
  t1->Write("", TObject::kOverwrite);
  
}



void addSeed(string filename, int Seed)
{
 

  TFile* file = new TFile (filename.c_str(), "update"); 
  TTree* t1 = (TTree*)file->Get("reconTree");
  
  if ( t1->FindBranch("Seed")!=nullptr){return;}
  TBranch *newBranch = t1->Branch("Seed", &Seed, "Seed/I");
  for (Int_t i = 0; i < t1->GetEntries(); i++)
    {newBranch->Fill();}
  
  t1->Write("", TObject::kOverwrite);
  
}


void AnalysisCo60()
{
  vector<nexo_runinfo> RL= Co60rl();
  TChain *Recon = Load_Runlist_Recon(RL);
  calib_results(RL,"Co60");
  save_tchain(Recon,"Co60_recon.root");

  TChain *SR = Freind_sens(RL);
  save_tchain(SR,"Co60_sens.root");


}

void AnalysisTh228()
{
  vector<nexo_runinfo> RL= Th228_Runlist();
  TChain *Recon = Load_Runlist_Recon(RL);
  calib_results(RL,"Th228");
  save_tchain(Recon,"Th228_recon.root");

  TChain *SR = Freind_sens(RL);
  save_tchain(SR,"Th228_sens.root");

}


void AnalysisRa226()
{
  vector<nexo_runinfo> RL= Ra226_Runlist();
  TChain *Recon = Load_Runlist_Recon(RL);
  calib_results(RL,"Ra226");
  save_tchain(Recon,"Ra226_recon.root");

  TChain *SR = Freind_sens(RL);
  save_tchain(SR,"Ra226_sens.root");

}

void AnalysisCs137()
{
  vector<nexo_runinfo> RL= Cs137_Runlist();
  TChain *Recon = Load_Runlist_Recon(RL);
  calib_results(RL,"Cs137");
  save_tchain(Recon,"Cs137_recon.root");

  TChain *SR = Freind_sens(RL);
  save_tchain(SR,"Cs137_sens.root");

}


void joinsenstrees(nexo_runinfo RI){
string filename=RI.filename_sensrecon;

TFile oldfile(filename.c_str());
//TTree Sens("Sens","sens");

vector<string> Treenames={"Standoff","Energy","ChargeQuanta","Photons","DNNTag"};


TTree * standoff;
oldfile.GetObject("Standoff",standoff);
TTree * energy;
oldfile.GetObject("Energy",energy);
TTree * CQ;
oldfile.GetObject("ChargeQuanta",CQ);
TTree * photons;
oldfile.GetObject("Photons",photons);
TTree * DNN;
oldfile.GetObject("DNNTag",DNN);

string tgt = RI.isotope;
int seed= RI.seed;

TFile newfile(Form("/home/wrkshp/nexo/analysis/calib_ana/%s_%d.root",tgt.c_str(),seed));


auto Sens=standoff->CloneTree();
Sens->AddFriend(energy->CloneTree());

Sens->Print();
newfile.Write();


}


TChain * Senscomb(vector<nexo_runinfo> runlist, string version ="",string newfilename="")
{

  if(newfilename.find("/")==std::string::npos){newfilename=simple_sens_path +"/"+ newfilename;}


 TChain *nexo= new TChain("nexo");

  for(int i=0; i<runlist.size();i++)
  {
    string filename=runlist[i].filename_simplesens;
    
    if(exists_test3(filename)){
cout << filename<< " "<<exists_test3(filename)<<endl;
      nexo->Add(filename.c_str());
      cout << nexo->GetEntries()<< endl;}
    else{
      cout<< "Missing ::" <<filename<<endl;
    }

  }


cout <<"Finished Loop: Opening new file"<<endl;

 TFile *file = TFile::Open(newfilename.c_str(),"RECREATE");
 cout << "Recreated new file" <<endl;

 nexo->CloneTree(-1,"fast"); 
 cout << "Cloned Tree"<<endl;
 file->Write(); 

 cout <<"writen to :"<< newfilename<<endl;
delete file;
cout <<"Deleted"<<endl;
  if(exists_test3(newfilename))
    {cout <<"no issue"<<endl;}
return nexo;
}



/*
void SensAna_all()
{


  for(auto i=0;i<)
  {
    filename=
    Senscomb(RL,filename);

  }
}
*/





