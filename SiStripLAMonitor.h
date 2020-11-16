#ifndef CalibTracker_SiStripLorentzAngle_SiStripLAMonitor_h
#define CalibTracker_SiStripLorentzAngle_SiStripLAMonitor_h 1

#include <string>
#include <vector>

#include "TChain.h"
#include "TH1.h" 
#include "TH2.h"
#include "TProfile.h"

int Init(int argc, char * argv[]);
void AnalyzeTheTree(std::map<int, double> LA_dict);
void ProcessTheEvent(std::map<int, double> LA_dict);
void ProcessTheModule(const unsigned int &, std::map<int, double> LA_dict);
void WriteOutputs(const bool & savehistos = true);
void WriteOutputsModules(TFile &, TH2F *);

void CalibTreeBranches(TTree *);
void InfoTreeBranches(TTree *);
std::string ModuleLocationType(const unsigned int &);
// void MeasureLorentzAngle(std::string method);

// number of events
int nevents_;

// B field
std::string bfield_;

// la info tree
std::vector<unsigned int> * inforawid_               = nullptr;
std::vector<float>        * infoglobalZofunitlocalY_ = nullptr;
std::vector<float>        * infolocalb_              = nullptr;
std::vector<float>        * infola_                  = nullptr;

std::map<unsigned int,int> orientation_;
std::map<unsigned int,float> la_db_;

// event data
unsigned int eventnumber_ = 0;
unsigned int runnumber_   = 0;
unsigned int luminumber_  = 0;

// calib data
std::vector<unsigned int>   * trackindex_ = nullptr;
std::vector<unsigned int>   * rawid_      = nullptr;
std::vector<unsigned short> * nstrips_    = nullptr;
std::vector<float>          * localdirx_  = nullptr;
std::vector<float>          * localdiry_  = nullptr;
std::vector<float>          * localdirz_  = nullptr;
std::vector<float>          * variance_   = nullptr;
std::vector<float>          * localx_     = nullptr;
std::vector<float>          * rhlocalx_   = nullptr;

// track data
std::vector<float>        * trackpt_        = nullptr;
std::vector<float>        * tracketa_       = nullptr;
std::vector<float>        * trackphi_       = nullptr;
std::vector<unsigned int> * trackhitsvalid_ = nullptr;
std::vector<float>        * trackchi2ndof_  = nullptr;

// histogramming
std::map<std::string, TH1F*> h1_;
std::map<std::string, TH2F*> h2_;

std::map<int, TH2F*> h2_ct_w_m_;
std::map<int, TH2F*> h2_ct_var2_m_;
std::map<int, TH2F*> h2_ct_var3_m_;

std::map<int, TH2F*> h2_t_w_m_;
std::map<int, TH2F*> h2_t_var2_m_;
std::map<int, TH2F*> h2_t_var3_m_;

std::map<std::string, TProfile*> hp_;

// info
std::map<std::string,int> nlayers_;
std::vector<std::string> modtypes_;

std::map<std::string,float> la_;



#endif
