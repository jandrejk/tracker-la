#include <iostream>
#include <fstream> 
#include <regex>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;


#include "SiStripLAMonitor.h"
#include "SiStripLAMonitorConfig.h"

#include "TOBDetId.h"
#include "TIBDetId.h"

#include "TFile.h" 
#include "TFileCollection.h"
#include "TVector3.h"
#include "TString.h"

#include <sstream>
#include <string>
#include <vector>

std::map<int, double> processList(std::string path){

   

   std::map<int, double> LA_dict;

   if (path.empty()) {
      return LA_dict;
   } 

   std::ifstream openFile(path);

   std::string keyString;
   std::string valueString;
   
   int counter = 0;
   while(true)
   {  
         
         if (!std::getline(openFile, keyString, ',')) {
            std::cout<<"csv file loaded"<<std::endl;
            break;
         }
         // std::cout<<"keyString: "<<keyString<<std::endl;
         std::getline(openFile, valueString, '\n');
         auto key = std::stoi(keyString);
         auto value = std::stod(valueString);
         LA_dict[key] = value;
         counter++;
         // if (counter == 10) {
         //    break;
         // }
   }
   std::cout<<"no of entries in csv file: "<<counter<<std::endl;
   return LA_dict;
}

// =============================================================================================   

int main(int argc, char * argv[])
{
   if ( Init(argc,argv) == -1 )
   {
      std::cout << "*** SiStripLAMonitor ***: -errors- Please check your configuration file" << std::endl;
      return -1;
   }
   std::map<int, double> LA_dict = processList(LA_csv_file_);
   AnalyzeTheTree(LA_dict);
   WriteOutputs(saveHistos_);
   
   std::cout << "-----" << std::endl;
   for ( auto & la : la_ )
   {
      std::cout << la.first << "  " << la.second << std::endl;
   }
   std::cout << "-----" << std::endl;
   
   std::cout << "SiStripLAMonitor finished!" << std::endl;
   
   std::ofstream finished;
   finished.open("finished.txt");
   finished << "SiStripLAMonitor finished!" << "\n";
   finished.close();
   
   
   return 0;
}

int Init(int argc, char * argv[])
{
   // read configuration
   if ( SiStripLAMonitorConfig(argc, argv) != 0 ) return -1;
   
   h1_["track_pt"]        = new TH1F ("track_pt","", 2000,0,1000);
   h1_["track_eta"]       = new TH1F ("track_eta","", 100,-4,4);
   h1_["track_phi"]       = new TH1F ("track_phi","", 80,-3.2,3.2);
   h1_["track_validhits"] = new TH1F ("track_validhits","", 50,0,50);
   h1_["track_chi2ndof"]  = new TH1F ("track_chi2ndof","", 100,0,5);
   h2_["track_chi2xhits"] = new TH2F ("track_chi2xhits_2d","", 100,0,5,50,0,50);
   h2_["track_ptxhits"]   = new TH2F ("track_ptxhits_2d","", 200,0,100,50,0,50);
   h2_["track_etaxhits"]  = new TH2F ("track_etaxhits_2d","", 60,-3,3,50,0,50);
   h2_["track_ptxchi2"]   = new TH2F ("track_ptxchi2_2d","", 200,0,100,100,0,5);
   h2_["track_ptxeta"]    = new TH2F ("track_ptxeta_2d","", 200,0,100,60,-3,3);
   h2_["track_etaxchi2"]  = new TH2F ("track_etaxchi2_2d","", 60,-3,3,100,0,5);
   
   //
   nlayers_["TIB"] = 4;
   nlayers_["TOB"] = 6;
   modtypes_.push_back("s");
   modtypes_.push_back("a");
   
   // prepare histograms
   for ( auto & layers : nlayers_)
   {
      std::string subdet = layers.first;
      for ( int l = 1; l <= layers.second; ++l )
      {
         for ( auto & t : modtypes_ )
         {
            std::string locationtype = Form("%s_L%d%s",subdet.c_str(),l,t.c_str());
            //std::cout << "preparing histograms for " << locationtype << std::endl;
            h1_[Form("%s_nstrips"     ,locationtype.c_str())]  = new TH1F (Form("%s_nstrips",locationtype.c_str()),     "", 20,0,20);
            h1_[Form("%s_tanthetatrk" ,locationtype.c_str())]  = new TH1F (Form("%s_tanthetatrk",locationtype.c_str()), "", 300,-1.5,1.5);
            h1_[Form("%s_cosphitrk"   ,locationtype.c_str())]  = new TH1F (Form("%s_cosphitrk",locationtype.c_str()), "", 40,-1,1);
            h1_[Form("%s_variance_w2" ,locationtype.c_str())]  = new TH1F (Form("%s_variance_w2",locationtype.c_str()),     "", 100,0,1);
            h1_[Form("%s_variance_w3" ,locationtype.c_str())]  = new TH1F (Form("%s_variance_w3",locationtype.c_str()),     "", 100,0,1);
            
            // local_x
            h1_[Form("%s_localx_vplus"       ,locationtype.c_str())]  = new TH1F (Form("%s_localx_vplus",locationtype.c_str()),     "", 200,-6,6);
            h1_[Form("%s_rhlocalx_vplus"     ,locationtype.c_str())]  = new TH1F (Form("%s_rhlocalx_vplus",locationtype.c_str()),     "", 200,-6,6);
            h1_[Form("%s_shift_vplus"        ,locationtype.c_str())]  = new TH1F (Form("%s_shift_vplus",locationtype.c_str()),     "", 200,-0.05,0.05);
            h2_[Form("%s_shift_tanthetatrk_minus_tanthetaLA_vplus",locationtype.c_str())] = new TH2F (Form("%s_shift_tanthetatrk_minus_tanthetaLA_vplus",locationtype.c_str()), "",12, -0.5, 0.5, 400, -0.01, 0.01);
            h2_[Form("%s_clusterwidth_tanthetatrk_minus_tanthetaLA_vplus",locationtype.c_str())] = new TH2F (Form("%s_clusterwidth_tanthetatrk_minus_tanthetaLA_vplus",locationtype.c_str()), "",12, -0.5, 0.5, 11, 0, 10);
            
            h1_[Form("%s_localx_vminus"       ,locationtype.c_str())]  = new TH1F (Form("%s_localx_vminus",locationtype.c_str()),     "", 200,-6,6);
            h1_[Form("%s_rhlocalx_vminus"     ,locationtype.c_str())]  = new TH1F (Form("%s_rhlocalx_vminus",locationtype.c_str()),     "", 200,-6,6);
            h1_[Form("%s_shift_vminus"        ,locationtype.c_str())]  = new TH1F (Form("%s_shift_vminus",locationtype.c_str()),     "", 200,-0.05,0.05);
            h2_[Form("%s_shift_tanthetatrk_minus_tanthetaLA_vminus",locationtype.c_str())] = new TH2F (Form("%s_shift_tanthetatrk_minus_tanthetaLA_vminus",locationtype.c_str()), "",12, -0.5, 0.5, 400, -0.01, 0.01);
            h2_[Form("%s_clusterwidth_tanthetatrk_minus_tanthetaLA_vminus",locationtype.c_str())] = new TH2F (Form("%s_clusterwidth_tanthetatrk_minus_tanthetaLA_vminus",locationtype.c_str()), "",12, -0.5, 0.5, 11, 0, 10);
            

            h2_[Form("%s_tanthcosphtrk_nstrip",locationtype.c_str())] = new TH2F (Form("%s_tanthcosphtrk_nstrip",locationtype.c_str()), "", 360, -0.9, 0.9, 20, 0, 20);
            h2_[Form("%s_thetatrk_nstrip",locationtype.c_str())]      = new TH2F (Form("%s_thetatrk_nstrip",locationtype.c_str())     , "", 360, -0.9, 0.9, 20, 0, 20);
            h2_[Form("%s_tanthcosphtrk_var2",locationtype.c_str())]   = new TH2F (Form("%s_tanthcosphtrk_var2",locationtype.c_str())  , "", 360, -0.9, 0.9, 50, 0,  1);
            h2_[Form("%s_tanthcosphtrk_var3",locationtype.c_str())]   = new TH2F (Form("%s_tanthcosphtrk_var3",locationtype.c_str())  , "", 360, -0.9, 0.9, 50, 0,  1);
            h2_[Form("%s_thcosphtrk_var2",locationtype.c_str())]      = new TH2F (Form("%s_thcosphtrk_var2",locationtype.c_str())     , "", 360, -0.9, 0.9, 50, 0,  1);
            h2_[Form("%s_thcosphtrk_var3",locationtype.c_str())]      = new TH2F (Form("%s_thcosphtrk_var3",locationtype.c_str())     , "", 360, -0.9, 0.9, 50, 0,  1);
         }
      }
   }
   
   // histograms per module require info avaialble in the info tree, namely the module id
   // in principle could use all modules but where can I find them easily?
   
   return 0;
}


void ProcessTheEvent(std::map<int, double> LA_dict)
{
   
   bool trk_done[1000] = { false };
   
   
   for ( size_t i = 0 ; i < rawid_->size(); ++i ) // loop over modules
   {
      // do whatever pre-selection needed
      int itrk = trackindex_->at(i);
      
   
      if ( itrk >= 1000 )
      {
         std::cout << "SiStripLAMonitor::ProcessTheEvent() *** WARNING *** More than 1000 tracks!!! Skipping!" << std::endl;
         break;
      }
      
      if ( ptmin_ > 0 && infolocalb_->at(0) > 1. && trackpt_->at(itrk) < ptmin_ ) continue; // if NO Bfield don't perform any pt selection 
      if ( ptmax_ > 0 && infolocalb_->at(0) > 1. && trackpt_->at(itrk) >= ptmax_ ) continue; // if NO Bfield don't perform any pt selection
      if ( tracketa_->at(itrk) < etamin_ ) continue;
      if ( tracketa_->at(itrk) >= etamax_ ) continue;
      if ( hitsvalmin_ > 0 && int(trackhitsvalid_->at(itrk)) < hitsvalmin_ ) continue;
      if ( hitsvalmax_ > 0 && int(trackhitsvalid_->at(itrk)) >= hitsvalmax_ ) continue;
      if ( chi2ndfmax_ > 0 && trackchi2ndof_->at(itrk) >= chi2ndfmax_ ) continue;
      if ( chi2ndfmin_ > 0 && trackchi2ndof_->at(itrk) < chi2ndfmin_ ) continue;
      
      // std::cout << trk_done[itrk] << std::endl;
      if ( ! trk_done[itrk] )
      {
         h1_["track_pt"]          -> Fill(trackpt_->at(itrk));
         h1_["track_eta"]         -> Fill(tracketa_->at(itrk));
         h1_["track_phi"]         -> Fill(trackphi_->at(itrk));
         h1_["track_validhits"]   -> Fill(trackhitsvalid_->at(itrk));
         h1_["track_chi2ndof"]    -> Fill(trackchi2ndof_->at(itrk));
         h2_["track_chi2xhits"]   -> Fill(trackchi2ndof_->at(itrk),trackhitsvalid_->at(itrk));
         h2_["track_ptxhits"]     -> Fill(trackpt_->at(itrk),trackhitsvalid_->at(itrk));
         h2_["track_etaxhits"]    -> Fill(tracketa_->at(itrk),trackhitsvalid_->at(itrk));
         h2_["track_ptxchi2"]     -> Fill(trackpt_->at(itrk),trackchi2ndof_->at(itrk));
         h2_["track_ptxeta"]      -> Fill(trackpt_->at(itrk),tracketa_->at(itrk));
         h2_["track_etaxchi2"]    -> Fill(tracketa_->at(itrk),trackchi2ndof_->at(itrk));
         
         trk_done[itrk] = true;
      }
      
//      if ( fabs(tracketa_->at(itrk)) > 0.2 ) continue;
      // process info
      // std::cout << "tracketa: " << tracketa_->at(itrk) << std::endl;
      ProcessTheModule(i, LA_dict);
   }
   
}


void ProcessTheModule(const unsigned int & i, std::map<int, double> LA_dict)
{
   // std::map<int, double> LA_dict = processList(LA_csv_file_);
   
   unsigned int mod = rawid_->at(i);
   std::string locationtype = ModuleLocationType(mod);
   
   if ( locationtype == "" ) return;
   
   la_[locationtype] = la_db_[mod];
   // std::cout << "locationtype: " << locationtype << std::endl;
   // std::cout << "mod: " << mod << std::endl;
   // std::cout << "la_db_ mod: " << la_db_[mod] << std::endl;
   // std::cout << "magnetic field: " << infolocalb_->at(0) << std::endl;
   
   TVector3 localdir(localdirx_->at(i),localdiry_->at(i),localdirz_->at(i));
   int sign = orientation_[mod];
   float tantheta = TMath::Tan(localdir.Theta());
   float cosphi   = TMath::Cos(localdir.Phi());
   float theta    = localdir.Theta();

   // std::cout << "x: " << localdirx_->at(i) << std::endl;
   // std::cout << "y: " << localdiry_->at(i) << std::endl;
   // std::cout << "z: " << localdirz_->at(i) << std::endl;
   // std::cout << "theta: " << theta << std::endl;
   
   
   unsigned short nstrips  = nstrips_->at(i);
   float variance = variance_->at(i);

   // std::cout << "nstrips: " << nstrips << std::endl;
   
   h1_[Form("%s_nstrips"    ,locationtype.c_str())] -> Fill(nstrips);
   h1_[Form("%s_tanthetatrk",locationtype.c_str())] -> Fill(sign*tantheta);
   h1_[Form("%s_cosphitrk"  ,locationtype.c_str())] -> Fill(cosphi);
   
   // addition for backplane correction: localx and rhlocalx 
   float tanLA_used = -1;
   if (infolocalb_->at(0) > 0.1) { // non-zero magnetic field --> use LA from data base
      // multiply LA from data base by B-field.
      // std::cout << "module ID: " << mod << std::endl;
      // std::cout << "corresponding LA: " << LA_dict[mod] << std::endl;
      if (LA_dict.empty()) { // use LA from DB
         tanLA_used = TMath::Tan(la_db_[mod] * infolocalb_->at(0));
      }
      else { // use LA from csv file
         tanLA_used = TMath::Tan(LA_dict[mod] * infolocalb_->at(0));
      }
   }
   else { // zero magnetic field --> LA = 0
      tanLA_used = TMath::Tan(0.0);
   }
   if (sign >= 0) {
      h1_[Form("%s_localx_vplus"    ,locationtype.c_str())] -> Fill(localx_->at(i));
      h1_[Form("%s_rhlocalx_vplus"    ,locationtype.c_str())] -> Fill(rhlocalx_->at(i));
      h1_[Form("%s_shift_vplus"       ,locationtype.c_str())] -> Fill(localx_->at(i) - rhlocalx_->at(i));
      h2_[Form("%s_shift_tanthetatrk_minus_tanthetaLA_vplus",locationtype.c_str())] -> Fill(sign*(tantheta-tanLA_used),localx_->at(i) - rhlocalx_->at(i));
      h2_[Form("%s_clusterwidth_tanthetatrk_minus_tanthetaLA_vplus",locationtype.c_str())] -> Fill(sign*(tantheta-tanLA_used),nstrips);
   }
   else {
      h1_[Form("%s_localx_vminus"    ,locationtype.c_str())] -> Fill(localx_->at(i));
      h1_[Form("%s_rhlocalx_vminus"    ,locationtype.c_str())] -> Fill(rhlocalx_->at(i));
      h1_[Form("%s_shift_vminus"       ,locationtype.c_str())] -> Fill(localx_->at(i) - rhlocalx_->at(i));
      h2_[Form("%s_shift_tanthetatrk_minus_tanthetaLA_vminus",locationtype.c_str())] -> Fill(sign*(tantheta-tanLA_used),localx_->at(i) - rhlocalx_->at(i)); 
      h2_[Form("%s_clusterwidth_tanthetatrk_minus_tanthetaLA_vminus",locationtype.c_str())] -> Fill(sign*(tantheta-tanLA_used),nstrips);
   }
   // nstrips
   h2_[Form("%s_tanthcosphtrk_nstrip",locationtype.c_str())] -> Fill(sign*cosphi*tantheta,nstrips);
   h2_[Form("%s_thetatrk_nstrip",locationtype.c_str())]      -> Fill(sign*theta*cosphi,nstrips);
   
   // variance for width == 2
   if ( nstrips == 2 )
   {
      h1_[Form("%s_variance_w2"       ,locationtype.c_str())]    -> Fill(variance);
      h2_[Form("%s_tanthcosphtrk_var2",locationtype.c_str())]    -> Fill(sign*cosphi*tantheta,variance);
      h2_[Form("%s_thcosphtrk_var2"   ,locationtype.c_str())]    -> Fill(sign*cosphi*theta,variance);
      if ( saveHistosMods_ )
      {
         h2_ct_var2_m_[mod] -> Fill(sign*cosphi*tantheta,variance);
         h2_t_var2_m_[mod]  -> Fill(sign*cosphi*theta,variance);
      }
   }
   // variance for width == 3
   if ( nstrips == 3 )
   {
      h1_[Form("%s_variance_w3"       ,locationtype.c_str())]    -> Fill(variance);
      h2_[Form("%s_tanthcosphtrk_var3",locationtype.c_str())]    -> Fill(sign*cosphi*tantheta,variance);
      h2_[Form("%s_thcosphtrk_var3"   ,locationtype.c_str())]    -> Fill(sign*cosphi*theta,variance);
      if ( saveHistosMods_ )
      {
         h2_ct_var3_m_[mod] -> Fill(sign*cosphi*tantheta,variance);
         h2_t_var3_m_[mod]  -> Fill(sign*cosphi*theta,variance);
      }
   }
   
   if ( saveHistosMods_ )
   {
      h2_ct_w_m_[mod] -> Fill(sign*cosphi*tantheta,nstrips);
      h2_t_w_m_[mod]  -> Fill(sign*cosphi*theta,nstrips);
   }
   
   
}

void AnalyzeTheTree(std::map<int, double> LA_dict)
{
   // std::cout << "test 0" << std::endl;
   int count_entries = 0;
   bool terminate = false;
   fs::path calibtree_path(calibTreeDir_);
   std::cout << "path: " << calibtree_path << std::endl;
   std::cout << run_ << std::endl;
   std::cout << nentriesmax_ << std::endl;

   if ( fs::exists(calibtree_path) )
   {
      if ( fs::is_directory(calibtree_path) )
      {
         for (fs::directory_entry& ls_file : fs::directory_iterator(calibtree_path)) // LOOP on RUNS!!!
         {
      
            std::string fileprefix = Form("calibTree_%d",run_);
            std::string filename = ls_file.path().filename().string();
            if ( filename.find(fileprefix) == std::string::npos || filename.find(".root") == std::string::npos ) continue;
            
            std::cout << "Working on file : " << ls_file.path().string() << std::endl;
            
            TFile * f = TFile::Open(ls_file.path().string().c_str(),"OLD");
            
            // INFO TREE - should be called once per run, one file is enough. Fill maps with information for each module
            std::string infotree_path = Form("%s/tree",infoTreePath_.c_str());
            if ( orientation_.empty() )
            {
               TTree * infotree = (TTree*) f->Get(infotree_path.c_str());
               InfoTreeBranches(infotree);
            }
            // CALIB TREE
            std::string tree_path = Form("gainCalibrationTree%s/tree",calibrationMode_.c_str());
            if ( infolocalb_->at(0) < 0.1 ) tree_path = Form("gainCalibrationTree%s0T/tree",calibrationMode_.c_str());
            TTree * tree = (TTree*) f->Get(tree_path.c_str());
            CalibTreeBranches(tree);
            
            // LOOP on EVENTS!!!
            unsigned int nentries = tree->GetEntries();
            for (unsigned int ientry = 0; ientry < nentries; ientry++)
            {
               ++count_entries;
//               if ( count_entries%100 == 0 ) std::cout << "Processed " << count_entries << "..." << std::endl;
               if ( nentriesmax_ > 0 && count_entries > nentriesmax_ )
               {
                  terminate = true;
                  break;
               }
               tree->GetEntry(ientry);
               ProcessTheEvent(LA_dict);
            } // end of events loop
            if ( terminate ) break;
                     
         } // end of file list loop
      }
   }
   
   
}

void WriteOutputs(const bool & savehistos)
{
   if ( ! savehistos ) return;
   // add the run number to the output file(s)
   if ( infolocalb_->at(0) < 0.1 )
      outputfile_ = std::regex_replace( outputfile_, std::regex(".root"), Form("_0T_%d.root",run_) );
   else
      outputfile_ = std::regex_replace( outputfile_, std::regex(".root"), Form("_4T_%d.root",run_) );

   TFile out(outputfile_.c_str(),"RECREATE");
   for ( auto h : h1_ )
   {
      if ( h.second -> GetEntries() == 0 ) continue;
      h.second -> Write();
   }
   for ( auto h : h2_ )
   {
      if ( h.second -> GetEntries() == 0 ) continue;
      h.second -> Write();
      if ( saveHistosProfile_ )
      {
         TProfile * hp = (TProfile*) h.second -> ProfileX();
         hp -> Write();
      }
   }
   if ( saveHistosMods_ )
   {
      for ( int i = 1 ; i <= nlayers_["TIB"]; ++i )
         out.mkdir(Form("modules/TIB/L%d",i));
      for ( int i = 1 ; i <= nlayers_["TOB"]; ++i )
         out.mkdir(Form("modules/TOB/L%d",i));
      for ( auto h : h2_ct_w_m_ )
      {
         if ( h.second -> GetEntries() == 0 ) continue;
         WriteOutputsModules(out,h.second);
      }
      for ( auto h : h2_t_w_m_ )
      {
         if ( h.second -> GetEntries() == 0 ) continue;
         WriteOutputsModules(out,h.second);
      }
      for ( auto h : h2_ct_var2_m_ )
      {
         if ( h.second -> GetEntries() == 0 ) continue;
         WriteOutputsModules(out,h.second);
      }
      for ( auto h : h2_t_var2_m_ )
      {
         if ( h.second -> GetEntries() == 0 ) continue;
         WriteOutputsModules(out,h.second);
      }
      for ( auto h : h2_ct_var3_m_ )
      {
         if ( h.second -> GetEntries() == 0 ) continue;
         WriteOutputsModules(out,h.second);
      }
   }
   out.Close();

}

void WriteOutputsModules(TFile & out, TH2F * h)
{
         if ( std::string(h->GetName()).find("TIB") != std::string::npos )
         {
            for ( int i = 1 ; i <= nlayers_["TIB"]; ++i )
            {
               if ( std::string(h->GetName()).find(Form("TIB_L%d",i)) != std::string::npos )
               out.cd(Form("modules/TIB/L%d",i));
            }
         }
         if ( std::string(h->GetName()).find("TOB") != std::string::npos )
         {
            for ( int i = 1 ; i <= nlayers_["TOB"]; ++i )
            {
               if ( std::string(h->GetName()).find(Form("TOB_L%d",i)) != std::string::npos )
               out.cd(Form("modules/TOB/L%d",i));
            }
         }
         h -> Write();
         if ( saveHistosProfile_ )
         {
            TProfile * hp = (TProfile*) h -> ProfileX();
            hp -> Write();
         }
   
}

void CalibTreeBranches(TTree * tree)
{
   // event data
   tree -> SetBranchAddress((eventPrefix_ + "event" + eventSuffix_).c_str(), &eventnumber_ );
   tree -> SetBranchAddress((eventPrefix_ + "run"   + eventSuffix_).c_str(), &runnumber_   );
   tree -> SetBranchAddress((eventPrefix_ + "lumi"  + eventSuffix_).c_str(), &luminumber_  );

   // calib data
   tree -> SetBranchAddress((calibPrefix_ + "trackindex" + calibSuffix_).c_str(), &trackindex_ );
   tree -> SetBranchAddress((calibPrefix_ + "rawid"      + calibSuffix_).c_str(), &rawid_      );
   tree -> SetBranchAddress((calibPrefix_ + "nstrips"    + calibSuffix_).c_str(), &nstrips_    );
   tree -> SetBranchAddress((calibPrefix_ + "localdirx"  + calibSuffix_).c_str(), &localdirx_  );
   tree -> SetBranchAddress((calibPrefix_ + "localdiry"  + calibSuffix_).c_str(), &localdiry_  );
   tree -> SetBranchAddress((calibPrefix_ + "localdirz"  + calibSuffix_).c_str(), &localdirz_  );
   tree -> SetBranchAddress((calibPrefix_ + "variance"   + calibSuffix_).c_str(), &variance_   );

   tree -> SetBranchAddress((calibPrefix_ + "localx"     + calibSuffix_).c_str(), &localx_     );
   tree -> SetBranchAddress((calibPrefix_ + "rhlocalx"   + calibSuffix_).c_str(), &rhlocalx_   );

   // track data
   tree -> SetBranchAddress((trackPrefix_ + "trackpt"        + trackSuffix_).c_str(), &trackpt_        );
   tree -> SetBranchAddress((trackPrefix_ + "tracketa"       + trackSuffix_).c_str(), &tracketa_       );
   tree -> SetBranchAddress((trackPrefix_ + "trackphi"       + trackSuffix_).c_str(), &trackphi_       );
   tree -> SetBranchAddress((trackPrefix_ + "trackhitsvalid" + trackSuffix_).c_str(), &trackhitsvalid_ );
   tree -> SetBranchAddress((trackPrefix_ + "trackchi2ndof"  + trackSuffix_).c_str(), &trackchi2ndof_  );
}

void InfoTreeBranches(TTree * tree)
{
   tree -> SetBranchAddress("rawid", &inforawid_);
   tree -> SetBranchAddress("globalZofunitlocalY", &infoglobalZofunitlocalY_);
   tree -> SetBranchAddress("localB", &infolocalb_);
   tree -> SetBranchAddress("lorentzAngle", &infola_);
   tree->GetEntry(0);
   for ( size_t im = 0; im < inforawid_->size() ; ++im )
   {
      int mod = inforawid_->at(im);
      orientation_[mod] = infoglobalZofunitlocalY_->at(im) < 0 ? -1 : 1;
      la_db_[mod] = infola_->at(im);
      // histograms for each module
      if ( saveHistosMods_ )
      {
         h1_[Form("%s_%d_nstrips"     ,ModuleLocationType(mod).c_str(),mod)]  = new TH1F (Form("%s_%d_nstrips"    ,ModuleLocationType(mod).c_str(),mod), "", 10,0,10);
         h1_[Form("%s_%d_tanthetatrk" ,ModuleLocationType(mod).c_str(),mod)]  = new TH1F (Form("%s_%d_tanthetatrk",ModuleLocationType(mod).c_str(),mod), "", 40,-1.,1.);
         h1_[Form("%s_%d_cosphitrk"   ,ModuleLocationType(mod).c_str(),mod)]  = new TH1F (Form("%s_%d_cosphitrk"  ,ModuleLocationType(mod).c_str(),mod), "", 40,-1,1);
         h1_[Form("%s_%d_variance_w2" ,ModuleLocationType(mod).c_str(),mod)]  = new TH1F (Form("%s_%d_variance_w2",ModuleLocationType(mod).c_str(),mod), "", 20,0,1);
         h1_[Form("%s_%d_variance_w3" ,ModuleLocationType(mod).c_str(),mod)]  = new TH1F (Form("%s_%d_variance_w3",ModuleLocationType(mod).c_str(),mod), "", 20,0,1);
         h2_ct_w_m_[mod]    = new TH2F (Form("ct_w_m_%s_%d"   ,ModuleLocationType(mod).c_str(),mod)   ,"", 90, -0.9, 0.9, 10, 0, 10); 
         h2_t_w_m_[mod]     = new TH2F (Form("t_w_m_%s_%d"    ,ModuleLocationType(mod).c_str(),mod)   ,"", 90, -0.9, 0.9, 10, 0, 10); 
         h2_ct_var2_m_[mod] = new TH2F (Form("ct_var2_m_%s_%d",ModuleLocationType(mod).c_str(),mod)   ,"", 90, -0.9, 0.9, 20, 0,  1); 
         h2_ct_var3_m_[mod] = new TH2F (Form("ct_var3_m_%s_%d",ModuleLocationType(mod).c_str(),mod)   ,"", 90, -0.9, 0.9, 20, 0,  1); 
         h2_t_var2_m_[mod]  = new TH2F (Form("t_var2_m_%s_%d" ,ModuleLocationType(mod).c_str(),mod)   ,"", 90, -0.9, 0.9, 20, 0,  1); 
         h2_t_var3_m_[mod]  = new TH2F (Form("t_var3_m_%s_%d" ,ModuleLocationType(mod).c_str(),mod)   ,"", 90, -0.9, 0.9, 20, 0,  1); 
      }
   }
}               


std::string ModuleLocationType(const unsigned int & mod)
{

  const SiStripDetId detid(mod);
  std::string subdet = "";
  unsigned int layer = 0;
  if ( detid.subDetector() == SiStripDetId::TIB )
  {
     subdet = "TIB";
     layer = TIBDetId(detid()).layer();
  }
  if ( detid.subDetector() == SiStripDetId::TOB )
  {
     subdet = "TOB";
     layer = TOBDetId(detid()).layer();
  }
  
  std::string type  = (detid.stereo() ? "s": "a");
  std::string d_l_t = Form("%s_L%d%s",subdet.c_str(),layer,type.c_str());
  
  if ( layer == 0 ) return subdet;
  
  return d_l_t;


}

