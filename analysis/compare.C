int compare(const string & det = "TIB")
{
   gStyle->SetOptTitle(0);
   TFile * f1 = new TFile("run2_ul/2018/Cosmics2018D/HLT_L1SingleMuCosmics_ptmin5_hitmin8_chi2max1p4/SiStripLAMonitor_4T.root","old");
   TFile * f2 = new TFile("run2_old/2018/Cosmics2018D/HLT_L1SingleMuCosmics_ptmin5_hitmin8_chi2max1p4/SiStripLAMonitor_4T.root","old");
   
   map<string,vector<string> > layers;
   
//   vector<string> layers{"TIB_L1s","TIB_L1a","TIB_L2s","TIB_L2a","TIB_L3a","TIB_L4a","TOB_L1s","TOB_L1a","TOB_L2s","TOB_L2a","TOB_L3a","TOB_L4a","TOB_L5a","TOB_L6a"};
   layers["TIB"] = {"TIB_L1s","TIB_L1a","TIB_L2s","TIB_L2a","TIB_L3a","TIB_L4a"};
   layers["TOB"] = {"TOB_L1s","TOB_L1a","TOB_L2s","TOB_L2a","TOB_L3a","TOB_L4a","TOB_L5a","TOB_L6a"};
   
   map<string, float> ymin;
   map<string, float> ymax;
         
   ymin["TIB"] = 2.10;
   ymax["TIB"] = 3.45;
   
   ymin["TOB"] = 2.35;
   ymax["TOB"] = 3.70;
   
   map<string, TProfile*> h1;
   map<string, TProfile*> h2;
   
   TCanvas * c = new TCanvas("c","",1920,870);
   c->Divide(4,2);
   
   for ( auto & l : layers[det] )
   {
      auto i = &l - &layers[det][0];
      h1[l] = (TProfile*) f1->Get(Form("%s_thetatrk_nstrip_pfx",l.c_str()));
      h1[l] -> SetName(Form("h1_%s_thetatrk_nstrip",l.c_str()));
      h1[l] -> SetTitle("2018D Legacy");
      h2[l] = (TProfile*) f2->Get(Form("%s_thetatrk_nstrip_pfx",l.c_str()));
      h2[l] -> SetName(Form("h2_%s_thetatrk_nstrip",l.c_str()));
      h2[l] -> SetTitle("2018D Private");
      h1[l] -> Rebin(5);
      h2[l] -> Rebin(5);
      
      h1[l] -> SetStats(0);
      h1[l] -> GetXaxis() -> SetRangeUser(-0.4,0.4);
      h1[l] -> GetYaxis() -> SetRangeUser(ymin[det],ymax[det]);
//      h1[l] -> GetXaxis() -> SetNdivisions(505);
//      h1[l] -> GetYaxis() -> SetNdivisions(505);
      h1[l] -> GetXaxis() -> SetTitle("#theta_{T}");
      h1[l] -> GetYaxis() -> SetTitle("<strips>");
      h1[l] -> GetYaxis() -> SetTitleSize(20);
      h1[l] -> GetYaxis() -> SetTitleFont(43);
      h1[l] -> GetYaxis() -> SetTitleOffset(2.5);
      h1[l] -> GetYaxis() -> SetLabelFont(43);
      h1[l] -> GetYaxis() -> SetLabelSize(20);
      h1[l] -> GetXaxis() -> SetTitleFont(43);
      h1[l] -> GetXaxis() -> SetTitleOffset(2.5);
      h1[l] -> GetXaxis() -> SetTitleSize(20);
      h1[l] -> GetXaxis() -> SetLabelFont(43);
      h1[l] -> GetXaxis() -> SetLabelSize(20);
      h1[l] -> SetMarkerStyle(20);
      h1[l] -> SetMarkerColor(1);
      h1[l] -> SetLineColor(1);
      h1[l] -> SetLineWidth(2);
      h2[l] -> SetLineColor(2);
      h2[l] -> SetLineWidth(2);
      
      c->cd(i+1);
      c->GetPad(i+1)->SetBottomMargin(0.12); 
      c->GetPad(i+1)->SetLeftMargin(0.12);
      c->GetPad(i+1)->SetTopMargin(0.08);
      
      h1[l] -> Draw("error same");
      h2[l] -> Draw("hist same");
      
      c->GetPad(i+1)-> BuildLegend(0.15,0.7,0.7,0.9,Form("%s: HLT_L1SingleMuCosmics",l.c_str()));
      
   }
   c -> SaveAs(Form("%s_MuCosmics.png",det.c_str()));
   
   return 0;
}
