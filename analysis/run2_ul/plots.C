#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <iterator>
#include <stdio.h>

#include <boost/tokenizer.hpp>

TCanvas * c;
vector<string> eras;
map<string,int> era_id;
vector<double> x;
vector<double> x_err;
vector<double> la;
vector<double> la_err;
vector<double> para;
vector<double> para_err;
vector<double> parb;
vector<double> parb_err;
vector<double> pars;
vector<double> pars_err;
vector<double> chi2;
vector<double> chi2_err;

int read_data(const string & filename,const string & color);

TGraphErrors * graph(const string &, const vector<double> &, const vector<double> &,const vector<double> &,const vector<double> &);

int plots(const string & subdet, const string & color="kGreen")
{
   
   eras = {"2015A","2015B","2015C","2015D","2015E","2016A","2016B","2016C","2016D","2016E","2016F","2016G","2016H","2017A","2017B","2017C","2017D","2017E","2017F","2017G","2017H","2018A","2018B","2018C","2018D"};
   
   for ( size_t i = 0; i<eras.size(); ++i )
      era_id[eras[i]] = i+1;
      
   using namespace boost;
   
   c = new TCanvas("c","",1200,500);
   c->SetBottomMargin(0.20);
   c->SetLeftMargin(0.13);
   c->SetRightMargin(0.05);
   
   read_data(Form("%s_Run2_UL_L1MuOpen_%s.csv",subdet.c_str(),"kGreen"),"kGreen");
   
   auto x_g = x;
   
   auto g_la_g = graph("la",x,la,x_err,la_err);
   g_la_g -> SetMarkerColor(kGreen);
   g_la_g -> SetLineColor(kGreen);
   auto g_para_g = graph("para",x,para,x_err,para_err);
   g_para_g -> SetMarkerColor(kGreen);
   g_para_g -> SetLineColor(kGreen);
   auto g_parb_g = graph("parb",x,parb,x_err,parb_err);
   g_parb_g -> SetMarkerColor(kGreen);
   g_parb_g -> SetLineColor(kGreen);
   auto g_chi2_g = graph("chi2",x,chi2,x_err,chi2_err);
   g_chi2_g -> SetMarkerColor(kGreen);
   g_chi2_g -> SetLineColor(kGreen);
  
   read_data(Form("%s_Run2_UL_L1MuOpen_%s.csv",subdet.c_str(),"kBlue"),"kBlue");
   
   auto g_la_b = graph("la",x,la,x_err,la_err);
   g_la_b -> SetMarkerColor(kBlue);
   g_la_b -> SetLineColor(kBlue);
   auto g_para_b = graph("para",x,para,x_err,para_err);
   g_para_b -> SetMarkerColor(kBlue);
   g_para_b -> SetLineColor(kBlue);
   auto g_parb_b = graph("parb",x,parb,x_err,parb_err);
   g_parb_b -> SetMarkerColor(kBlue);
   g_parb_b -> SetLineColor(kBlue);
   auto g_chi2_b = graph("chi2",x,chi2,x_err,chi2_err);
   g_chi2_b -> SetMarkerColor(kBlue);
   g_chi2_b -> SetLineColor(kBlue);
   
   
   
   
   
   read_data(Form("%s_Run2_UL_L1MuOpen_%s.csv",subdet.c_str(),"kRed"),"kRed");
   
   auto g_la_r = graph("la",x,la,x_err,la_err);
   g_la_r -> SetMarkerColor(kRed);
   g_la_r -> SetLineColor(kRed);
   auto g_para_r = graph("para",x,para,x_err,para_err);
   g_para_r -> SetMarkerColor(kRed);
   g_para_r -> SetLineColor(kRed);
   auto g_parb_r = graph("parb",x,parb,x_err,parb_err);
   g_parb_r -> SetMarkerColor(kRed);
   g_parb_r -> SetLineColor(kRed);
   auto g_chi2_r = graph("chi2",x,chi2,x_err,chi2_err);
   g_chi2_r -> SetMarkerColor(kRed);
   g_chi2_r -> SetLineColor(kRed);
   
   
   //================================================
   
   TMultiGraph * mg_la = new TMultiGraph();
   
   mg_la -> Add(g_la_g);
   if ( color != "kGreen" )
   {
      mg_la -> Add(g_la_b);
      mg_la -> Add(g_la_r);
      mg_la -> Draw("apl");
   }
   else
   {
      g_la_g -> SetMarkerColor(kBlack);
      g_la_g -> SetLineColor(kBlack);
      mg_la -> Draw("ap");
      auto la_fit = g_la_g -> Fit("pol1","S");
      auto lf = g_la_g->GetListOfFunctions();
      
      TF1 * ff = (TF1*)lf->At(0);
      ff->SetLineWidth(2); 
   }
   mg_la -> SetTitle(Form("%s",subdet.c_str()));
   mg_la -> GetYaxis() -> SetRangeUser(-0.12,-0.04);
   mg_la -> GetYaxis() -> SetTitle("lorentz angle");
   mg_la -> GetXaxis() -> SetLabelSize(0.07);
   mg_la -> GetYaxis() -> SetLabelSize(0.05);
   mg_la -> GetYaxis() -> SetTitleSize(0.06);
   mg_la -> GetYaxis() -> SetTitleOffset(0.9);
   mg_la -> GetHistogram()->GetXaxis()->SetRangeUser(0,x_g.back()+2);
   
   
   TAxis* a_la = mg_la->GetXaxis();
   for ( int i = 0; i < int(x_g.back()); ++i )
   {
      int bin = a_la->FindBin(float(i+1));
      a_la -> SetBinLabel(bin,eras[i].c_str());
   }
   
   gPad->Modified();
   gPad->Update();   
   c -> SaveAs(Form("plots/%s_LA_Run2_UL_L1MuOpen_%s.png",subdet.c_str(),color.c_str()));
   
//   return 0;
   
   //================================================
   
   TMultiGraph * mg_para = new TMultiGraph();
   
   mg_para -> Add(g_para_g);
   if ( color != "kGreen" )
   {
      mg_para -> Add(g_para_b);
      mg_para -> Add(g_para_r);
      mg_para -> Draw("apl");
   }
   else
   {
      g_para_g -> SetMarkerColor(kBlack);
      g_para_g -> SetLineColor(kBlack);
      mg_para -> Draw("apl");
   }
   
   mg_para -> SetTitle(Form("%s",subdet.c_str()));
   mg_para -> GetYaxis() -> SetRangeUser(300,700);
   mg_para -> GetYaxis() -> SetTitle("parameter a");
   mg_para -> GetXaxis() -> SetLabelSize(0.07);
   mg_para -> GetYaxis() -> SetLabelSize(0.05);
   mg_para -> GetYaxis() -> SetTitleSize(0.06);
   mg_para -> GetYaxis() -> SetTitleOffset(0.9);
   mg_para -> GetHistogram()->GetXaxis()->SetRangeUser(0,x_g.back()+2);
   TAxis* a_para = mg_para->GetXaxis();
   for ( int i = 0; i < int(x_g.back()); ++i )
   {
      int bin = a_para->FindBin(float(i+1));
      a_para -> SetBinLabel(bin,eras[i].c_str());
   }
   
   gPad->Modified();
   gPad->Update();   
   c -> SaveAs(Form("plots/%s_ParA_Run2_UL_L1MuOpen_%s.png",subdet.c_str(),color.c_str()));
   
   
   //================================================
   
   TMultiGraph * mg_parb = new TMultiGraph();
   
   mg_parb -> Add(g_parb_g);
   if ( color != "kGreen" )
   {
      mg_parb -> Add(g_parb_b);
      mg_parb -> Add(g_parb_r); 
      mg_parb -> Draw("apl");
   }
   else
   {
      g_parb_g -> SetMarkerColor(kBlack);
      g_parb_g -> SetLineColor(kBlack);
      mg_parb -> Draw("apl");
   }
   mg_parb -> SetTitle(Form("%s",subdet.c_str()));
   mg_parb -> GetYaxis() -> SetRangeUser(700,1050);
   mg_parb -> GetYaxis() -> SetTitle("parameter b");
   mg_parb -> GetXaxis() -> SetLabelSize(0.07);
   mg_parb -> GetYaxis() -> SetLabelSize(0.05);
   mg_parb -> GetYaxis() -> SetTitleSize(0.06);
   mg_parb -> GetYaxis() -> SetTitleOffset(0.9);
   mg_parb -> GetHistogram()->GetXaxis()->SetRangeUser(0,x_g.back()+2);
   TAxis* a_parb = mg_parb->GetXaxis();
   for ( int i = 0; i < int(x_g.back()); ++i )
   {
      int bin = a_parb->FindBin(float(i+1));
      a_parb -> SetBinLabel(bin,eras[i].c_str());
   }
   
   gPad->Modified();
   gPad->Update();   
   c -> SaveAs(Form("plots/%s_ParB_Run2_UL_L1MuOpen_%s.png",subdet.c_str(),color.c_str()));
   
   
   //================================================
   
   TMultiGraph * mg_chi2 = new TMultiGraph();
   
   mg_chi2 -> Add(g_chi2_g);
   if ( color != "kGreen" )
   {
      mg_chi2 -> Add(g_chi2_b);
      mg_chi2 -> Add(g_chi2_r);
      mg_chi2 -> Draw("apl");
   }
   else
   {
      g_chi2_g -> SetMarkerColor(kBlack);
      g_chi2_g -> SetLineColor(kBlack);
      mg_chi2 -> Draw("apl");
   }
   
   mg_chi2 -> SetTitle(Form("%s",subdet.c_str()));
   mg_chi2 -> GetYaxis() -> SetRangeUser(0,3);
   mg_chi2 -> GetYaxis() -> SetTitle("chi2/ndf");
   mg_chi2 -> GetXaxis() -> SetLabelSize(0.07);
   mg_chi2 -> GetYaxis() -> SetLabelSize(0.05);
   mg_chi2 -> GetYaxis() -> SetTitleSize(0.06);
   mg_chi2 -> GetYaxis() -> SetTitleOffset(0.9);
   mg_chi2 -> GetHistogram()->GetXaxis()->SetRangeUser(0,x_g.back()+2);
   TAxis* a_chi2 = mg_chi2->GetXaxis();
   for ( int i = 0; i < int(x_g.back()); ++i )
   {
      int bin = a_chi2->FindBin(float(i+1));
      a_chi2 -> SetBinLabel(bin,eras[i].c_str());
   }
   
   gPad->Modified();
   gPad->Update();   
   c -> SaveAs(Form("plots/%s_Chi2_Run2_UL_L1MuOpen_%s.png",subdet.c_str(),color.c_str()));
   
   
   return 0;
   
}

TGraphErrors * graph(const string & var, const vector<double> &x , const vector<double> & y,const vector<double> & xe,const vector<double> &ye)
{
   
   int n = int(x.size());
   
   TGraphErrors * g = new TGraphErrors(n,&x[0],&y[0],&xe[0],&ye[0]);
   g -> SetMarkerStyle(20);
   g -> SetMarkerSize(1.2);
   g -> SetLineWidth(2);
   g -> GetXaxis() -> SetNdivisions(int(x.back()));
   TAxis* a = g->GetXaxis();
   for ( int i = 0; i < int(x.back()); ++i )
   {
      int bin = a->FindBin(float(i+1));
      a -> SetBinLabel(bin,eras[i].c_str());
   }
//   g -> Draw("APL");

   return g;
}


int read_data(const string & filename, const string & color)
{
   

   using namespace boost;
   
   x.clear();
   x_err.clear();
   la.clear();
   la_err.clear();
   para.clear();
   para_err.clear();
   parb.clear();
   parb_err.clear();
   pars.clear();
   pars_err.clear();
   chi2.clear();
   chi2_err.clear();
   
   string data(filename);
   
   ifstream in(data.c_str());
   
   if (!in.is_open()) return 1;
   
   typedef tokenizer< escaped_list_separator<char> > Tokenizer;
   
   vector< string > vec;
   string line;
   
   while (getline(in,line))
   {
      Tokenizer tok(line);
      vec.assign(tok.begin(),tok.end());
      string era    = vec[0];
      double xi = double(era_id[era]);
      if ( color == "kBlue" ) xi -= 0.1;
      if ( color == "kRed"  ) xi += 0.1;
      x.push_back(xi);
      x_err.push_back(0.);
      la.push_back(stod(vec[1])); 
      la_err.push_back(stod(vec[2])); 
      para.push_back(stod(vec[3])); 
      para_err.push_back(stod(vec[4])); 
      parb.push_back(stod(vec[5])); 
      parb_err.push_back(stod(vec[6])); 
      pars.push_back(stod(vec[7])); 
      pars_err.push_back(stod(vec[8])); 
      chi2.push_back(stod(vec[9])); 
      chi2_err.push_back(0.); 
   }
   return 0;
}
