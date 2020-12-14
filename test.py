import ROOT as R
import glob
import re
import time
import os
R.gROOT.SetBatch(True)
from CalibTrees import CalibTrees
from payloadInspector import dump
from datetime import date
import argparse

argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--year',           		action='store',      default='2015_MuOpen', choices=["2015_MuOpen"], help="Choose cosmic run")
args = argParser.parse_args()


s = time.time()

# Collect files in PEAK and DECO mode with and without magnetic field
for mode in ["PEAK", "DECO"] :
    for Bfield in ["0","4"] :
        PreProcess(mode,Bfield)



RunIDmodes = {}
for CalibTree in CalibTrees[:1] :
    print("looping only over 2015 calib trees!")
    files = os.listdir(CalibTree["path"])
    runIDs = list(dict.fromkeys([f.split("_")[1] for f in files]))
    print "{} run IDs found".format(len(runIDs))
    for r in runIDs :
        m_ID = min(dump["plot_data"], key=lambda x:( int(r)-x["x"] if x["x"]<int(r) else 10**6 ))["y"]
        
        if m_ID == 37 :
            RunIDmodes[r] = "DECO"
        elif m_ID == 47 :
            RunIDmodes[r] = "PEAK"
        else :
            print("no mode ?!? exiting...")
            exit(0)

print date.today()        

exit(0)

directory=CalibTrees[0]["tag"]
stamp=date.today()

bin_angle = 0.5
bin_angle2 = 0.5

bin_shift = 0.01
bin_clusterwidth = 1.0
file_names =  glob.glob("{}/SiStripLAMonitor_*.root".format(directory))

layers = [
        "TIB_L1a",
        "TIB_L1s",
        "TIB_L2a",
        "TIB_L2s",
        "TIB_L3a",
        "TIB_L4a",
        "TOB_L1a",  
        "TOB_L1s",  
        "TOB_L2a",  
        "TOB_L2s",  
        "TOB_L3a",  
        "TOB_L4a",  
        "TOB_L5a",  
        "TOB_L6a",  
    ]
incl_layers = [
    "TIB_L1",
    "TIB_L2",
    "TOB_L1",    
    "TOB_L2",    
    "TIB",
    "TOB",  
]
orientations = [
    "vplus",
    "vminus",
]

N_orienations = len(orientations)

def AddInclusivePlots (mode="PEAK", Bfield = "0",direction="shift") :
    Tfile = R.TFile("{0}/{1}_{2}T.root".format(directory,mode,Bfield),"UPDATE")
    
    for o in orientations:
        for t in ["TIB","TOB"] :
            for l in ["L1","L2"] :
                h2d_s = Tfile.Get("{0}_{1}s_{2}_tanthetatrk_minus_tanthetaLA_{3}".format(t,l,direction,o))
                h2d_a = Tfile.Get("{0}_{1}a_{2}_tanthetatrk_minus_tanthetaLA_{3}".format(t,l,direction,o))
                h2d_s.SetDirectory(0)
                h2d_a.SetDirectory(0)

                h2d = h2d_s.Clone("{0}_{1}_{2}_tanthetatrk_minus_tanthetaLA_{3}".format(t,l,direction,o))
                h2d.Add(h2d_a,1)
                h2d.Write()
    
    TB_L_dict = {
        "TIB" : ["L1a", "L1s", "L2a", "L2s", "L3a", "L4a"],
        "TOB" : ["L1a", "L1s", "L2a", "L2s", "L3a", "L4a","L5a","L6a"],   
    }

    

    for o in orientations:
        for barrel in TB_L_dict.keys() :
            for il,l in enumerate(TB_L_dict[barrel]) :
               
                h2d_in = Tfile.Get("{0}_{1}_{2}_tanthetatrk_minus_tanthetaLA_{3}".format(barrel,l,direction,o))
                h2d_in.SetDirectory(0)
                h_shift = Tfile.Get("{0}_{1}_shift_{2}".format(barrel,l,o))
                h_shift.SetDirectory(0)
                if (il == 0) :
                    h2d_out     = h2d_in.Clone("{0}_{1}_tanthetatrk_minus_tanthetaLA_{2}".format(barrel,direction,o))
                    h_shift_out = h_shift.Clone("{0}_shift_{1}".format(barrel,o))
                else :
                    h2d_out.Add(h2d_in)
                    h_shift_out.Add(h_shift)
            h2d_out.Write()
            h_shift_out.Write()


    Tfile.Close()
    

def ProcessFile(name,mode,B) :
    f_run_number = re.search(r'T_*[0-9]{6}.root',name).group()[2:8]
    f_Bfield = re.search(r'[0-9]T',name).group()[0]

    
    f_mode = RunIDmodes[f_run_number]
    
    # if (f_run_number in PEAK_modes) :
    #     f_mode = "PEAK"
    
    if (mode != f_mode or Bfield != f_Bfield) :
        return False
    else :
        return True


def PreProcess (mode, Bfield,directory,orientations,layers) :
    
    output = R.TFile("{0}/{1}_{2}T.root".format(directory,mode,Bfield), "RECREATE")

    h_out_localx = [R.TH1F("{layer}_{variable}_{orientation}".format(layer=l,variable="localx",orientation=o),"",200,-6,6) for l in layers for o in orientations]
    h_out_rhlocalx = [R.TH1F("{layer}_{variable}_{orientation}".format(layer=l,variable="rhlocalx",orientation=o),"",200,-6,6) for l in layers for o in orientations]
    h_out_shift = [R.TH1F("{layer}_{variable}_{orientation}".format(layer=l,variable="shift",orientation=o),"", 200,-0.05,0.05) for l in layers for o in orientations]
    h_out_shift_tanthetatrk_minus_tanthetaLA = [R.TH2F("{layer}_{variable}_{orientation}".format(layer=l,variable="shift_tanthetatrk_minus_tanthetaLA",orientation=o),"", 12, -bin_angle, bin_angle, 400, -bin_shift, bin_shift) for l in layers for o in orientations]
    h_out_clusterwidth_tanthetatrk_minus_tanthetaLA = [R.TH2F("{layer}_{variable}_{orientation}".format(layer=l,variable="clusterwidth_tanthetatrk_minus_tanthetaLA",orientation=o),"", 12, -bin_angle, bin_angle, 11, 0, 10) for l in layers for o in orientations]

    for f in file_names : 
        if (ProcessFile(name=f,mode=mode,B=Bfield) == False) :
            continue
        print "file {f} will be processed".format(f=f)
        
        
        Tfile = R.TFile("{}".format(f))
        keylist = Tfile.GetListOfKeys()
        for il, l in enumerate(layers) :
            for io, o in enumerate(orientations) :

                if keylist.Contains("{0}_localx_{1}".format(l,o)) == False :
                    # print "no hits in layer: {0} with orientation: {1}".format(l,o)
                    continue

                localx   = Tfile.Get("{0}_localx_{1}".format(l,o))
                rhlocalx   = Tfile.Get("{0}_rhlocalx_{1}".format(l,o))
                shift   = Tfile.Get("{0}_shift_{1}".format(l,o))
                shift_tanthetatrk_minus_tanthetaLA   = Tfile.Get("{0}_shift_tanthetatrk_minus_tanthetaLA_{1}".format(l,o))
                clusterwidth_tanthetatrk_minus_tanthetaLA   = Tfile.Get("{0}_clusterwidth_tanthetatrk_minus_tanthetaLA_{1}".format(l,o))
                
                h_out_localx[N_orienations*il+io].Add(localx)
                h_out_rhlocalx[N_orienations*il+io].Add(rhlocalx)
                h_out_shift[N_orienations*il+io].Add(shift)
                h_out_shift_tanthetatrk_minus_tanthetaLA[N_orienations*il+io].Add(shift_tanthetatrk_minus_tanthetaLA)
                h_out_clusterwidth_tanthetatrk_minus_tanthetaLA[N_orienations*il+io].Add(clusterwidth_tanthetatrk_minus_tanthetaLA)
                
        Tfile.Close()

    output.cd()
    for h1 in [h_out_localx, h_out_rhlocalx, h_out_shift, h_out_shift_tanthetatrk_minus_tanthetaLA, h_out_clusterwidth_tanthetatrk_minus_tanthetaLA] :
        for h2 in h1 :
            h2.Write()
    output.Close()
    AddInclusivePlots(mode=mode,Bfield=Bfield)
    # AddInclusivePlots(mode=mode,Bfield=Bfield,direction="clusterwidth")



def Plot(directory,Bfield = "0",) :
    if not os.path.exists("{}/BP_plots".format(directory)):
        os.makedirs("{}/BP_plots".format(directory))

    Tfile_PEAK = R.TFile("{0}/PEAK_{1}T.root".format(directory,Bfield))
    Tfile_DECO = R.TFile("{0}/DECO_{1}T.root".format(directory,Bfield))
    #shift = DECO - PEAK
    for o in orientations :
        for l in layers+incl_layers :
        
            h2d_DECO = Tfile_DECO.Get("{0}_shift_tanthetatrk_minus_tanthetaLA_{1}".format(l,o))
            h2d_PEAK = Tfile_PEAK.Get("{0}_shift_tanthetatrk_minus_tanthetaLA_{1}".format(l,o))
            
            h1_DECO = h2d_DECO.RebinX(1).ProfileX()
            h1_DECO.SetDirectory(0)
            h1_PEAK = h2d_PEAK.RebinX(1).ProfileX()
            h1_PEAK.SetDirectory(0)
            # h1_DECO.Sumw2()
            # h1_PEAK.Sumw2()


            h1_PEAK.Scale(10**3)
            h1_DECO.Scale(10**3)

            
            fit_function_DECO = R.TF1("lin_fit_deco","pol1",-bin_angle2,bin_angle2)
            fit_function_PEAK = R.TF1("lin_fit_peak","pol1",-bin_angle2,bin_angle2)
            
            h1_DECO.Fit(fit_function_DECO,"Q0")
            h1_PEAK.Fit(fit_function_PEAK,"Q0")
            
                
            c = R.TCanvas()
            R.gStyle.SetOptStat(0)
            h1_PEAK.GetYaxis().SetTitle("<#Deltau_{mode}> [#mum]")
            h1_PEAK.GetXaxis().SetTitle("<tan(#theta_{trk}) - tan(#theta_{LA})>")
            h1_PEAK.SetLineColor(2)
            h1_PEAK.Draw("E1")
            
            fit_function_PEAK.SetLineColor(2)
            fit_function_PEAK.SetLineStyle(2)
            fit_function_PEAK.Draw("same")

            h1_DECO.SetLineColor(4)
            h1_DECO.Draw("same E1")
            fit_function_DECO.SetLineColor(4)
            fit_function_DECO.SetLineStyle(2)
            fit_function_DECO.Draw("same")

            
            paves = R.TPaveText(0.1,0.9,0.5,1.0,"NDC")
            paves.SetFillColor(0) # text is black on white
            paves.SetTextSize(0.03) 
            paves.SetTextAlign(12)
            
            if ("_" in l) :
                barrelID, layer = l.split("_")
                layer_ID = layer[1]
            else :
                barrelID = l
                layer = "incl"
                layer_ID = "incl"

            module_type = "incl"
            module_orientation = ""
            if (len(layer) == 3) :
                if (layer[2] == 'a') :
                    module_type = "analog"
                elif (layer[2] == "s") :
                    module_type = "stereo"

            if (o == 'vplus') :
                module_orientation = "parallel to z-axis"
            elif (o == "vminus") :
                module_orientation = "anti-parallel to z-axis"
            
            Text = "B-field: {} T - barrel type: {} - layer: {} - module type: {} - module orientation: {}".format(Bfield, barrelID, layer_ID, module_type, module_orientation)            

            paves.AddText(Text)
            paves.Draw("same")


            leg = R.TLegend(.1,.7,.5,.9)
            leg.SetTextSize(0.04)
            leg.SetFillColor(0)
            leg.AddEntry(h1_PEAK,"PEAK: #Deltaw = {:.3} #pm {:.3}".format(fit_function_PEAK.GetParameter(1),fit_function_PEAK.GetParError(1)),"l")
            leg.AddEntry(h1_DECO,"DECO: #Deltaw = {:.3} #pm {:.3}".format(fit_function_DECO.GetParameter(1),fit_function_DECO.GetParError(1)),"l")
            leg.Draw("same")

            c.SaveAs("{}/BP_plots/{}T_{}_L{}MT_{}_O{}.png".format(directory,Bfield, barrelID, layer_ID, module_type, o) )

            cDiff = R.TCanvas()
            cDiff.cd()



            h_DECOminusPEAK = R.TH1F("{0}_shift_tanthetatrk_minus_tanthetaLA_{1}_DECO_minus_PEAK".format(l,o),"",h1_DECO.GetNbinsX(),-bin_angle,bin_angle)
            for b in xrange(h_DECOminusPEAK.GetNbinsX()) :
                
                h_DECOminusPEAK.SetBinContent(b+1,h1_DECO.GetBinContent(b+1)-h1_PEAK.GetBinContent(b+1))
                err = R.TMath.Sqrt(h1_DECO.GetBinError(b+1)**2 + h1_PEAK.GetBinError(b+1)**2)
                h_DECOminusPEAK.SetBinError(b+1,err)
            
            #     print "DECO bin content: {}".format(h1_DECO.GetBinContent(b+1))
            #     print "PEAK bin content: {}".format(h1_PEAK.GetBinContent(b+1))
            #     print "error: {}".format(err)
            #     print "-"*90
            
            # exit(0)
            # h_DECOminusPEAK.Sumw2()

            
            # h_DECOminusPEAK.SetDirectory(0)
            # h_DECOminusPEAK.Sumw2()
            # h_DECOminusPEAK.Add(h1_PEAK,1)

            fit_function_DECOminusPEAK = R.TF1("lin_fit_decominuspeak","pol1",-bin_angle2,bin_angle2)
            h_DECOminusPEAK.Fit(fit_function_DECOminusPEAK,"Q0")

            
            R.gStyle.SetOptStat(0)
            h_DECOminusPEAK.GetYaxis().SetTitle("<#Deltau_{DECO} - #Deltau_{PEAK} > [#mum]")
            h_DECOminusPEAK.GetXaxis().SetTitle("<tan(#theta_{trk}) - tan(#theta_{LA})>")
            h_DECOminusPEAK.SetLineColor(6)
            h_DECOminusPEAK.Draw("E1")
            
            fit_function_DECOminusPEAK.SetLineColor(6)
            fit_function_DECOminusPEAK.SetLineStyle(2)
            fit_function_DECOminusPEAK.Draw("same")

            
            paves = R.TPaveText(0.1,0.9,0.5,1.0,"NDC")
            paves.SetFillColor(0) # text is black on white
            paves.SetTextSize(0.03) 
            paves.SetTextAlign(12)
            
            
            Text = "B-field: {} T - barrel type: {} - layer: {} - module type: {} - module orientation: {}".format(Bfield, barrelID, layer_ID, module_type, module_orientation)            
            paves.AddText(Text)
            paves.Draw("same")


            leg = R.TLegend(.1,.7,.7,.9)
            leg.SetTextSize(0.04)
            leg.SetFillColor(0)
            leg.AddEntry(h_DECOminusPEAK,"DECO-PEAK: #Deltaw = {:.3} #pm {:.3}".format(fit_function_DECOminusPEAK.GetParameter(1),fit_function_DECOminusPEAK.GetParError(1)),"l")
            leg.Draw("same")

            cDiff.SaveAs("{}/BP_plots/{}T_{}_L{}MT_{}_O{}_diff.png".format(directory,Bfield, barrelID, layer_ID, module_type, o) )

            txtoutput = open("{}/BP_plots/{}T_output.txt".format(directory,Bfield),"a")
            txtoutput.write("{}_{}: \t {} \t +/- \t {} \n".format(l,o,fit_function_DECOminusPEAK.GetParameter(1),fit_function_DECOminusPEAK.GetParError(1)))
            txtoutput.close()
    
           
    Tfile_PEAK.Close()
    Tfile_DECO.Close()








#---------------
s = time.time()
for mode in ["PEAK", "DECO"] :
    for Bfield in ["0","4"] :
        PreProcess(mode,Bfield)



for Bfield in ["0","4"] :
    Plot(directory=directory,Bfield=Bfield)
print "execution time: {}".format(time.time()-s)

peak = R.TFile("{}/PEAK_4T.root".format(directory))
deco = R.TFile("{}/DECO_4T.root".format(directory))

for i,f in enumerate([peak, deco]) :
    for b in ["TIB", "TOB"] :
        h_shift_minus = f.Get("{}_shift_vminus".format(b))
        h_shift_plus  = f.Get("{}_shift_vplus".format(b))


        out = R.TCanvas()
        h_shift_minus.Draw()
        h_shift_minus.GetXaxis().SetTitle("#Deltau [cm]")
        
        
        h_shift_plus.SetLineColor(2)
        h_shift_plus.Draw("same")

        leg = R.TLegend(.1,.7,.45,.9)
        leg.SetTextSize(0.03)
        leg.SetFillColor(0)
        leg.AddEntry(h_shift_minus,"anti-parallel to z-axis","l")
        leg.AddEntry(h_shift_plus,"parallel to z-axis","l")
        leg.Draw("same")

        if (i==0) :
            mode = "PEAK"
        else :
            mode = "DECO"
        h_shift_minus.SetTitle("4T - {} - {} mode".format(b,mode))
        out.SaveAs("{}/BP_plots/shift_{}_{}.png".format(directory,mode,b) )
peak.Close()
deco.Close()
