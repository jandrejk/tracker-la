import ROOT as R
import glob
import re
import time
R.gROOT.SetBatch(True)

file_names =  glob.glob("SiStripLAMonitor_*.root")

PEAK_modes = ["318050", "318051", "318112", "318116", "318118", "318517", "318519", "318520", "318521"]
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
orientations = [
    "vplus",
    "vminus",
]

N_orienations = len(orientations)



def ProcessFile(name,mode,B) :
    f_run_number = re.findall(r'\d{6}',name)[0]
    f_Bfield = re.findall(r'\d{1}',name)[0]
    
    f_mode = "DECO"
    if (f_run_number in PEAK_modes) :
        f_mode = "PEAK"
    
    if (mode != f_mode or Bfield != f_Bfield) :
        return False
    else :
        return True


def PreProcess (mode, Bfield) :
    output = R.TFile("{0}_{1}T.root".format(mode,Bfield), "RECREATE")

    h_out_localx = [R.TH1F("{layer}_{variable}_{orientation}".format(layer=l,variable="localx",orientation=o),"",200,-6,6) for l in layers for o in orientations]
    h_out_rhlocalx = [R.TH1F("{layer}_{variable}_{orientation}".format(layer=l,variable="rhlocalx",orientation=o),"",200,-6,6) for l in layers for o in orientations]
    h_out_shift = [R.TH1F("{layer}_{variable}_{orientation}".format(layer=l,variable="shift",orientation=o),"", 200,-0.05,0.05) for l in layers for o in orientations]
    h_out_shift_tanthetatrk_minus_tanthetaLA = [R.TH2F("{layer}_{variable}_{orientation}".format(layer=l,variable="shift_tanthetatrk_minus_tanthetaLA",orientation=o),"", 400, -0.001, 0.001, 360, -0.9, 0.9) for l in layers for o in orientations]

    for f in file_names : 

        if (ProcessFile(name=f,mode=mode,B=Bfield) == False) :
            continue
        print "file {f} will be processed".format(f=f)

        Tfile = R.TFile("{}".format(f))
        keylist = Tfile.GetListOfKeys()
        for il, l in enumerate(layers) :
            for io, o in enumerate(orientations) :

                if keylist.Contains("{0}_localx_{1}".format(l,o)) == False :
                    print "no hits in layer: {0} with orientation: {1}".format(l,o)
                    continue

                localx   = Tfile.Get("{0}_localx_{1}".format(l,o))
                rhlocalx   = Tfile.Get("{0}_rhlocalx_{1}".format(l,o))
                shift   = Tfile.Get("{0}_shift_{1}".format(l,o))
                shift_tanthetatrk_minus_tanthetaLA   = Tfile.Get("{0}_shift_tanthetatrk_minus_tanthetaLA_{1}".format(l,o))
                
                h_out_localx[N_orienations*il+io].Add(localx)
                h_out_rhlocalx[N_orienations*il+io].Add(rhlocalx)
                h_out_shift[N_orienations*il+io].Add(shift)
                h_out_shift_tanthetatrk_minus_tanthetaLA[N_orienations*il+io].Add(shift_tanthetatrk_minus_tanthetaLA)
                
        Tfile.Close()

    output.cd()
    for h1 in [h_out_localx, h_out_rhlocalx, h_out_shift, h_out_shift_tanthetatrk_minus_tanthetaLA] :
        for h2 in h1 :
            h2.Write()
    output.Close()


def Plot() :
    Tfile_PEAK = R.TFile("PEAK_0T.root")
    Tfile_DECO = R.TFile("DECO_0T.root")
    #shift = DECO - PEAK
    for l in layers :
        for o in orientations :

            h2d_DECO = Tfile_DECO.Get("{0}_shift_tanthetatrk_minus_tanthetaLA_{1}".format(l,o))
            h2d_PEAK = Tfile_PEAK.Get("{0}_shift_tanthetatrk_minus_tanthetaLA_{1}".format(l,o))

            
            h1_DECO = h2d_DECO.ProfileY().Rebin(36)
            h1_DECO.SetDirectory(0)
            h1_PEAK = h2d_PEAK.ProfileY().Rebin(36)
            h1_PEAK.SetDirectory(0)



            h1_PEAK.Scale(10**6)
            h1_DECO.Scale(10**6)
            
            fit_function_DECO = R.TF1("lin_fit_deco","pol1",-0.9,0.9)
            fit_function_PEAK = R.TF1("lin_fit_peak","pol1",-0.9,0.9)
            
            h1_DECO.Fit(fit_function_DECO,"0")
            h1_PEAK.Fit(fit_function_PEAK,"0")
            
                
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
            barrelID, layer = l.split("_")
            module_type = ""
            module_orientation = ""
            if (layer[2] == 'a') :
                module_type = "analog"
            elif (layer[2] == "s") :
                module_type = "stereo"

            if (o == 'vplus') :
                module_orientation = "parallel to z-axis"
            elif (o == "vminus") :
                module_orientation = "anti-parallel to z-axis"
            
            Text = "{} - layer: {} - module type: {} - module orientation: {}".format(barrelID, layer[1], module_type, module_orientation)            
            # text.Add("{}")
            paves.AddText(Text)
            paves.Draw("same")


            leg = R.TLegend(.1,.7,.5,.9)
            leg.SetTextSize(0.04)
            leg.SetFillColor(0)
            leg.AddEntry(h1_PEAK,"PEAK: #Deltaw = {:.3} #pm {:.3}".format(fit_function_PEAK.GetParameter(1),fit_function_PEAK.GetParError(1)),"l")
            leg.AddEntry(h1_DECO,"DECO: #Deltaw = {:.3} #pm {:.3}".format(fit_function_DECO.GetParameter(1),fit_function_DECO.GetParError(1)),"l")
            leg.Draw("same")

            c.SaveAs("{}/{}_L{}MT_{}_O{}.png".format("BP_plots",barrelID, layer[1], module_type, o) )

            cDiff = R.TCanvas()
            cDiff.cd()

            h_DECOminusPEAK = h1_DECO.Clone("")
            h_DECOminusPEAK.Add(h1_PEAK,-1)
            print h_DECOminusPEAK.GetNbinsX()
            fit_function_DECOminusPEAK = R.TF1("lin_fit_decominuspeak","pol1",-0.9,0.9)
            print fit_function_DECOminusPEAK
            h_DECOminusPEAK.Fit(fit_function_DECOminusPEAK,"0")

            
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
            
            
            Text = "{} - layer: {} - module type: {} - module orientation: {}".format(barrelID, layer[1], module_type, module_orientation)            
            paves.AddText(Text)
            paves.Draw("same")


            leg = R.TLegend(.1,.7,.6,.9)
            leg.SetTextSize(0.04)
            leg.SetFillColor(0)
            leg.AddEntry(h_DECOminusPEAK,"DECO-PEAK: #Deltaw = {:.3} #pm {:.3}".format(fit_function_DECOminusPEAK.GetParameter(1),fit_function_DECOminusPEAK.GetParError(1)),"l")
            leg.Draw("same")

            cDiff.SaveAs("{}/{}_L{}MT_{}_O{}_diff.png".format("BP_plots",barrelID, layer[1], module_type, o) )

            txtoutput = open("out.txt","a")
            txtoutput.write("{}_{}: {} +/- {} \n".format(l,o,fit_function_DECOminusPEAK.GetParameter(1),fit_function_DECOminusPEAK.GetParError(1)))
            txtoutput.close()


           
    Tfile_PEAK.Close()
    Tfile_DECO.Close()

s = time.time()
# for mode in ["PEAK", "DECO"] :
#     for Bfield in ["0","4"] :
#         PreProcess(mode,Bfield)

Plot()
print "execution time: {}".format(time.time()-s)
