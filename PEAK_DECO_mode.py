from CalibTrees import CalibTrees
from payloadInspector import dump
import csv
import os


for CalibTree in CalibTrees :
    RunIDmodes = {}
    
    outputFileName = CalibTree["tag"].replace("/","_")
    
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
    
    w = csv.writer(open("{}.csv".format(outputFileName), "w"))
    for key, val in dict.items():
        w.writerow([key, val])