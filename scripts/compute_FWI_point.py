# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 11:51:17 2019

@author: guillaume
"""
from fwi_compute_index import FWICLASS

def main():
    ffmc0 = 85.0
    dmc0 = 6.0
    dc0 = 15.0
    infile = open('data.txt','r')
    lines = infile.readlines()[1:]
    outfile = open('fwioutput.txt','w')
    try:
        for line in lines:
            mth,day,temp,rhum,wind,prcp=[float(field) for field in line.strip().split()]
            if rhum>100.0:
                rhum = 100.0
            mth = int(mth)
            fwisystem = FWICLASS(temp,rhum,wind,prcp)
            ffmc = fwisystem.FFMCcalc(ffmc0)
            dmc = fwisystem.DMCcalc(dmc0,mth)
            dc = fwisystem.DCcalc(dc0,mth)
            isi = fwisystem.ISIcalc(ffmc)
            bui = fwisystem.BUIcalc(dmc,dc)
            fwi = fwisystem.FWIcalc(isi,bui)
            ffmc0 = ffmc
            dmc0 = dmc
            dc0 = dc
            print(fwi)
            outfile.write('%s %s %s %s %s %s\n'%(str(ffmc),str(dmc),str(dc),str(isi),str(bui),str(fwi)))
    finally:
        infile.close()
        outfile.close()
main()
