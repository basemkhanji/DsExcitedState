#!/usr/bin/env python
import ROOT
from ROOT import *
from math import *
import sys
import os
import threading ,subprocess, multiprocessing
import multiprocessing as mp
import time
from os.path import isfile
#rom termcolor import colored
###########################################################################
def PrepareTrees( mychain , Mycuts ):
    mychain_new = mychain.CopyTree( Mycuts_t , "" , 10000000000 , 0 )
    return mychain_new
###########################################################################
def EffError(a,b):
    r ,  sr = 1 , 1
    if( b!=0):
        r = a/b
        sr = sqrt(r*(1-r)/b)
    return sr
###########################################################################
def Colorprint(mytext , color):
    if (color in 'red' )   :print("\033[2;31;107m " + mytext + " \033[0m" )
    if (color in 'blue' )  :print("\033[2;34;107m " + mytext + " \033[0m" )
    if (color in 'green') :print("\033[2;32;107m " + mytext + " \033[0m" )
###########################################################################
def GetFilesinDir( mypath  ):
    list_ofnames = []
    for root, dirs, files in os.walk( mypath ):#os.walk("./"):
        for file in files:
            if file.endswith(".root"):
                if file not in list_ofnames:
                    print '============================='
                    print 'new file name is found', file , 'Appending ... '
                    print '============================='
                    list_ofnames.append(  file)
                    
    return list_ofnames
###########################################################################
## KK_M :
K1_PE = "sqrt(493.67**2 + (K1_PX**2 + K1_PY**2 + K1_PZ**2))"  ;
K2_PE = "sqrt(493.67**2 + (K2_PX**2 + K2_PY**2 + K2_PZ**2))" ;
K1_PE_K2_PE = K1_PE+"*"+K2_PE ;

## KKpi_M under Kpipi_M hypothesis :
K2asPi_PE = "sqrt(139.57**2 + (K2_PX**2 + K2_PY**2 + K2_PZ**2))"  ;
K1_PE_K2asPi_PE = K1_PE+"*"+ K2asPi_PE ;

KK_M        = "sqrt( 2*(493.67**2)         + 2*(" + K1_PE_K2_PE      + " - ( K1_PX*K2_PX + K1_PY*K2_PY+ K1_PZ*K2_PZ ) ) )"
K_KasPi_M   = "sqrt( 139.57**2 + 493.67**2 + 2*(" + K1_PE_K2asPi_PE  + " - ( K1_PX*K2_PX + K1_PY*K2_PY+ K1_PZ*K2_PZ ) ) )"

K_KasPi_PX  = "(K1_PX + K2_PX)"
K_KasPi_PY  = "(K1_PY + K2_PY)"
K_KasPi_PZ  = "(K1_PZ + K2_PZ)"
K_KasPi_PE  = "sqrt(" + K_KasPi_M +"**2 + "+ K_KasPi_PX+"**2 + " + K_KasPi_PY+"**2 + " + K_KasPi_PZ + "**2)" 
## Now Kpipi_M =
Pi_PE            = "sqrt(139.57**2 + Pi_PX**2 + Pi_PY**2 + Pi_PZ**2 )" ;
K_KasPi_PE_Pi_PE =  K_KasPi_PE+"*"+Pi_PE ;
K_KasPi_P_Pi_P   = "( " + K_KasPi_PX+"*Pi_PX +"+ K_KasPi_PY+"*Pi_PY +" + K_KasPi_PZ + "*Pi_PZ)"
K_KasPi_Pi_M     = "sqrt(139.57**2 +"+ K_KasPi_M + "**2 +2*( "+ K_KasPi_PE_Pi_PE + " - " + K_KasPi_P_Pi_P + ") )"   
###########################################################################
eos_dir_DATATuple  = "/eos/lhcb/user/b/bkhanji/data/Alldata/"
PATH_TO_File       =  eos_dir_DATATuple + "*.root";
Output_PATH        =  eos_dir_DATATuple + "TrimmedRunII.root"

MyCuts = " (B_Hlt2TopoMu2BodyDecision_TOS==1 || B_Hlt2TopoMu3BodyDecision_TOS==1 || B_Hlt2TopoMu4BodyDecision_TOS==1 || B_Hlt2Topo2BodyDecision_TOS==1 || B_Hlt2Topo3BodyDecision_TOS==1 || B_Hlt2Topo4BodyDecision_TOS==1) && Mu_ProbNNmu>0.5 && sqrt(D_PX**2+D_PY**2)>2000 && B_MCORR < 6000 && B_ENDVERTEX_CHI2<15 && D_ENDVERTEX_CHI2<15 && Dsstr_ENDVERTEX_CHI2<15 && B_DTF_TAU>0.1 && sqrt(K1_PX**2 + K1_PY**2)>600 && sqrt(K2_PX**2 + K2_PY**2)>600 && K1_IPCHI2_OWNPV>10 && K2_IPCHI2_OWNPV>10 && K1_PIDK>5 && K1_PIDK>5 && sqrt(Dsstr_PX**2+Dsstr_PY**2)>2500 && Mu_IPCHI2_OWNPV>10 && Pi1_Dsstr_IsoMinBDT>-0.95 && Pi2_Dsstr_IsoMinBDT>-0.95 && K1_IsoMinBDT>-0.95 && K2_IsoMinBDT>-0.95 && Pi_IsoMinBDT>-0.95 && nCandidate==0"

#gROOT.cd();
print "Now Trimming the content of the TChain using the following cuts : "
print '-------------------------------------------------------------------'
print '========================  Cuts applied  ===========================' 
import re
cuts_string = filter(None , re.split(r"[&&||]" , MyCuts ))
#print cuts_string
for i_string in cuts_string:
    Colorprint('| ' + i_string.ljust( len(max( cuts_string , key = len )) ) + '|' , 'blue')
print '=================================================================='
print '------------------------------------------------------------------'
# MULTIPROCESSING :
Output_PATH = '/eos/lhcb/user/b/bkhanji/data/TrimmedData/'

def RunMyJobs( PATH_TO_File  , sema ):
    print 'Worker is acquired now:'
    print '------------------------------------------------------------------'
    sema.acquire()
    outputfilename = 'Trimmed_' + PATH_TO_File 
    PATH_TO_File = eos_dir_DATATuple  + PATH_TO_File
    Colorprint ("Adding all the files in the directory : " + eos_dir_DATATuple + " in a TChain object"  , 'blue')
    InputFile             = TFile( PATH_TO_File , "OPEN") 
    chain_data            = InputFile.Get("Ds1Tuple/DecayTree");
    chain_data_SS         = InputFile.Get("Ds1Tuple_SS/DecayTree")
    chain_data_Ds1SS      = InputFile.Get("Ds1Tuple_Ds1SS/DecayTree") ;
    chain_data_Ds1SS_oppD = InputFile.Get("Ds1Tuple_Ds1SS_oppD/DecayTree")
    #print "(FYI : The TTree has "+ str(chain_data.GetNbranches()) + " branches inside ...)"
    #print "(FYI : The TTree has "+ str(chain_data.GetEntries()  ) + " of events ...)"
    RunII_trimmed = chain_data.CopyTree( MyCuts );  RunII_trimmed_SS = chain_data_SS.CopyTree( MyCuts );
    RunII_trimmed_Ds1SS = chain_data_Ds1SS.CopyTree( MyCuts );  RunII_trimmed_Ds1SS_oppD = chain_data_Ds1SS_oppD.CopyTree( MyCuts );
    Colorprint('Trimming file : '+  PATH_TO_File + ' is done.' , 'green')
    Colorprint ('Now write the file to : ' + Output_PATH , 'blue' )
    outfile = TFile( Output_PATH +  outputfilename  , "RECREATE");
    outfile.cd();
    RunII_trimmed.SetName('t');
    RunII_trimmed_SS.SetName('t_SS');
    RunII_trimmed_Ds1SS.SetName('t_SS_Ds1SS');
    RunII_trimmed_Ds1SS_oppD.SetName('t_SS_Ds1SS_oppD');
    RunII_trimmed.CloneTree().Write(); 
    RunII_trimmed_SS.CloneTree().Write();
    RunII_trimmed_Ds1SS.CloneTree().Write(); 
    RunII_trimmed_Ds1SS_oppD.CloneTree().Write();
    
    outfile.Close()
    InputFile.Close()
    Colorprint('Trimming reduced the file size by factor : ' + str(float( InputFile.GetSize()/outfile.GetSize())) , 'red')
    Colorprint('Writting ' + outputfilename + ' is done.' , 'green')
    sema.release()
    del chain_data , chain_data_SS , chain_data_Ds1SS  , chain_data_Ds1SS_oppD ,InputFile, outfile ; 
    print 'Worker is released now .... moving to the next task!'
    print '------------------------------------------------------------------'
start_time  =  time.time()
print '------------------------------------------------------------------'
print 'Nb. of CPUs on my machine are : ', multiprocessing.cpu_count() , ' CPUs'
print 'restricting myself to ', str(multiprocessing.cpu_count() -1) , ' CPUs'
print '------------------------------------------------------------------'
sema = threading.Semaphore( multiprocessing.cpu_count() -1)
i_file = 0
files_list = GetFilesinDir(eos_dir_DATATuple)

for Root_files in files_list:
    i_file +=1
    Colorprint( ' File  ' + str(i_file) + ' out of ' + str(len(files_list) ) , 'green')
    #RunMyJobs( Root_files  , 1)
    t = threading.Thread(target= RunMyJobs , args =( Root_files , sema ,) )
    t.start()
                     
print("--- %s Minutes ---" % (time.time() - start_time)/60.)

        
