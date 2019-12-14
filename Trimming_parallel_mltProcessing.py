#!/usr/bin/env python
import ROOT , sys , os , threading ,subprocess, multiprocessing , time
from ROOT import *
from math import *
from array import array
from multiprocessing import Pool
from os.path import isfile

discriminatingvariables  = [
    "log(Pi_IPChi2)",                
    "log(K1_IPChi2)",                
    "log(K2_IPChi2)",                
    "sqrt(B_PX**2+B_PY**2)",         
    "sqrt(K1_PX**2+K1_PY**2)",       
    "sqrt(Dsstr_PX**2+Dsstr_PY**2)", 
    "sqrt(K2_PX**2+K2_PY**2)",       
    "D_FDChi2",
    "Mu_ORIVX_CHI2",                 
    "K1_PIDK",                       
    "K1_PIDp",                       
    "K2_PIDp",                       
    "K2_PIDK",                             
    "sqrt(2*(139.570**2)+2*(sqrt(139.570**2+(Pi1_Dsstr_PX**2+Pi1_Dsstr_PY**2+Pi1_Dsstr_PZ**2))*sqrt(139.570**2+(Pi2_Dsstr_PX**2+Pi2_Dsstr_PY**2+Pi2_Dsstr_PZ**2))-(Pi1_Dsstr_PX*Pi2_Dsstr_PX+Pi1_Dsstr_PY*Pi2_Dsstr_PY+Pi1_Dsstr_PZ*Pi2_Dsstr_PZ)))",
    "sqrt((493.67**2)+(139.570**2)+2*(sqrt(493.67**2+(K1_PX**2+K1_PY**2+K1_PZ**2))*sqrt(139.570**2+(Pi_PX**2+Pi_PY**2+Pi_PZ**2))-(Pi_PX*K1_PX+Pi_PY*K1_PY+Pi_PZ*K1_PZ)))",
    "K2_IsoMinBDT",        
    "K1_IsoMinBDT",        
    "Pi1_Dsstr_IsoMinBDT", 
    "Pi2_Dsstr_IsoMinBDT", 
    "Pi_IsoMinBDT",        
    "Dsstr_ENDVERTEX_CHI2",
    "sqrt(2*(493.67**2)+2*(sqrt(493.67**2+(K1_PX**2+K1_PY**2+K1_PZ**2))*sqrt(493.67**2+(K2_PX**2+K2_PY**2+K2_PZ**2))-(K1_PX*K2_PX+K1_PY*K2_PY+K1_PZ*K2_PZ)))",
    "log(B_FDChi2)" ,
    "B_ENDVERTEX_CHI2",
    "log(B_IPCHI2_OWNPV)" 
    ]

discriminatingvariables_treecommand  = [
    "log(MyTrimmedTree.Pi_IPChi2)",                
    "log(MyTrimmedTree.K1_IPChi2)",                
    "log(MyTrimmedTree.K2_IPChi2)",                
    "sqrt(MyTrimmedTree.B_PX**2  + MyTrimmedTree.B_PY**2 )",         
    "sqrt(MyTrimmedTree.K1_PX**2 + MyTrimmedTree.K1_PY**2)",       
    "sqrt(MyTrimmedTree.Dsstr_PX**2+MyTrimmedTree.Dsstr_PY**2)", 
    "sqrt(MyTrimmedTree.K2_PX**2+MyTrimmedTree.K2_PY**2)",       
    "MyTrimmedTree.D_FDChi2",
    "MyTrimmedTree.Mu_ORIVX_CHI2",                 
    "MyTrimmedTree.K1_PIDK",                       
    "MyTrimmedTree.K1_PIDp",                       
    "MyTrimmedTree.K2_PIDp",                       
    "MyTrimmedTree.K2_PIDK",                             
    "sqrt(2*(139.570**2)+2*(sqrt(139.570**2+(MyTrimmedTree.Pi1_Dsstr_PX**2+MyTrimmedTree.Pi1_Dsstr_PY**2+MyTrimmedTree.Pi1_Dsstr_PZ**2))*sqrt(139.570**2+(MyTrimmedTree.Pi2_Dsstr_PX**2+MyTrimmedTree.Pi2_Dsstr_PY**2+MyTrimmedTree.Pi2_Dsstr_PZ**2))-(MyTrimmedTree.Pi1_Dsstr_PX*MyTrimmedTree.Pi2_Dsstr_PX+MyTrimmedTree.Pi1_Dsstr_PY*MyTrimmedTree.Pi2_Dsstr_PY+MyTrimmedTree.Pi1_Dsstr_PZ*MyTrimmedTree.Pi2_Dsstr_PZ)))",
    "sqrt((493.67**2)+(139.570**2)+2*(sqrt(493.67**2+(MyTrimmedTree.K1_PX**2+MyTrimmedTree.K1_PY**2+MyTrimmedTree.K1_PZ**2))*sqrt(139.570**2+(MyTrimmedTree.Pi_PX**2+MyTrimmedTree.Pi_PY**2+MyTrimmedTree.Pi_PZ**2))-(MyTrimmedTree.Pi_PX*MyTrimmedTree.K1_PX+MyTrimmedTree.Pi_PY*MyTrimmedTree.K1_PY+MyTrimmedTree.Pi_PZ*MyTrimmedTree.K1_PZ)))",
    "MyTrimmedTree.K2_IsoMinBDT",        
    "MyTrimmedTree.K1_IsoMinBDT",        
    "MyTrimmedTree.Pi1_Dsstr_IsoMinBDT", 
    "MyTrimmedTree.Pi2_Dsstr_IsoMinBDT", 
    "MyTrimmedTree.Pi_IsoMinBDT",        
    "MyTrimmedTree.Dsstr_ENDVERTEX_CHI2",
    "sqrt(2*(493.67**2)+2*(sqrt(493.67**2+(MyTrimmedTree.K1_PX**2+MyTrimmedTree.K1_PY**2+MyTrimmedTree.K1_PZ**2))*sqrt(493.67**2+(MyTrimmedTree.K2_PX**2+MyTrimmedTree.K2_PY**2+MyTrimmedTree.K2_PZ**2))-(MyTrimmedTree.K1_PX*MyTrimmedTree.K2_PX+MyTrimmedTree.K1_PY*MyTrimmedTree.K2_PY+MyTrimmedTree.K1_PZ*MyTrimmedTree.K2_PZ)))",
    "log(MyTrimmedTree.B_FDChi2)" ,
    "MyTrimmedTree.B_ENDVERTEX_CHI2",
    "log(MyTrimmedTree.B_IPCHI2_OWNPV)" 
    ] 
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
    print 'Looking for all *.root files in : ' , mypath
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
def RunMyJobs( PATH_TO_File ):
    #sema.acquire()
    outputfilename = 'Trimmed_' + PATH_TO_File 
    PATH_TO_File = eos_dir_DATATuple  + PATH_TO_File
    Colorprint ("Opening the file: " + PATH_TO_File + " "  , 'blue')
    InputFile             = TFile( PATH_TO_File , "OPEN") 
    chain_data            = InputFile.Get("Ds1Tuple/DecayTree");
    chain_data_SS         = InputFile.Get("Ds1Tuple_SS/DecayTree")
    chain_data_Ds1SS      = InputFile.Get("Ds1Tuple_Ds1SS/DecayTree") ;
    chain_data_Ds1SS_oppD = InputFile.Get("Ds1Tuple_Ds1SS_oppD/DecayTree")
    #print "(FYI : The TTree has "+ str(chain_data.GetNbranches()) + " branches inside ...)"
    #print "(FYI : The TTree has "+ str(chain_data.GetEntries()  ) + " of events ...)"
    RunII_trimmed = chain_data.CopyTree( MyCuts );  RunII_trimmed_SS = chain_data_SS.CopyTree( MyCuts );
    RunII_trimmed_Ds1SS = chain_data_Ds1SS.CopyTree( MyCuts );  RunII_trimmed_Ds1SS_oppD = chain_data_Ds1SS_oppD.CopyTree( MyCuts );
    RunII_trimmed.SetName('t');
    RunII_trimmed_SS.SetName('t_SS');
    RunII_trimmed_Ds1SS.SetName('t_SS_Ds1SS');
    RunII_trimmed_Ds1SS_oppD.SetName('t_SS_Ds1SS_oppD');
    Colorprint('Trimming file : '+  PATH_TO_File + ' is done.' , 'green')
    # To add the BDT: First, clone the trimmed tree (not CPU-expensive)
    tree_list = [ RunII_trimmed, RunII_trimmed_SS, RunII_trimmed_Ds1SS, RunII_trimmed_Ds1SS_oppD ]
    BDTtrees_list = []
    Colorprint ('Now write the file to : ' + Output_PATH , 'blue' )
    outfile = TFile( Output_PATH +  outputfilename  , "RECREATE");
    outfile.cd();
    for MyTrimmedTree in tree_list:
        BDTreader = TMVA.Reader("!Color:!Silent")
        Vars = []
        n = 0
        for var in discriminatingvariables:
            exec('var'+str(n)+' = array(\'f\',[0])')
            #print 'variable name : ' , var
            exec('BDTreader.AddVariable("'+var+'",var'+str(n)+')')
            exec('Vars.append(var'+str(n)+')')
            n += 1

        BDTreader.BookMVA('BDTG','./TMVAClassification_BDTG.weights.xml')
        MyBDTTree = MyTrimmedTree.CopyTree('0')
        
        BDTOutput = array( 'f', [ 0 ] ) #numpy.zeros(1 , dtype=float)
        MyBDTTree.Branch('BDTOutput', BDTOutput,' BDTOutput/F')
        N = MyTrimmedTree.GetEntries()
        for n in range(N):
            #print 'Entry : ' , n , ' in Tree: ' , MyTrimmedTree.GetName()
            MyTrimmedTree.GetEntry(n)
            a = 0
            for var in discriminatingvariables_treecommand:
                #print 'var expression ' , var
                exec('var'+str(a)+'[0] = '+var)
                a += 1
            BDTOutput[0] = BDTreader.EvaluateMVA('BDTG')
            #print 'BDT value = ' , BDTOutput[0]
            MyBDTTree.Fill()
        BDTtrees_list.append(MyBDTTree)    
        MyBDTTree.Write("", TObject.kOverwrite)

    outfile.Close()
    InputFile.Close()
    Colorprint('Trimming reduced the file size by factor : ' + str(float( InputFile.GetSize()/outfile.GetSize())) , 'red')
    Colorprint('Writting ' + outputfilename + ' is done.' , 'green')
    for bdttree_i in BDTtrees_list:
        del bdttree_i
    for chain_data_i in tree_list:
        del chain_data_i
    del chain_data , chain_data_SS , chain_data_Ds1SS  , chain_data_Ds1SS_oppD ,InputFile, outfile ;
    print 'Worker is released now .... moving to the next task!'
    print '------------------------------------------------------------------'
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
Output_PATH        =  eos_dir_DATATuple + "Trimmed.root"

#MyCuts = "(B_Hlt2TopoMu2BodyDecision_TOS==1 || B_Hlt2TopoMu3BodyDecision_TOS==1 || B_Hlt2TopoMu4BodyDecision_TOS==1 || B_Hlt2Topo2BodyDecision_TOS==1 || B_Hlt2Topo3BodyDecision_TOS==1 || B_Hlt2Topo4BodyDecision_TOS==1) && Mu_ProbNNmu>0.5 && sqrt(D_PX**2+D_PY**2)>2000 && B_MCORR < 6000 && B_ENDVERTEX_CHI2<15 && D_ENDVERTEX_CHI2<15 && Dsstr_ENDVERTEX_CHI2<15 && B_DTF_TAU>0.1 && sqrt(K1_PX**2 + K1_PY**2)>600 && sqrt(K2_PX**2 + K2_PY**2)>600 && K1_IPCHI2_OWNPV>10 && K2_IPCHI2_OWNPV>10 && K1_PIDK>5 && K1_PIDK>5 && sqrt(Dsstr_PX**2+Dsstr_PY**2)>2500 && Mu_IPCHI2_OWNPV>10 && Pi1_Dsstr_IsoMinBDT>-0.95 && Pi2_Dsstr_IsoMinBDT>-0.95 && K1_IsoMinBDT>-0.95 && K2_IsoMinBDT>-0.95 && Pi_IsoMinBDT>-0.95 && nCandidate==0"

MyCuts = "Mu_ProbNNmu>0.5 && sqrt(D_PX**2+D_PY**2)>2000 && B_MCORR < 6000 && B_ENDVERTEX_CHI2<15 && D_ENDVERTEX_CHI2<15 && Dsstr_ENDVERTEX_CHI2<15 && B_DTF_TAU>0.1 && sqrt(K1_PX**2 + K1_PY**2)>600 && sqrt(K2_PX**2 + K2_PY**2)>600 && K1_IPCHI2_OWNPV>10 && K2_IPCHI2_OWNPV>10 && K1_PIDK>5 && K1_PIDK>5 && sqrt(Dsstr_PX**2+Dsstr_PY**2)>2500 && Mu_IPCHI2_OWNPV>10 && Pi1_Dsstr_IsoMinBDT>-0.95 && Pi2_Dsstr_IsoMinBDT>-0.95 && K1_IsoMinBDT>-0.95 && K2_IsoMinBDT>-0.95 && Pi_IsoMinBDT>-0.95 && nCandidate==0"

print '-------------------------------------------------------------------'
Colorprint( '========================  Cuts applied  ===========================' , 'red')
import re
cuts_string = filter(None , re.split(r"[&&||]" , MyCuts ))
for i_string in cuts_string:
    Colorprint('| ' + i_string.ljust( len(max( cuts_string , key = len )) ) + '|' , 'blue')
Colorprint('==================================================================', 'red')
print '------------------------------------------------------------------'
# I/O
Output_PATH = '/eos/lhcb/user/b/bkhanji/data/Trimmed/'
print '------------------------------------------------------------------'
print 'Nb. of CPUs on my machine are : ', multiprocessing.cpu_count() , ' CPUs'
print 'restricting myself to ', str(multiprocessing.cpu_count() -1) , ' CPUs'
print '------------------------------------------------------------------'
pool = Pool(processes=  multiprocessing.cpu_count() -1 )
i_file = 0
                     
# 1- Get the Tuples:  
files_list = GetFilesinDir(eos_dir_DATATuple)
#2 - Trimm in a loop allowing Multiprocessing
start_time = time.time()
for Root_files in files_list:
    i_file +=1
    Colorprint( ' File  ' + str(i_file) + ' out of ' + str(len(files_list) ) , 'green')
    pool.apply_async( RunMyJobs , ( Root_files , ) )
pool.close()
pool.join()
print("Job took :" , ( time.time() - start_time )/3600 , " Hours")

