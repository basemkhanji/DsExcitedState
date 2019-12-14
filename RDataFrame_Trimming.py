#!/usr/bin/env python
import ROOT
from ROOT import *
from math import *
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
PATH_TO_File       =  eos_dir_DATATuple + "11*2015*Up*.root";
Output_PATH        =  eos_dir_DATATuple + "TrimmedRunII_TDF.root"

#MyCuts = " (B_Hlt2TopoMu2BodyDecision_TOS==1 | B_Hlt2TopoMu3BodyDecision_TOS==1 | B_Hlt2TopoMu4BodyDecision_TOS==1 | B_Hlt2Topo2BodyDecision_TOS==1 | B_Hlt2Topo3BodyDecision_TOS==1 | B_Hlt2Topo4BodyDecision_TOS==1) & Mu_ProbNNmu>0.5 & sqrt(pow(D_PX,2)+ pow(D_PY,2))>2000 & B_MCORR < 6000 & B_ENDVERTEX_CHI2<15 & D_ENDVERTEX_CHI2<15 & Dsstr_ENDVERTEX_CHI2<15 & B_DTF_TAU>0.1 & sqrt(pow(K1_PX , 2) + pow(K1_PY,2))>600 & sqrt(pow(K2_PX,2) + pow(K2_PY,2))>600 & K1_IPCHI2_OWNPV>10 & K2_IPCHI2_OWNPV>10 & K1_PIDK>5 & K1_PIDK>5 & sqrt(pow(Dsstr_PX,2)+ pow(Dsstr_PY,2))>2500 & Mu_IPCHI2_OWNPV>10 & Pi1_Dsstr_IsoMinBDT>-0.95 & Pi2_Dsstr_IsoMinBDT>-0.95 & K1_IsoMinBDT>-0.95 & K2_IsoMinBDT>-0.95 & Pi_IsoMinBDT>-0.95 & nCandidate==0"

MyCuts = "nCandidate==0"

Colorprint ("Adding all the files in the directory : " + eos_dir_DATATuple + " in a TChain object"  , 'blue')
ROOT.ROOT.EnableImplicitMT()
RDataFrame_RunIIdata    = ROOT.ROOT.RDataFrame( "Ds1Tuple/DecayTree" , '/eos/lhcb/user/b/bkhanji/data/Alldata/119_DTT_2017_Reco17Strip29r2_Up_SEMILEPTONIC.DST.root' )
#print RDataFrame_RunIIdata.Count().GetValue()
#RDataFrame_RunIIdata.Add( PATH_TO_File + "/Ds1Tuple/DecayTree")
print "Now Trimming the content of the TChain using the following cuts : "
print '-------------------------------------------------------------------'
Colorprint( '======================== ' + "Cuts applied" + ' ===========================' , 'blue')
import re
cuts_string = filter(None , re.split(r"[&|]" , MyCuts ))
#print cuts_string
for i_string in cuts_string:
    Colorprint('| ' + i_string.ljust( len(max( cuts_string , key = len )) ) + '|' , 'blue')
print '=================================================================='
print '------------------------------------------------------------------'
import time
start_time = time.time()

RunII_trimmed_bool = RDataFrame_RunIIdata.Filter( MyCuts  )
print 'Filtering done!'
RunII_trimmed_bool.Snapshot( 'Ds1Tuple' , Output_PATH  )
Colorprint('Copying is done.' , 'green')
print("--- %s seconds ---" % (time.time() - start_time))


