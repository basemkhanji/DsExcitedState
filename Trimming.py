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
Pi_PE            = "sqrt(139.57**2 + Pi2_Dsstr_PX**2 + Pi2_Dsstr_PY**2 + Pi2_Dsstr_PZ**2 )" ;
K_KasPi_PE_Pi_PE =  K_KasPi_PE+"*"+Pi_PE ;
K_KasPi_P_Pi_P   = "( " + K_KasPi_PX+"*Pi2_Dsstr_PX +"+ K_KasPi_PY+"*Pi2_Dsstr_PY +" + K_KasPi_PZ + "*Pi2_Dsstr_PZ)"
K_KasPi_Pi_M     = "sqrt(139.57**2 +"+ K_KasPi_M + "**2 +2*( "+ K_KasPi_PE_Pi_PE + " - " + K_KasPi_P_Pi_P + ") )"   

# Now mak ethe D pipi mass :
KK_PE            = "sqrt(" + KK_M +"**2 + "+ K_KasPi_PX+"**2 + " + K_KasPi_PY+"**2 + " + K_KasPi_PZ + "**2)" 
KK_PE_Pi_PE =  KK_PE+"*"+Pi_PE ;
KK_P_Pi_P   = "( " + K_KasPi_PX+"*Pi_PX +"+ K_KasPi_PY+"*Pi_PY +" + K_KasPi_PZ + "*Pi_PZ)"
KK_Pi_M     = "sqrt(139.57**2 +"+ KK_M + "**2 +2*( "+ KK_PE_Pi_PE + " - " + KK_P_Pi_P + ") )"   

## Now Kpipipi_M =
'''
Pi1Dsstr_PE       = "sqrt(139.57**2 + Pi1_Dsstr_PX**2 + Pi1_Dsstr_PY**2 + Pi1_Dsstr_PZ**2 )" ;
K2Pi_PE          = "sqrt("+ K_KasPi_Pi_M +"**2 + D_PX**2 + D_PY**2 + D_PZ**2)"
K2Pi_PE_Pi_PE    =  K2Pi_PE+"*"+Pi1Dsstr_PE ;
K2Pi_P_Pi1Dsstr_P= "( D_PX*Pi1_Dsstr_PX + D_PY*Pi1_Dsstr_PY + D_PZ*Pi1_Dsstr_PZ)"
K2Pi_Pi1Dsstr_M  = "sqrt(139.57**2 +"+ K_KasPi_Pi_M + "**2 +2*( "+ K2Pi_PE_Pi_PE + " - " + K2Pi_P_Pi1Dsstr_P + ") )"   
'''
Pi1Dsstr_PE      = "sqrt(139.57**2 + Pi1_Dsstr_PX**2 + Pi1_Dsstr_PY**2 + Pi1_Dsstr_PZ**2 )" ;
K2Pi_PE          = "sqrt("+ K_KasPi_Pi_M +"**2 + D_PX**2 + D_PY**2 + D_PZ**2)"
K2Pi_PE_Pi_PE    =  K2Pi_PE+"*"+Pi1Dsstr_PE ;
K2Pi_P_Pi1Dsstr_P= "( D_PX*Pi1_Dsstr_PX + D_PY*Pi1_Dsstr_PY + D_PZ*Pi1_Dsstr_PZ)"
K2Pi_Pi1Dsstr_M  = "sqrt(139.57**2 +"+ K_KasPi_Pi_M + "**2 +2*( "+ K2Pi_PE_Pi_PE + " - " + K2Pi_P_Pi1Dsstr_P + ") )"   
## Now K3pi + pi mass : 
Pi2Dsstr_PE    = "sqrt(493.67**2 + Pi2_Dsstr_PX**2 + Pi2_Dsstr_PY**2 + Pi2_Dsstr_PZ**2 )" ;
K3Pi_PX        = "D_PX + Pi1_Dsstr_PX"
K3Pi_PY        = "D_PY + Pi1_Dsstr_PY"
K3Pi_PZ        = "D_PZ + Pi1_Dsstr_PZ"
K3Pi_PE        = "sqrt("+ K2Pi_Pi1Dsstr_M +"**2 + "+ K3Pi_PX + "**2 + "+K3Pi_PY+"**2 + "+ K3Pi_PZ +"**2)"
K3Pi_PE_Pi_PE  = K3Pi_PE +"*" + Pi2Dsstr_PE
K2Pi_P_Pi1Dsstr_P= "( "+ K3Pi_PX  +"*Pi2_Dsstr_PX + "+ K3Pi_PY+"*Pi2_Dsstr_PY + "+ K3Pi_PZ +"*Pi2_Dsstr_PZ)"
K3PiPi_M       = "sqrt(493.67**2 +"+ K2Pi_Pi1Dsstr_M + "**2 +2*( "+ K3Pi_PE_Pi_PE  + " - " +  K2Pi_P_Pi1Dsstr_P + ") )"   
# pi+pi- M
Pi2P2_PE          = Pi1Dsstr_PE + "*" + Pi2Dsstr_PE
pipi_M            = "sqrt(2*(139.570**2) +2*(" + Pi2P2_PE  + " - " + " ( Pi1_Dsstr_PX*Pi2_Dsstr_PX + Pi1_Dsstr_PY*Pi2_Dsstr_PY + Pi1_Dsstr_PZ*Pi2_Dsstr_PZ ) ) )" 
#
###########################################################################
# prepare the data and make the cuts :
PATH_TO_File  = "/eos/lhcb/user/b/bkhanji/data/Trimmed/";
chain_data    = TChain()
chain_data.Add( PATH_TO_File + "Trimmed_*2015*.root/t")
chain_data.Add( PATH_TO_File + "Trimmed_*2016*.root/t")
chain_data.Add( PATH_TO_File + "Trimmed_*2017*.root/t")

chain_data_SS    = TChain()
chain_data_SS.Add( PATH_TO_File + "Trimmed_2015.root/t_SS")
chain_data_SS.Add( PATH_TO_File + "Trimmed_2016.root/t_SS")
chain_data_SS.Add( PATH_TO_File + "Trimmed_2017.root/t_SS")

chain_data_Ds1SS = TChain()
chain_data_Ds1SS.Add( PATH_TO_File + "Trimmed_2015.root/t_SS_Ds1SS")
chain_data_Ds1SS.Add( PATH_TO_File + "Trimmed_2016.root/t_SS_Ds1SS")
chain_data_Ds1SS.Add( PATH_TO_File + "Trimmed_2017.root/t_SS_Ds1SS")
chain_data_Ds1SS.Add( PATH_TO_File + "Trimmed_2015.root/t_SS_Ds1SS_oppD")
chain_data_Ds1SS.Add( PATH_TO_File + "Trimmed_2016.root/t_SS_Ds1SS_oppD")
chain_data_Ds1SS.Add( PATH_TO_File + "Trimmed_2017.root/t_SS_Ds1SS_oppD")

MyCuts   = "abs(D_M - 1968.34) <30 && BDTOutput > 0.2 && Dsstr_M<2400 &&" + KK_M + ">1030 &&" +  K2Pi_Pi1Dsstr_M +">1900 "
#(B_Hlt2TopoMu2BodyDecision_TOS==1 || B_Hlt2TopoMu3BodyDecision_TOS==1 || B_Hlt2TopoMu4BodyDecision_TOS==1 || B_Hlt2Topo2BodyDecision_TOS==1 || B_Hlt2Topo3BodyDecision_TOS==1 || B_Hlt2Topo4BodyDecision_TOS==1) && abs(D_M - 1968.34) <30 && BDTOutput > 0. && Dsstr_M<2550 && " + KK_M + " >1030"
#abs("+ K2Pi_Pi1Dsstr_M + "-1986)>30" 
#abs("+ K2Pi_Pi1Dsstr_M + "-1864)>30 && abs("+K_KasPi_Pi_M +"- 1869.65)>30"   
print chain_data.GetEntries(MyCuts)
#########################################################################
import time
start_time = time.time()
#------------------------------------------------------------------------
v_min = 1950 ; v_max = 1990 
D_M_h            = TH1D("D_M_h"           , "D_M_h" , 100 , v_min , v_max) ; D_M_h.SetLineColor(kRed);
D_M_h_Ds1SS      = TH1D("D_M_h_Ds1SS"     , "D_M_h_Ds1SS" , 100 , v_min , v_max) ; D_M_h_Ds1SS.SetLineColor(kBlue);
D_M_h_Ds1SS_oppD = TH1D("D_M_h_Ds1SS_oppD"     , "D_M_h_Ds1SS_oppD" , 100 , v_min , v_max) ; D_M_h_Ds1SS_oppD.SetLineColor(kBlack);
c_D_M = TCanvas('c_D_M' , 'c_D_M' )
c_D_M.cd()
chain_data      .Draw('D_M>>D_M_h' , MyCuts  , "norm")
chain_data_Ds1SS.Draw('D_M>>D_M_h_Ds1SS' , MyCuts  , "normsame")
#chain_data      .Draw('D_M>>D_M_h_Ds1SS_oppD' ,   , "norm")
#------------------------------------------------------------------------

#------------------------------------------------------------------------
v_min = 1630 ; v_max = 2050  
kpipi_M_h            = TH1D("kpipi_M_h"           , "kpipi_M_h" , 100 , v_min , v_max) ; kpipi_M_h.SetLineColor(kRed);
kpipi_M_h_Ds1SS      = TH1D("kpipi_M_h_Ds1SS"     , "kpipi_M_h_Ds1SS" , 100 , v_min , v_max) ; kpipi_M_h_Ds1SS.SetLineColor(kBlue);
c_kpipi_M = TCanvas('c_kpipi_M' , 'c_kpipi_M' )
c_kpipi_M.cd()
chain_data      .Draw( K_KasPi_Pi_M + '>>kpipi_M_h' , MyCuts  , "norm")
chain_data_Ds1SS.Draw( K_KasPi_Pi_M + '>>kpipi_M_h_Ds1SS' , MyCuts  , "normsame")
#------------------------------------------------------------------------

#------------------------------------------------------------------------
v_min = 1550 ; v_max = 2400  
kpipipi_M_h     = TH1D("kpipipi_M_h" , "kpipipi_M_h" , 100 , v_min , v_max) ; kpipipi_M_h.SetLineColor(kRed);
kpipipi_M_h_Ds1SS     = TH1D("kpipipi_M_h_Ds1SS" , "kpipipi_M_h_Ds1SS" , 100 , v_min , v_max) ; kpipipi_M_h_Ds1SS.SetLineColor(kBlue);
kpipipi_M_h_Ds1SS_oppD = TH1D("kpipipi_M_h_Ds1SS_oppD" , "kpipipi_M_h_Ds1SS_oppD" , 100 , v_min , v_max) ; kpipipi_M_h_Ds1SS_oppD.SetLineColor(kBlack);
c_kpipipi_M     = TCanvas('c_kpipipi_M' , 'c_kpipipi_M' )
c_kpipipi_M.cd()
chain_data      .Draw( K2Pi_Pi1Dsstr_M + '>>kpipipi_M_h' , MyCuts  , "norm")
chain_data_Ds1SS.Draw( K2Pi_Pi1Dsstr_M + '>>kpipipi_M_h_Ds1SS' , MyCuts  , "normsame")
#chain_data_Ds1SS_oppD.Draw( K2Pi_Pi1Dsstr_M + '>>kpipipi_M_h_Ds1SS_oppD' , MyCuts  , "normsame")
#------------------------------------------------------------------------

#------------------------------------------------------------------------
v_min = 950 ; v_max = 1900  
KK_M_h     = TH1D("KK_M_h" , "KK_M_h" , 100 , v_min , v_max) ; KK_M_h.SetLineColor(kRed);
KK_M_h_Ds1SS= TH1D("KK_M_h_Ds1SS" , "KK_M_h_Ds1SS" , 100 , v_min , v_max) ; KK_M_h_Ds1SS.SetLineColor(kBlue);
c_KK_M     = TCanvas('c_KK_M' , 'c_KK_M' )
c_KK_M.cd()
chain_data      .Draw( KK_M + '>>KK_M_h' , MyCuts  , "norm")
chain_data_Ds1SS.Draw( KK_M + '>>KK_M_h_Ds1SS' , MyCuts  , "normsame")
#------------------------------------------------------------------------

#------------------------------------------------------------------------
v_min = 450 ; v_max = 1900  
KPi_M_h     = TH1D("KPi_M_h" , "KPi_M_h" , 100 , v_min , v_max) ; KPi_M_h.SetLineColor(kRed);
KPi_M_h_Ds1SS = TH1D("KPi_M_h_Ds1SS" , "KPi_M_h_Ds1SS" , 100 , v_min , v_max) ; KPi_M_h_Ds1SS.SetLineColor(kBlue);
c_KPi_M     = TCanvas('c_KPi_M' , 'c_KPi_M' )
c_KPi_M.cd()
chain_data      .Draw( K_KasPi_M + '>>KPi_M_h' , MyCuts  , "norm")
chain_data_Ds1SS.Draw( K_KasPi_M + '>>KPi_M_h_Ds1SS' , MyCuts  , "normsame")
#------------------------------------------------------------------------

#------------------------------------------------------------------------
v_min = 2000 ; v_max = 4000  
kpipipipi_M_h     = TH1D("kpipipipi_M_h" , "kpipipipi_M_h" , 100 , v_min , v_max) ; kpipipipi_M_h.SetLineColor(kRed);
kpipipipi_M_h_Ds1SS     = TH1D("kpipipipi_M_h_Ds1SS" , "kpipipipi_M_h_Ds1SS" , 100 , v_min , v_max) ; kpipipipi_M_h_Ds1SS.SetLineColor(kBlue);
c_kpipipipi_M     = TCanvas('c_kpipipipi_M' , 'c_kpipipipi_M' )
c_kpipipipi_M.cd()
chain_data      .Draw( K3PiPi_M + '>>kpipipipi_M_h' , MyCuts  , "norm")
chain_data_Ds1SS.Draw( K3PiPi_M + '>>kpipipipi_M_h_Ds1SS' , MyCuts  , "normsame")
#------------------------------------------------------------------------

#------------------------------------------------------------------------
v_min = 2000 ; v_max = 3000  
Dsstr_M_h     = TH1D("Dsstr_M_h" , "Dsstr_M_h" , 100 , v_min , v_max) ; Dsstr_M_h.SetLineColor(kRed);
Dsstr_M_h_Ds1SS     = TH1D("Dsstr_M_h_Ds1SS" , "Dsstr_M_h_Ds1SS" , 100 , v_min , v_max) ; Dsstr_M_h_Ds1SS.SetLineColor(kBlue);
Dsstr_M_h_Ds1SS_oppD     = TH1D("Dsstr_M_h_Ds1SS_oppD" , "Dsstr_M_h_Ds1SS_oppD" , 100 , v_min , v_max) ; Dsstr_M_h_Ds1SS_oppD.SetLineColor(kBlack);
c_Dsstr_M     = TCanvas('c_Dsstr_M' , 'c_Dsstr_M' )
c_Dsstr_M.cd()
chain_data      .Draw( 'Dsstr_M>>Dsstr_M_h' , MyCuts  , "norm")
chain_data_Ds1SS.Draw( 'Dsstr_M>>Dsstr_M_h_Ds1SS' , MyCuts  , "normsame")
#chain_data_Ds1SS_oppD.Draw( 'Dsstr_M>>Dsstr_M_h_Ds1SS_oppD' , MyCuts  , "normsame")
#------------------------------------------------------------------------

#------------------------------------------------------------------------
v_min = 400 ; v_max = 800  
pipi_M_h     = TH1D("pipi_M_h" , "pipi_M_h" , 100 , v_min , v_max) ; pipi_M_h.SetLineColor(kRed);
pipi_M_h_Ds1SS     = TH1D("pipi_M_h_Ds1SS" , "pipi_M_h_Ds1SS" , 100 , v_min , v_max) ; pipi_M_h_Ds1SS.SetLineColor(kBlue);
pipi_M_h_Ds1SS_oppD     = TH1D("pipi_M_h_Ds1SS_oppD" , "pipi_M_h_Ds1SS_oppD" , 100 , v_min , v_max) ; pipi_M_h_Ds1SS_oppD.SetLineColor(kBlack);
c_pipi_M     = TCanvas('c_pipi_M' , 'c_pipi_M' )
c_pipi_M.cd()
chain_data      .Draw( pipi_M+'>>pipi_M_h' , MyCuts  , "norm")
chain_data_Ds1SS.Draw( pipi_M+'>>pipi_M_h_Ds1SS' , MyCuts  , "normsame")
#chain_data_Ds1SS_oppD.Draw( pipi_M+'>>pipi_M_h_Ds1SS_oppD' , MyCuts  , "normsame")
#------------------------------------------------------------------------

#------------------------------------------------------------------------
v_min = 3000 ; v_max = 5500  
B_MCORR_h     = TH1D("B_MCORR_h" , "B_MCORR_h" , 100 , v_min , v_max) ; B_MCORR_h.SetLineColor(kRed);
B_MCORR_h_Ds1SS     = TH1D("B_MCORR_h_Ds1SS" , "B_MCORR_h_Ds1SS" , 100 , v_min , v_max) ; B_MCORR_h_Ds1SS.SetLineColor(kBlue);
B_MCORR_h_Ds1SS_oppD     = TH1D("B_MCORR_h_Ds1SS_oppD" , "B_MCORR_h_Ds1SS_oppD" , 100 , v_min , v_max) ; B_MCORR_h_Ds1SS_oppD.SetLineColor(kBlack);
c_B_MCORR     = TCanvas('c_B_MCORR' , 'c_B_MCORR' )
c_B_MCORR.cd()
chain_data      .Draw( 'B_MCORR>>B_MCORR_h' , MyCuts  , "norm")
chain_data_Ds1SS.Draw( 'B_MCORR>>B_MCORR_h_Ds1SS' , MyCuts  , "normsame")
#chain_data_Ds1SS_oppD.Draw( 'B_MCORR>>B_MCORR_h_Ds1SS_oppD' , MyCuts  , "normsame")
#------------------------------------------------------------------------

#------------------------------------------------------------------------
v_min = 3000 ; v_max = 5500  
B_M_h     = TH1D("B_M_h" , "B_M_h" , 100 , v_min , v_max) ; B_M_h.SetLineColor(kRed);
B_M_h_Ds1SS     = TH1D("B_M_h_Ds1SS" , "B_M_h_Ds1SS" , 100 , v_min , v_max) ; B_M_h_Ds1SS.SetLineColor(kBlue);
B_M_h_Ds1SS_oppD     = TH1D("B_M_h_Ds1SS_oppD" , "B_M_h_Ds1SS_oppD" , 100 , v_min , v_max) ; B_M_h_Ds1SS_oppD.SetLineColor(kBlack);
c_B_M     = TCanvas('c_B_M' , 'c_B_M' )
c_B_M.cd()
chain_data      .Draw( 'B_M>>B_M_h' , MyCuts  , "norm")
chain_data_Ds1SS.Draw( 'B_M>>B_M_h_Ds1SS' , MyCuts  , "normsame")
#chain_data_Ds1SS_oppD.Draw( 'B_M>>B_M_h_Ds1SS_oppD' , MyCuts  , "normsame")
#------------------------------------------------------------------------

#------------------------------------------------------------------------
outfile = TFile( 'helpfullhistos.root' , "RECREATE");
outfile.cd();
Dsstr_M_h.Write();
kpipi_M_h.Write();
kpipipi_M_h.Write();
outfile.Close()
print("--- %s seconds ---" % (time.time() - start_time))
#########################################################################
#chain_data_Ds1SS    = TChain()
#chain_data_Ds1SS.Add( PATH_TO_File + "/Ds1Tuple/DecayTree")
#
#chain_data_Ds1SS_oppD    = TChain()
#chain_data_Ds1SS_oppD.Add( PATH_TO_File + "/Ds1Tuple/DecayTree")
#c_Dsstr_M.BuildLegend(0.8 , 0.8 , 1.0 , 1.0 )
'''
'''

