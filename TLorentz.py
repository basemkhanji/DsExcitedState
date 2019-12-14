#!/usr/bin/env python
import ROOT
from ROOT import *
from math import *
import sys
import time
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
PATH_TO_File  = "/eos/lhcb/user/b/bkhanji/data/Trimmed/";
chain_data    = TChain()
chain_data.Add( PATH_TO_File + "Trimmed_2015.root/t")
chain_data.Add( PATH_TO_File + "Trimmed_2016.root/t")
chain_data.Add( PATH_TO_File + "Trimmed_2017.root/t")

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

MyCuts   = "abs(D_M - 1968.34) <30 && BDTOutput > 0.2 && Dsstr_M<2400" 
#(B_Hlt2TopoMu2BodyDecision_TOS==1 || B_Hlt2TopoMu3BodyDecision_TOS==1 || B_Hlt2TopoMu4BodyDecision_TOS==1 || B_Hlt2Topo2BodyDecision_TOS==1 || B_Hlt2Topo3BodyDecision_TOS==1 || B_Hlt2Topo4BodyDecision_TOS==1) 
print chain_data.GetEntries(MyCuts)
#########################################################################
import time
start_time = time.time()
print 'Looping over events'

Default_l = [
    ('K1'        , 493.677 ,  +1 ) ,
    ('K2'        , 493.677 ,  -1 ) ,
    ('Pi'        , 139.577 ,  +1 ) ,
    ('Pi1_Dsstr' , 139.570 ,  +1 ) ,
    ('Pi2_Dsstr' , 139.570 ,  -1 ) ,
    # MisID particles 
    ]

h_l    = []
Combinedparticle_list = []

from itertools import combinations
import itertools
for i in range( 2 , len( Default_l  ) +1 ):
    Combinedparticle_list.extend( itertools.combinations( Default_l , i) )

for r in list(Combinedparticle_list):
    sum_charges =  sum(j[2] for j in list(r))
    if abs(sum_charges) > 1: 
        Combinedparticle_list.remove(r)



for i_tpl in Combinedparticle_list:
    h_name   = ''.join( particle[0] for particle in list(i_tpl)) + '_M'
    sum_mass =  sum(mass_i[1] for mass_i in list(i_tpl))
    print 'histo_name is ' , h_name
    print 'mass limits : ', sum_mass , ' , ' , 2.5*sum_mass 
    h_invm = TH1D( h_name , h_name , 100, sum_mass , 2.5*sum_mass)
    h_l.append(h_invm)

pool = Pool(processes=  multiprocessing.cpu_count() -1 )
nentries = chain_data.GetEntries()
MyDp_M_h = TH1D( 'MyDp_M' , 'MyDp_M' , 100, 1800 , 1900 )
Dp_M_h = TH1D( 'Dp_M' , 'Dp_M' , 100, 1800 , 1900 )    
for i in range( nentries ):
    #print 'Hi I am here : ' ,i
    chain_data.GetEntry(i)
    v_K1     = TLorentzVector()
    v_K2     = TLorentzVector()
    v_phi    = TLorentzVector()
    v_Pi     = TLorentzVector()
    v_Dp     = TLorentzVector()
    
    v_K1.SetXYZM(chain_data.K1_PX,chain_data.K1_PY,chain_data.K1_PZ, 493.677 )
    v_K2.SetXYZM(chain_data.K2_PX,chain_data.K2_PY,chain_data.K2_PZ, 493.677 )
    v_Pi.SetXYZM(chain_data.Pi_PX,chain_data.Pi_PY,chain_data.Pi_PZ, 139.570 )
    v_phi.SetXYZM(chain_data.Pi_PX,chain_data.Pi_PY,chain_data.Pi_PZ, 139.570 )
    v_phi= v_K1 + v_K2
    v_Dp = v_K1 + v_K2 + v_Pi
    loopprogress = round(float(i)/nentries ,1)
    #print loopprogress
    sys.stdout.write('\rprogress : '+str(loopprogress)+' % ( '+str(i)+' out of '+str(nentries) + ' )')
    sys.stdout.flush()
        
    if (abs(chain_data.D_M -1968.34) <30 and \
        chain_data.BDTOutput > 0.1 and \
        chain_data.Dsstr_M <3000 and \
        abs(v_phi.M() -1020)<5 ) :
        
        v_K1.SetXYZM(chain_data.K1_PX,chain_data.K1_PY,chain_data.K1_PZ, 493.677 )
        v_K2.SetXYZM(chain_data.K2_PX,chain_data.K2_PY,chain_data.K2_PZ, 139.570 )
        v_Pi.SetXYZM(chain_data.Pi_PX,chain_data.Pi_PY,chain_data.Pi_PZ, 139.570 )
        v_Dp = v_K1 + v_K2 + v_Pi
        
        MyDp_M_h.Fill(v_Dp.M())
        Dp_M_h.Fill(chain_data.D_M)
        
        
c1 = TCanvas()
MyDp_M_h.SetLineColor(kRed) ; MyDp_M_h.Draw()
Dp_M_h.SetLineColor(kBlue) ;Dp_M_h.Draw("same")


print 'closing the loop'

'''
    for particles1 in Default_l:
        for particles2 in Default_l:
            if particles2 == particles1:
                'same particle , skipping ...'
            print particles1
            print particles2
            exec('v_particle.SetXYZM(chain_data.'+ particles1  + '_PX ,  chain_data.' +  particles1 +'_PY , chain_data.' + particles1 + '_PZ , ' + mass_particle + ' ) ')
            exec('v_particle.SetXYZM(chain_data.'+ particles2  + '_PX ,  chain_data.' +  particles2 +'_PY , chain_data.' + particles2 + '_PZ , ' + mass_particle + ' ) ')
            
            print v_particle.E()
            v_particle



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


