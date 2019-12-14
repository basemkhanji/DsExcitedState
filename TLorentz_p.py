#!/usr/bin/env python
import ROOT , sys , os , threading ,subprocess, multiprocessing , time
from ROOT import *
from math import *
from array import array
from multiprocessing import Pool
from os.path import isfile
###########################################################################
PATH_TO_File  = "/eos/lhcb/user/b/bkhanji/data/Trimmed/Merged/";
chain_data    = TChain()
#chain_data.Add( PATH_TO_File + "Trimmed_9_DTT_2015_Reco15aStrip24r1_Up_SEMILEPTONIC.DST.root/t")
#chain_data.Add( PATH_TO_File + "Trimmed_*_2015_*.root/t")
#chain_data.Add( PATH_TO_File + "Trimmed_2011_Up.root/t")
#chain_data.Add( PATH_TO_File + "Trimmed_2011_Down.root/t")

chain_data.Add( PATH_TO_File + "Trimmed_*2015*.root/t")
chain_data.Add( PATH_TO_File + "Trimmed_*2016*.root/t")
chain_data.Add( PATH_TO_File + "Trimmed_*2017*.root/t")

print 'added all files ... \n'

'''
chain_data_SS    = TChain()
chain_data_SS.Add( PATH_TO_File + "Trimmed_2015.root/t_SS")
chain_data_SS.Add( PATH_TO_File + "Trimmed_2011.root/t_SS")
chain_data_SS.Add( PATH_TO_File + "Trimmed_2017.root/t_SS")

chain_data_Ds1SS = TChain()
chain_data_Ds1SS.Add( PATH_TO_File + "Trimmed_2015.root/t_SS_Ds1SS")
chain_data_Ds1SS.Add( PATH_TO_File + "Trimmed_2016.root/t_SS_Ds1SS")
chain_data_Ds1SS.Add( PATH_TO_File + "Trimmed_2017.root/t_SS_Ds1SS")
chain_data_Ds1SS.Add( PATH_TO_File + "Trimmed_2015.root/t_SS_Ds1SS_oppD")
chain_data_Ds1SS.Add( PATH_TO_File + "Trimmed_2016.root/t_SS_Ds1SS_oppD")
chain_data_Ds1SS.Add( PATH_TO_File + "Trimmed_2017.root/t_SS_Ds1SS_oppD")
'''

#MyCuts   = "abs(D_M - 1968.34) <30 && BDTOutput > 0.2 && Dsstr_M<2400"
MyCuts   = "abs(D_M - 1968.34) <30 && BDTOutput > 0.2 && Dsstr_M>2400 && Dsstr_M<2600" # cuts fpor Marco
#(B_Hlt2TopoMu2BodyDecision_TOS==1 || B_Hlt2TopoMu3BodyDecision_TOS==1 || B_Hlt2TopoMu4BodyDecision_TOS==1 || B_Hlt2Topo2BodyDecision_TOS==1 || B_Hlt2Topo3BodyDecision_TOS==1 || B_Hlt2Topo4BodyDecision_TOS==1) 
#print chain_data.GetEntries(MyCuts)
#########################################################################
import time
start_time = time.time()
print 'Looping over events'

Default_l = [
    ('K1'        , 493.677 ,  +1 ) ,
    ('K2'        , 493.677 ,  -1 ) ,
    ('Pi'        , 139.577 ,  -1 ) ,
    ('Pi1_Dsstr' , 139.570 ,  +1 ) ,
    ('Pi2_Dsstr' , 139.570 ,  -1 ) ,
    ('Mu'        , 139.570 ,  +1 ) ,
    # MisID particles 
    ]

#--------------------------------
# Make combinatorix 
Combinedparticle_list = []
from itertools import combinations
import itertools
for i in range( 2 , len( Default_l  ) +1 ):
    Combinedparticle_list.extend( itertools.combinations( Default_l , i) )

for r in list(Combinedparticle_list):
    sum_charges =  sum(j[2] for j in list(r))
    if abs(sum_charges) > 1: 
        Combinedparticle_list.remove(r)
        
nentries = chain_data.GetEntries()

def loopOverMyEvents( i , chain_data  = chain_data  , Combinedparticle_list = Combinedparticle_list):
    #print 'Hi I am here : ' ,i
    
    chain_data.GetEntry(i)
    v_K1           = TLorentzVector()
    v_K2           = TLorentzVector()
    v_K2_asPi      = TLorentzVector()
    v_Pi2_Dsstr    = TLorentzVector()
    v_Pi1_Dsstr    = TLorentzVector()
    v_Pi           = TLorentzVector()
    v_Mu           = TLorentzVector()
    #v_K1.SetXYZM(chain_data.K1_PX,chain_data.K1_PY,chain_data.K1_PZ  , 139.570 )
    v_K1.SetXYZM(chain_data.K1_PX,chain_data.K1_PY,chain_data.K1_PZ  , 493.677 )
    #v_K2.SetXYZM(chain_data.K2_PX,chain_data.K2_PY,chain_data.K2_PZ  , 493.677 )
    #v_K2.SetXYZM(chain_data.K2_PX,chain_data.K2_PY,chain_data.K2_PZ  ,  139.570 )
    v_K2.SetXYZM(chain_data.K2_PX,chain_data.K2_PY,chain_data.K2_PZ  , 938.272 )
    v_Pi.SetXYZM(chain_data.Pi_PX,chain_data.Pi_PY,chain_data.Pi_PZ  , 139.570 )
    v_Pi1_Dsstr.SetXYZM(chain_data.Pi1_Dsstr_PX,chain_data.Pi1_Dsstr_PY,chain_data.Pi1_Dsstr_PZ , 139.570 )
    v_Pi2_Dsstr.SetXYZM(chain_data.Pi2_Dsstr_PX,chain_data.Pi2_Dsstr_PY,chain_data.Pi2_Dsstr_PZ , 139.570 )
    v_Mu.SetXYZM(chain_data.Mu_PX,chain_data.Mu_PY,chain_data.Mu_PZ  , 105.658 )

    v_K2_asPi.SetXYZM(chain_data.K2_PX,chain_data.K2_PY,chain_data.K2_PZ  , 139.570 )
    v_Dstar = v_K1 + v_K2_asPi + v_Pi2_Dsstr + v_Pi1_Dsstr + v_Pi
    v_Dzero = v_K1             + v_Pi2_Dsstr + v_Pi1_Dsstr + v_Pi
    Dstar_M = v_Dstar.M()
    Dzero_M = v_Dzero.M()
    v_phi= v_K1 + v_K2
    loopprogress = 100*round(float(i)/nentries ,1)
    #print loopprogress
    sys.stdout.write('\rprogress : '+str(loopprogress)+' % ( '+str(i)+' out of '+str(nentries) + ' )')
    sys.stdout.flush()
    var_l = []    
    if (abs(chain_data.D_M -1968.34) <30 and \
        chain_data.BDTOutput > 0.3 ) and \
        chain_data.Dsstr_M > 2520 ) and \
        chain_data.Dsstr_M < 2560 ): #and \
        #abs(Dstar_M-2010) > 25 and \
        #abs(Dzero_M-1864) > 25 and \
        #chain_data.K2_ProbNNp<0.4 ):
        for r in Combinedparticle_list:
            SumTvector = ''.join( particle[0] for particle in list(r)) + '_v '
            InvMassvar = ''.join( particle[0] for particle in list(r)) + '_m '
            sum_r   = SumTvector + ' = ' +  '+'.join( 'v_' + particle[0] for particle in list(r)) 
            exec( sum_r) 
            exec( InvMassvar +" = " + SumTvector+".M()" )
            exec( "var_l.append("+ InvMassvar  +") ")
    
    return var_l
    
starttime = time.time()        

pool     = Pool(processes=  multiprocessing.cpu_count() - 5 )
result_l = pool.map( loopOverMyEvents , range( nentries ) )
pool.close()
pool.join()

#for i in range( nentries ):
#    loopOverMyEvents(i)

outfile = TFile( 'helpfullhistos.root' , "RECREATE");
outfile.cd();
# filter out empty entries:
#result_l = [x for x in list(result_l) if x is not None]
result_l = filter(None, result_l)

for i_tpl in Combinedparticle_list:
    print '------------------------------------------------------start' 
    h_name   = ''.join( particle[0] for particle in list(i_tpl)) + '_M'
    sum_mass =  sum(mass_i[1] for mass_i in list(i_tpl))
    print 'histo_name is ' , h_name
    print 'mass limits : ', sum_mass , ' , ' , 2.5*sum_mass 
    canvas_e = TCanvas(h_name + '_c' , h_name + '_c') 
    canvas_e.cd()
    h_e =  TH1D( h_name , h_name , 100, 0 , 0 ) 
    #canvas_l.append(TCanvas(h_name + '_c' , h_name + '_c')  )
    #h_l.append(TH1D( h_name , h_name , 100, sum_mass , 2.5*sum_mass))
    for a in result_l:
        #print Combinedparticle_list.index(i_tpl)
        #print a[Combinedparticle_list.index(i_tpl)]
        #print '------------------------------------------------------end'
        mymass = a[Combinedparticle_list.index(i_tpl)]
        h_e.Fill( a[Combinedparticle_list.index(i_tpl)] )

    h_e.GetXaxis().SetTitle(h_name + ' #[MeV/C^2]')
    h_e.Draw()
    canvas_e.Modified()
    canvas_e.Update()
    h_e.Write()
    canvas_e.Write()
    
outfile.Close()    

endtime = time.time()
#c1 = TCanvas()
#h_l[1].Draw()

print 'code took : ' , endtime - starttime , ' seconds'
#MyDp_M_h.SetLineColor(kRed) ; MyDp_M_h.Draw()
#Dp_M_h.SetLineColor(kBlue) ;Dp_M_h.Draw("same")




'''
print 'closing the loop'




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


