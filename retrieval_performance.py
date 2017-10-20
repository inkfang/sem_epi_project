# -*- coding: utf-8 -*-
"""
Created on Tue Apr 15 11:37:35 2014

@author: jing
"""

##############################################################################################
# This code generate retrieval result for two object with diff noise level, noise drift,     #
# and jumping within and between seq.                                                        #
##############################################################################################

# impact of number of storing seq#
import pypar
myid = pypar.rank()

from numpy import *
from pylab import *
import mdp
import scipy as sy
import inpobject as lsau
import time
#import randomwalkgenerator as rw
#import altiinput as al


def hierachynet(bo,recf,ovl):  #correct
    switchboard = mdp.hinet.Rectangular2dSwitchboard(in_channels_xy = (bo,bo),
                                                 field_channels_xy=(recf,recf),
                                                 field_spacing_xy =(ovl,ovl))

    sfa_dim = 48
    sfa_lower_out= 32 #64

    sfanode = mdp.nodes.SFANode(input_dim = switchboard.out_channel_dim, output_dim = sfa_dim)
    sfa2node = mdp.nodes.QuadraticExpansionNode(input_dim=sfa_dim)
    #noi_node = mdp.nodes.NoiseNode(input_dim = sfa2node.output_dim,noise_args=(0,sqrt(0.0005)))   #test#
    sfanode2 = mdp.nodes.SFANode(input_dim = sfa2node.output_dim,output_dim = sfa_lower_out)

    flownode = mdp.hinet.FlowNode(mdp.Flow([sfanode,sfa2node,sfanode2]))
    #flownode = mdp.hinet.FlowNode(mdp.Flow([sfanode,sfa2node,noi_node,sfanode2])) #test#

    sfalayer = mdp.hinet.CloneLayer(flownode, n_nodes = switchboard.output_channels)
    flow = mdp.Flow([switchboard, sfalayer])

    sfa_upper_bo = 6
    sfa_upper_recf = 4
    sfa_upper_out = 32 
    sfa_top_out = 10
    ovl2 = 2  

    switchboard2 =  mdp.hinet.Rectangular2dSwitchboard(in_channels_xy = (sfa_upper_bo,sfa_upper_bo),
                                                 field_channels_xy=(sfa_upper_recf,sfa_upper_recf),
                                                 field_spacing_xy = (ovl2,ovl2),  # new adding
                                                 in_channel_dim = sfa_lower_out)                        
#                                                
    sfa_uppernode = mdp.nodes.SFANode(input_dim = switchboard2.out_channel_dim, output_dim = sfa_dim)
    sfa_upperexp = mdp.nodes.QuadraticExpansionNode(input_dim = sfa_dim)
    #noi_node = mdp.nodes.NoiseNode(input_dim = sfa_upperexp.output_dim,noise_args=(0,sqrt(0.0005)))   #test#

    sfa_uppernode2 = mdp.nodes.SFANode(input_dim = sfa_upperexp.output_dim, output_dim = sfa_upper_out)
    upper_flownode = mdp.hinet.FlowNode(mdp.Flow([sfa_uppernode,sfa_upperexp,sfa_uppernode2]))
    #upper_flownode = mdp.hinet.FlowNode(mdp.Flow([sfa_uppernode,sfa_upperexp,noi_node,sfa_uppernode2])) #test#
    upper_sfalayer = mdp.hinet.CloneLayer(upper_flownode, n_nodes = switchboard2.output_channels)

    
    #sfa_top_node = mdp.nodes.SFANode(input_dim = switchboard2.out_channel_dim, output_dim = sfa_top_out) #mistake
    sfa_top_node = mdp.nodes.SFANode(input_dim = upper_sfalayer.output_dim, output_dim = sfa_dim)
    sfa_topexp =mdp.nodes.QuadraticExpansionNode(input_dim = sfa_dim)
    #noi_node = mdp.nodes.NoiseNode(input_dim = sfa_topexp.output_dim,noise_args=(0,sqrt(0.0005)))   #test#

    sfa_topnode2 =mdp.nodes.SFANode(input_dim = sfa_topexp.output_dim,output_dim = sfa_top_out)
    sfa_over_node = mdp.hinet.FlowNode(mdp.Flow([sfa_top_node,sfa_topexp,sfa_topnode2]))
    #sfa_over_node = mdp.hinet.FlowNode(mdp.Flow([sfa_top_node,sfa_topexp,noi_node,sfa_topnode2])) #test#
    
    network = mdp.Flow([switchboard,sfalayer,switchboard2,upper_sfalayer,sfa_over_node])

    return network                                                 



def noseqN2(stoep,lseq,lseg,q,di,sd):

    reseq = zeros((lseq,di))
    p = q+normal(0,(sd+1e-20),di)


    for j in range(lseq-1):

        p=q+normal(0,(sd+1e-20),di)

        #disma = sy.spatial.distance.cdist(qq,stoep[:,:di])[0]
        disma = sqrt(sum((stoep[:,:di]-p)**2,1))
        mind = where(disma == min(disma))[0]
        shuffle(mind)
        ind = mind[0]

        reseq[j] = stoep[ind,:di]
        q = stoep[ind,lseg*di-di::]
        #qq = (q+normal(0,(sd+1e-20),di))*ones([len(stoep),len(q)])

    reseq[lseq-1] = stoep[ind,lseg*di-di::]
    return reseq



# generate pieces #
def genpiece(slow,lseg):
    seg = mdp.nodes.TimeFramesNode(lseg)                 # only work for 2D data
    stoep = seg(slow)
    
    return stoep

    
def distancewithoutk(re,slow,stoep,icue,lseg,di,nosd):
    #l = zeros(re)
    reseq = noseqN2(stoep,re,lseg,slow[0],di,nosd)
    l = sqrt(sum((slow-reseq)**2,1))
    
    return l
### drift ###    
def noseqwithoutseg(stoep,lseq,lseg,q,di,sd):

    reseq = zeros((lseq,di))
    p = q+normal(0,(sd+1e-20),di)

    for j in range(lseq):
        disma = sqrt(sum((stoep[:,:di]-p)**2,1))
        mind = where(disma == min(disma))[0]
        shuffle(mind)
        ind = mind[0]

        reseq[j] = stoep[ind,:di]
        q = stoep[ind,:di]
        p = q+normal(0,(sd+1e-20),di)

    return reseq

def distdrift(re,slow,stoep,icue,lseg,di,nosd):
    reseq = noseqwithoutseg(stoep,re,lseg,slow[0],di,nosd)
    l = sqrt(sum((reseq - slow[0])**2,1))
    #cl = mean(sqrt(sum((reseq[1::,:] - reseq[:49,:])**2,1)))

    return l#,cl

### jump in seq ###
def noseqnew(stoep,lseq,lseg,q,di,sd):

    reseq = zeros((lseq,di))

    #qq = (q+normal(0,(sd+1e-20),di))*ones([len(stoep),len(q)])
    p = q+normal(0,(sd+1e-20),di)
    reind = zeros(lseq-1)

    for j in range(lseq-1):

        p=q+normal(0,(sd+1e-20),di)

        disma = sqrt(sum((stoep[:,:di]-p)**2,1))
        mind = where(disma == min(disma))[0]
        shuffle(mind)
        ind = mind[0]

        reseq[j] = stoep[ind,:di]
        q = stoep[ind,lseg*di-di::]
        reind[j] =  ind
    reseq[lseq-1] = stoep[ind,lseg*di-di::]
    return reseq,reind
    
### jump am seq ##
def jumpamseq(re,slow,stoep,icue,lseg,di,nosd):
      reseq,reind = noseqnew(stoep,re,lseg,slow[0],di,nosd)
      errori = zeros(2)
      for e in range(1,re-1):# errind[~((errind==49)|(errind==0))]:           #eliminate the last element if it is error
          if (int(reind[e-1]/49)!= int(reind[e]/49)):# and ((reind[e]-1)!=reind[e-1]):
               errori[0] += 1    # jump among seq
          else:
               errori[1] += 1    # no jump or jump in seq
      return errori

def jumpinseq(re,slow,stoep,icue,lseg,di,nosd):
      reseq,reind = noseqnew(stoep,re,lseg,slow[0],di,nosd)
      errori = zeros(2)
      for e in range(1,re-1):# errind[~((errind==49)|(errind==0))]:           #eliminate the last element if it is error
          if (int(reind[e-1]/49)== int(reind[e]/49)) and ((reind[e]-1)!=reind[e-1]):
               errori[0] += 1    # jump within seq
          else:
               errori[1] += 1    # no jump or jump out of seq
      return errori


#  break each episode respectively then put them into the pool #
#  so that there is no conection between two episodes          #
def segpool(pooldata,ls,nst):
    for i in range(nst):                                   
        segslow = genpiece(pooldata[i],ls)
        if i== 0:
            allep = segslow
        else:
            allep = concatenate((allep,segslow),0)
    return allep
    
    
def whiten(refseq):
    pcanode = mdp.nodes.PCANode()
    pcanode(refseq)    
    D = diag(pcanode.d**-.5)

    return pcanode,D
    
def lowdi(highdi,lendata,nbin,bo):
    lowout  = zeros((lendata,nbin*nbin))
    for t in range(lendata):
        pattern = reshape(highdi[t],(bo,bo))
        indx = nonzero(pattern)[0]
        indy = nonzero(pattern)[1]
        H,tlx,tly = histogram2d(indx,indy,nbin,[[0,bo],[0,bo]],normed= False)
        lowout[t] = reshape(H,(nbin*nbin))
    return lowout

def TvsL(nstorp,networka,networkb,nois,reflen):    
    storinp_T = zeros((nstorp,re,nbin**2))
    storpa_aT = zeros((nstorp,re,di)) 
    storpa_iT = zeros((nstorp,re,di))
    storinp_L = zeros((nstorp,re,nbin**2))
    storpa_aL = zeros((nstorp,re,di))
    storpa_iL = zeros((nstorp,re,di))
    
    refinp_T = lsau.Umore(reflen,p,q,omega2,bo,bar_l2,bar_l2,bar_w,no)
    reft = lowdi(refinp_T,reflen,nbin,bo)
    del refinp_T

    refinp_L = lsau.Arrmore(reflen,p,q,omega2,bo,bar_l2,bar_w,no)
    refl = lowdi(refinp_L,reflen,nbin,bo)
    del refinp_L

    ref_aT = networka(reft)[:,:4]
    ref_iT = networkb(reft)[:,:4]
    ref_aL = networkb(refl)[:,:4]
    ref_iL = networka(refl)[:,:4]    
    
    pca_aT,D_aT = whiten(ref_aT)
    pca_iT,D_iT = whiten(ref_iT)
    pca_aL,D_aL = whiten(ref_aL)
    pca_iL,D_iL = whiten(ref_iL)
    
    
    for i in range(nstorp):
        sT = lsau.Umore(re,p,q,omega2,bo,bar_l2,bar_l2,bar_w,no)

        storinp_T[i] = lowdi(sT,re,nbin,bo)
        sL = lsau.Arrmore(re,p,q,omega2,bo,bar_l2,bar_w,no)
        storinp_L[i] = lowdi(sL,re,nbin,bo)        

        saT = networka(storinp_T[i])[:,:4]
        siT = networkb(storinp_T[i])[:,:4]
        saL= networkb(storinp_L[i])[:,:4]
        siL = networka(storinp_L[i])[:,:4]

        storpa_aT[i] = inner(pca_aT(saT),D_aT)
        storpa_iT[i] = inner(pca_iT(siT),D_iT)
        storpa_aL[i] = inner(pca_aL(saL),D_aL)
        storpa_iL[i] = inner(pca_iL(siL),D_iL)

    
    data_aT = [storpa_aT[r] for r in range(nstorp)]
    data_iT = [storpa_iT[r] for r in range(nstorp)]
    
    data_aL = [storpa_aL[r] for r in range(nstorp)]
    data_iL = [storpa_iL[r] for r in range(nstorp)]    
    
    allsto_aT = segpool(data_aT,ls,nstorp)
    allsto_iT = segpool(data_iT,ls,nstorp)
    allsto_aL = segpool(data_aL,ls,nstorp)
    allsto_iL= segpool(data_iL,ls,nstorp)
    
    reerr = zeros((4,nstorp,re))
    drift = zeros((4,nstorp,re))
    jins = zeros((4,nstorp,2))
    jams = zeros((4,nstorp,2))

    for j in range(nstorp):
#### retrieval error ###
        reerr[0,j] = distancewithoutk(re,storpa_aT[j],allsto_aT,icue,ls,di,nois)#[0])   #normalized noise according to the randam level
        reerr[1,j] = distancewithoutk(re,storpa_iT[j],allsto_iT,icue,ls,di,nois)#[1])
        reerr[2,j] = distancewithoutk(re,storpa_aL[j],allsto_aL,icue,ls,di,nois)#[2])    
        reerr[3,j] = distancewithoutk(re,storpa_iL[j],allsto_iL,icue,ls,di,nois)#[3])
#### noise drift ###
        drift[0,j] = distdrift(re,storpa_aT[j],allsto_aT,icue,ls,di,nois)
        drift[1,j] = distdrift(re,storpa_iT[j],allsto_iT,icue,ls,di,nois)
        drift[2,j] = distdrift(re,storpa_aL[j],allsto_aL,icue,ls,di,nois)
        drift[3,j] = distdrift(re,storpa_iL[j],allsto_iL,icue,ls,di,nois)
#### in seq jump ####
        jins[0,j] = jumpinseq(re,storpa_aT[j],allsto_aT,icue,ls,di,nois)
        jins[1,j] = jumpinseq(re,storpa_iT[j],allsto_iT,icue,ls,di,nois)
        jins[2,j] = jumpinseq(re,storpa_aL[j],allsto_aL,icue,ls,di,nois)
        jins[3,j] = jumpinseq(re,storpa_iL[j],allsto_iL,icue,ls,di,nois)
#### am seq jump ####
        jams[0,j] = jumpamseq(re,storpa_aT[j],allsto_aT,icue,ls,di,nois)
        jams[1,j] = jumpamseq(re,storpa_iT[j],allsto_iT,icue,ls,di,nois)
        jams[2,j] = jumpamseq(re,storpa_aL[j],allsto_aL,icue,ls,di,nois)
        jams[3,j] = jumpamseq(re,storpa_iL[j],allsto_iL,icue,ls,di,nois)
        
    return reerr,drift,jins,jams


re = 50                                                # retrival length
icue = 0                                               # index of the cue
ls = 2                                                   # length segment
nosd  = 0.1                                        # std of the adding noise
di = 4


bo = 300
step = 45
bar_l1 = 24.5 # for shape H
#bar_l1 = 3.5
bar_l2 = 44.5# 4.5
bar_wh = 6.5
bar_v = 3.5
#bar_v = 6.5
bar_w = 14.5#6.5
no = 0
v = 1
sudturnp = 0
noro = 0


p = .02*pi#.02*pi
q = .06#.06
omega1 = .009*e   

omega2 =.025*e #.015*e#.018*e

trainstep = 10000 #10000
reflen = 3000
nois =[.0,.1,.2,.3,.4,.5,.6,.7,.8,.9] #.2#

nbin = 30
nstorp = 30
#nstorp = [500,1000,1500,2000] #30
outname = ['no','01','02','03','04','05','06','07','08','09']
#outname = ['500','1000','1500','2000']
outname = ['02']
for t in range(len(outname)):
  atraindata = lsau.Umore(trainstep,p,q,omega2,bo,bar_l2,bar_l2,bar_w,no)
  traina = lowdi(atraindata,trainstep,nbin,bo)
  del atraindata

  btraindata = lsau.Arrmore(trainstep,p,q,omega2,bo,bar_l2,bar_w,no)
  trainb = lowdi(btraindata,trainstep,nbin,bo)
  del btraindata

  networka = hierachynet(30,15,3)
  networkb = hierachynet(30,15,3)
  networka(traina)
  networkb(trainb)

  reerr,drift,jins,jams= TvsL(nstorp,networka,networkb,nois[t],reflen)
  err = mean(reerr,1)
  dri = mean(drift,1)
  jin = mean(jins,1)
  jam = mean(jams,1)

  savez("UA_30_opp_"+outname[t]+"_%d.npz"%(myid),err = err,dri=dri,jins=jin,jams=jam)

pypar.finalize()




