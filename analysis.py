# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 14:39:52 2014

@author: jing
"""
###############################################################################################
# This code generate two object comparing data which has no need of retrieval process         #
# including pairwise distance distribution (single seq and all seq),mean concecutive distance,#
# distance from each element to first element, sequence distance divergence,mean random       #
# retrieval error                                                                             #
###############################################################################################
from numpy import *
from pylab import *
import mdp
import scipy as sy
import inpobject as lsau
def hierachynettest(bod,recf,ovl): 
    switchboard = mdp.hinet.Rectangular2dSwitchboard(in_channels_xy = (bod,bod),
                                                 field_channels_xy=(recf,recf),
                                                 field_spacing_xy =(ovl,ovl))

    sfa_dim = 48
    sfa_lower_out= 32

    sfanode = mdp.nodes.SFANode(input_dim = switchboard.out_channel_dim, output_dim = sfa_dim)
    sfa2node = mdp.nodes.QuadraticExpansionNode(input_dim=sfa_dim)
    sfanode2 = mdp.nodes.SFANode(input_dim = sfa2node.output_dim,output_dim = sfa_lower_out)
    flownode = mdp.hinet.FlowNode(mdp.Flow([sfanode,sfa2node,sfanode2]))
    sfalayer = mdp.hinet.CloneLayer(flownode, n_nodes = switchboard.output_channels)
    flow = mdp.Flow([switchboard, sfalayer])

    sfa_upper_bo = 6
    sfa_upper_recf = 4
    sfa_upper_out = 32 
    sfa_top_out = 10
    ovl2 = 2  
    switchboard2 =  mdp.hinet.Rectangular2dSwitchboard(in_channels_xy = (sfa_upper_bo,sfa_upper_bo),
                                                 field_channels_xy=(sfa_upper_recf,sfa_upper_recf),
                                                 field_spacing_xy = (ovl2,ovl2),  
                                                 in_channel_dim = sfa_lower_out)

    sfa_uppernode = mdp.nodes.SFANode(input_dim = switchboard2.out_channel_dim, output_dim = sfa_dim)
    sfa_upperexp = mdp.nodes.QuadraticExpansionNode(input_dim = sfa_dim)
    sfa_uppernode2 = mdp.nodes.SFANode(input_dim = sfa_upperexp.output_dim, output_dim = sfa_upper_out)
    upper_flownode = mdp.hinet.FlowNode(mdp.Flow([sfa_uppernode,sfa_upperexp,sfa_uppernode2]))
    upper_sfalayer = mdp.hinet.CloneLayer(upper_flownode, n_nodes = switchboard2.output_channels)
    sfa_top_node = mdp.nodes.SFANode(input_dim = upper_sfalayer.output_dim, output_dim = sfa_dim)
    sfa_topexp =mdp.nodes.QuadraticExpansionNode(input_dim = sfa_dim)
    sfa_topnode2 =mdp.nodes.SFANode(input_dim = sfa_topexp.output_dim,output_dim = sfa_top_out)
    sfa_over_node = mdp.hinet.FlowNode(mdp.Flow([sfa_top_node,sfa_topexp,sfa_topnode2]))

    network = mdp.Flow([switchboard,sfalayer,switchboard2,upper_sfalayer,sfa_over_node])
    return network
                                                                                             

  
def noseqN2(thr,stoep,lseq,lseg,q,di,sd):

    reseq = zeros((lseq,di))
    p = q+normal(0,(sd+1e-20),di)

    for j in range(lseq-1):

        p=q+normal(0,(sd+1e-20),di)
        disma = sqrt(sum((stoep[:,:di]-p)**2,1))
        mind = where(disma == min(disma))[0]
        shuffle(mind)
        ind = mind[0]

        reseq[j] = stoep[ind,:di]
        q = stoep[ind,lseg*di-di::]

    reseq[lseq-1] = stoep[ind,lseg*di-di::]
    return reseq
# generate pieces #
def genpiece(slow,lseg):
    seg = mdp.nodes.TimeFramesNode(lseg)                 # only work for 2D data
    stoep = seg(slow)

    return stoep

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

def calforseq(storseq):
    #distance to first element
    disttofirst = mean(sqrt(sum((storseq - storseq[:,0,newaxis])**2,2)),0)
    #concecutive distance
    condist = mean(sqrt(sum((storseq[:,1::,:]-storseq[:,:49,:])**2,2)))
    #pairwise distance for all storing ep
    z = shape(storseq)
    seqall = reshape(storseq,(z[0]*z[1],z[2]))
    distall = sy.spatial.distance.pdist(seqall,'euclidean')
    return disttofirst, condist,distall
## random level 
def distancerand(re,nstorp,slow,storseq):
     
    rinds = randint(0,nstorp,re)
    rinde = randint(0,re,re)
    reseq = storseq[rinds,rinde]
    l = mean(sqrt(sum((slow-reseq)**2,1)))
    return l
    
def compdist(nstorp,ns,ne,seq):
    indall = range(nstorp)
    indall.remove(ns)
    disat = sum((seq[indall,:48]- seq[ns,ne])**2,2) 
    cloind = unravel_index(disat.argmin(),disat.shape)
    lentest = min(49 - cloind[1],10,49-ne)
    seqdist = sum((seq[ns,ne:(ne+lentest)]- seq[cloind[0],cloind[1]:(cloind[1]+lentest)])**2,1)
    seqdist = seqdist - seqdist[0]
    return lentest,seqdist

        
def allresult(nstorp,networka,networkb,nois,reflen):
    storinp_T = zeros((nstorp,re,nbin**2))
    storpa_aT = zeros((nstorp,re,di))
    storpa_iT = zeros((nstorp,re,di))
    storinp_L = zeros((nstorp,re,nbin**2))
    storpa_aL = zeros((nstorp,re,di))
    storpa_iL = zeros((nstorp,re,di))

    refinp_T = lsau.Umore(reflen,p,q,omega2,bo,bar_l2,bar_l2,bar_w,no)
    reft = lowdi(refinp_T,reflen,nbin,bo)
    refinp_L = lsau.Arrmore(reflen,p,q,omega2,bo,bar_l2,bar_w,no)
    refl = lowdi(refinp_L,reflen,nbin,bo)
    ref_aT = networka(reft)[:,:di]
    ref_iT = networkb(reft)[:,:di]
    ref_aL = networkb(refl)[:,:di]
    ref_iL = networka(refl)[:,:di]

    pca_aT,D_aT = whiten(ref_aT)
    pca_iT,D_iT = whiten(ref_iT)
    pca_aL,D_aL = whiten(ref_aL)
    pca_iL,D_iL = whiten(ref_iL)
    
    dist_inseq = zeros((4,nstorp,re*(re-1)/2))
    dist_allseq = zeros((4,nstorp*re*(nstorp*re -1)/2))
    dist_first = zeros((4,re))
    dist_con = zeros(4)
    rl = zeros((4,nstorp))
    nlen = zeros((4,10))
    ctl = zeros((4,10))
    for i in range(nstorp):

        sT = lsau.Umore(re,p,q,omega2,bo,bar_l2,bar_l2,bar_w,no)
        storinp_T[i] = lowdi(sT,re,nbin,bo)
        sL = lsau.Arrmore(re,p,q,omega2,bo,bar_l2,bar_w,no)

        storinp_L[i] = lowdi(sL,re,nbin,bo)
        saT = networka(storinp_T[i])[:,:di]
        siT = networkb(storinp_T[i])[:,:di]
        saL= networkb(storinp_L[i])[:,:di]
        siL = networka(storinp_L[i])[:,:di]

        storpa_aT[i] = inner(pca_aT(saT),D_aT)
        storpa_iT[i] = inner(pca_iT(siT),D_iT)
        storpa_aL[i] = inner(pca_aL(saL),D_aL)
        storpa_iL[i] = inner(pca_iL(siL),D_iL)

        dist_inseq[0,i] = sy.spatial.distance.pdist(storpa_aT[i],'euclidean')       #largest distance 
        dist_inseq[1,i] = sy.spatial.distance.pdist(storpa_iT[i],'euclidean')
        dist_inseq[2,i] = sy.spatial.distance.pdist(storpa_aL[i],'euclidean')
        dist_inseq[3,i] = sy.spatial.distance.pdist(storpa_iL[i],'euclidean')
        rl[0,i] = distancerand(re,nstorp,storpa_aT[i],storpa_aT)
        rl[1,i] = distancerand(re,nstorp,storpa_iT[i],storpa_iT)
        rl[2,i] = distancerand(re,nstorp,storpa_aL[i],storpa_aL)
        rl[3,i] = distancerand(re,nstorp,storpa_iL[i],storpa_iL)  

    for j in range(210):

        ne = randint(0,(re-2),4)   
        nlen = zeros((4,10))
        ctl = zeros((4,10))
        c1,aT = compdist(nstorp,j%30,ne[0],storpa_aT)
        c2,iT = compdist(nstorp,j%30,ne[1],storpa_iT)
        c3,aL = compdist(nstorp,j%30,ne[2],storpa_aL)
        c4,iL = compdist(nstorp,j%30,ne[3],storpa_iL)
        nlen[0,:c1] += ones(c1)
        nlen[1,:c2] += ones(c2)
        nlen[2,:c3] += ones(c3)
        nlen[3,:c4] += ones(c4)
        ctl[0,:c1] += aT
        ctl[1,:c2] += iT
        ctl[2,:c3] += aL
        ctl[3,:c4] += iL
        nlen+=1e-20
        ctl+=1e-20
    mdist = ctl/nlen  
        
    dist_first[0],dist_con[0],dist_allseq[0] = calforseq(storpa_aT)
    dist_first[1],dist_con[1],dist_allseq[1] = calforseq(storpa_iT)
    dist_first[2],dist_con[2],dist_allseq[2] = calforseq(storpa_aL)
    dist_first[3],dist_con[3],dist_allseq[3] = calforseq(storpa_iL)
    dist_inseq_all = reshape(dist_inseq,(4,size(dist_inseq[0])))
    rl_all = mean(rl,1)
    
    return dist_inseq_all,dist_first,dist_con,dist_allseq,rl_all,mdist

def lowdi(highdi,lendata,nbin,bo):
    lowout  = zeros((lendata,nbin*nbin))
    for t in range(lendata):
        pattern = reshape(highdi[t],(bo,bo))
        indx = nonzero(pattern)[0]
        indy = nonzero(pattern)[1]
        H,tlx,tly = histogram2d(indx,indy,nbin,[[0,bo],[0,bo]],normed= False)
        lowout[t] = reshape(H,(nbin*nbin))
    return lowout

v = 1
re = 50                                                # retrival length
icue = 0                                               # index of the cue
ls = 2                                                   # length segment
nosd  = 0.1                                        # std of the adding noise
di = 4


bo = 300
step = 45
bar_l1 = 24.5
bar_l2 = 44.5
bar_l3 = 3.5
bar_v = 3.5
bar_w = 14.5#6.5
bar_wh = 6.5
no = 0
v = 1
sudturnp = 0
noro = 0


p = .02*pi #.02*pi only for HX
q = .06 #.06 only for HX
omega1 = .009*e #.0197 correct for T
#omega1 = .032*e

omega2 =.025*e #.015*e#.018*e
#omega2 = .0197*e  #.032 correct for L

nstorp = 30 # number of storing episodes
Nr = 60
nois = .0
trainstep = 10000 #10000
reflen = 3000

dist_first = zeros((Nr,4,re))
dist_con = zeros((Nr,4))
ranl = zeros((Nr,4))
divseq = zeros((Nr,4,10))
nbin = 30
for t in range(Nr):
    
    atraindata = lsau.Umore(trainstep,p,q,omega2,bo,bar_l2,bar_l2,bar_w,no)
    traina = lowdi(atraindata,trainstep,nbin,bo)
    del atraindata
    btraindata = lsau.Arrmore(trainstep,p,q,omega2,bo,bar_l2,bar_w,no)
    trainb = lowdi(btraindata,trainstep,nbin,bo)
    del btraindata
    networka = hierachynettest(30,15,3)
    networkb = hierachynettest(30,15,3)
    networka(traina)
    networkb(trainb)

    dist_inseq_all,dist_first[t],dist_con[t],dist_allseq,ranl[t],divseq[t]= allresult(nstorp,networka,networkb,nois,reflen)
dist_first_all = mean(dist_first,0) 
dist_con_all = mean(dist_con,0)    
rl_all = mean(ranl,0)
divseq_all = mean(divseq,0)

np.savez("UA_alldata_opp.npz",distinseq = dist_inseq_all,distallseq = dist_allseq,distfall = dist_first_all,distcon = dist_con_all,rlall = 
rl_all,divall = divseq_all)


