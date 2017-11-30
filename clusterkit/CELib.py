#!/usr/bin/env python
##########################################################################
#Nanokit module for Python
#Last editted by Dr. Bin Ouyang 2017.02.17
##########################################################################
#This is the lib containing all the methods for cluster expansion
#******dismatcreate: Creating the distance matrix for each pari of atoms
#******clustercount: Count the number for each cluster
#******clustergen: Generate random cluster at specific condition
#******clusterfit: Fit for the cluster coefficient on base of input
##########################################################################
#Format of definition
#functions are all non-capitalized
#dictionaries are all capitalized
#Other variables are first letter capitalized
#Please note that for keys in some dictionary, if the name is the same
#with INCAR tags, it will be all capital

__author__ = 'Bin Ouyang'
import Vasp as vp
import numpy as np
import math
import random as rd
import copy
import MathKit


#######################
def clustercount1(Clusterdes,POS,TS=0.2):
    '''
    enumerate and count clusters in a given lattce, this version is cleaner
    and more robust than the method below: clustercount

    Args:
        Clusterdes: Cluster description, in the format of list, somthing
                    like [[[0,1,2],[2.6,2.7,2.8]],[[1,1],[2.5]],[[2]]]
        POS: Dictionary containing the position information, in the format of POSCAR
        TS: Allowed variation of cluster bond length
        Outputs: ClusterLst, which is a list with all the description of indentified
                 clusters as specified in Clusterdes
    '''
    ClusterNum = len(Clusterdes);
    ClusterLst = [[] for i in range(ClusterNum)];
    SumLst = [sum(POS['AtomNum'][0:i]) for i in range(POS['EleNum'])];

    for CInd, Cluster in enumerate(Clusterdes):
        CSize = len(Cluster[0]); #Cluster Size
        if CSize == 1:
            for i in range(POS['AtomNum'][Cluster[0][0]]):
                GInd = i + sum(POS['AtomNum'][0:Cluster[0][0]]);
                ClusterLst[CInd].append([GInd]);
        else:
            IndLst = [0]*CSize; IndLstMax = [];
            GIndLst = [0]*CSize;
            GrpLst = MathKit.findGrpLst(Cluster[0]);
            PermuteLst = MathKit.grporderPermute(Cluster[0],GrpLst);
            DistRef = MathKit.grpSort(Cluster[1],PermuteLst);
            for Ele in Cluster[0]:
                IndLstMax.append(POS['AtomNum'][Ele] - 1);
            while (IndLst[-1] <= IndLstMax[-1]):
                for i, Ind in enumerate(IndLst):
                    GIndLst[i] = Ind + SumLst[Cluster[0][i]];
                Dist = [];
                for i in range(CSize):
                    for j in range(i+1,CSize):
                        Dist.append(POS['dismat'][GIndLst[i]][GIndLst[j]]);
                flag = 1;
                Dist = MathKit.grpSort(Dist,PermuteLst);
                for Dind,D in enumerate(Dist):
                    if abs(D - DistRef[Dind]) > TS:
                        flag = 0;
                        break;
                if flag:
                    ClusterLst[CInd].append(list(GIndLst));
                MathKit.lstGrpAdd(IndLst,IndLstMax,GrpLst);
    return ClusterLst;

#######################
def clustercount(Clusterdes,POS):
    '''
    enumerate and count clusters in a given lattce

    Args:
        Clusterdes: Cluster description, in the format of list, somthing
                    like [[[0,1,2],[2.6,2.7,2.8]],[[1,1],[2.5]],[[2]]]
        POS: Dictionary containing the position information, in the format of POSCAR
        Outputs: Clusterlst, which is a list with all the description of indentified
                 clusters as specified in Clusterdes
    '''
    TS1 = 0.9; TS2 = 1.1; #The range of bond length variation
    ClusterNum = len(Clusterdes);
    ClusterLst = [[] for i in range(ClusterNum)]; 
    for ClusterInd, Cluster in enumerate(Clusterdes): 
        ClusterSize = len(Cluster[0]); 
        for Ind0 in range(POS['AtomNum'][Cluster[0][0]]):
            if (ClusterSize == 1):
                GInd0 = Ind0 + sum(POS['AtomNum'][0:Cluster[0][0]]);
                ClusterLst[ClusterInd].append([GInd0]);
            elif (ClusterSize == 2):
                for Ind1 in range(POS['AtomNum'][Cluster[0][1]]):
                    GInd0 = Ind0 + sum(POS['AtomNum'][0:Cluster[0][0]]);
                    GInd1 = Ind1 + sum(POS['AtomNum'][0:Cluster[0][1]]);
                    Dis_01 = POS['dismat'][GInd0][GInd1];
                    if (Dis_01*TS1<Cluster[1][0]) & (Dis_01*TS2>Cluster[1][0]):
                        #print 'yes';
                        noexist = 1;
                        for i in range(len(ClusterLst[ClusterInd])):
                            if (sorted([GInd0,GInd1])==sorted(ClusterLst[ClusterInd][i])):
                                noexist = 0;
                                #print GInd0, GInd1, ClusterLst[ClusterInd][i]
                        if noexist:
                            ClusterLst[ClusterInd].append([GInd0,GInd1]);
                       # Count = Count + 1;
            elif (ClusterSize ==3):
                for Ind1 in range(POS['AtomNum'][Cluster[0][1]]):
                    for Ind2 in range(POS['AtomNum'][Cluster[0][2]]):
                        GInd0 = Ind0 + sum(POS['AtomNum'][0:Cluster[0][0]]);
                        GInd1 = Ind1 + sum(POS['AtomNum'][0:Cluster[0][1]]);
                        GInd2 = Ind2 + sum(POS['AtomNum'][0:Cluster[0][2]]);
                        Dis_01 = POS['dismat'][GInd0][GInd1];
                        Dis_02 = POS['dismat'][GInd0][GInd2];
                        Dis_12 = POS['dismat'][GInd1][GInd2];
                        #print Dis_01, Dis_02, Dis_03
                        flag_01 = (Dis_01*TS1<Cluster[1][0]) & \
                                (Dis_01*TS2>Cluster[1][0]);
                        flag_02 = (Dis_02*TS1<Cluster[1][1]) & \
                                (Dis_02*TS2>Cluster[1][1]);
                        flag_12 = (Dis_12*TS1<Cluster[1][2]) & \
                                (Dis_12*TS2>Cluster[1][2]);
                        if (flag_01*flag_02*flag_12):
                            noexist = 1;
                            #print "cluster found"
                            for i in range(len(ClusterLst[ClusterInd])):
                                if (sorted([GInd0,GInd1,GInd2])==\
                                        sorted(ClusterLst[ClusterInd][i])):
                                    noexist = 0;
                            if noexist:
                                ClusterLst[ClusterInd].append([GInd0,GInd1,GInd2]);
                            #Count = Count + 1;
            elif (ClusterSize == 4):
                for Ind1 in range(POS['AtomNum'][Cluster[0][1]]):
                    for Ind2 in range(POS['AtomNum'][Cluster[0][2]]):
                        for Ind3 in range(POS['AtomNum'][Cluster[0][3]]):
                            GInd0 = Ind0 + sum(POS['AtomNum'][0:Cluster[0][0]]);
                            GInd1 = Ind1 + sum(POS['AtomNum'][0:Cluster[0][1]]);
                            GInd2 = Ind2 + sum(POS['AtomNum'][0:Cluster[0][2]]);
                            GInd3 = Ind3 + sum(POS['AtomNum'][0:Cluster[0][3]]);
                            Dis_01 = POS['dismat'][GInd0][GInd1];
                            Dis_02 = POS['dismat'][GInd0][GInd2];
                            Dis_03 = POS['dismat'][GInd0][GInd3];
                            Dis_12 = POS['dismat'][GInd1][GInd2];
                            Dis_13 = POS['dismat'][GInd1][GInd3];
                            Dis_23 = POS['dismat'][GInd2][GInd3];
                            flag_01 = (Dis_01*TS1<Cluster[1][0]) & \
                                    (Dis_01*TS2>Cluster[1][0]);
                            flag_02 = (Dis_02*TS1<Cluster[1][1]) & \
                                    (Dis_02*TS2>Cluster[1][1]);
                            flag_03 = (Dis_03*TS1<Cluster[1][2]) & \
                                    (Dis_03*TS2>Cluster[1][2]);
                            flag_12 = (Dis_12*TS1<Cluster[1][3]) & \
                                    (Dis_12*TS2>Cluster[1][3]);
                            flag_13 = (Dis_13*TS1<Cluster[1][4]) & \
                                    (Dis_13*TS2>Cluster[1][4]);
                            flag_23 = (Dis_23*TS1<Cluster[1][5]) & \
                                    (Dis_23*TS2>Cluster[1][5]);
                            if (flag_01*flag_02*flag_03*flag_12*flag_13*flag_23):
                                noexist = 1;
                                #print "cluster found"
                                for i in range(len(ClusterLst[ClusterInd])):
                                    if (sorted([GInd0,GInd1,GInd2,GInd3])==\
                                            sorted(ClusterLst[ClusterInd][i])):
                                        noexist = 0;
                                        #print('yes')
                                if noexist:
                                    ClusterLst[ClusterInd].append([GInd0,GInd1,GInd2,GInd3]);
                                #Count = Count + 1;
            else:
                print "Okay, I don not know what to do with this cluster \
                        size "+str(CluserSize);
    return ClusterLst
#######################

#######################
def clusterswap1(ClusterDes,POS,ClusterLst,Atom1,Atom2,Ind1,Ind2,TS=0.2):
    '''
    Update the cluster information after swapping atoms
    This is a cleaner and more robust version of clusterswap method below

    Args:
        ClusterDes: Description about clusters
        POS: POSCAR dictionary
        ClusterLst: Cluster information
        Atom1, Atom2: Atom sublattice
        Ind1, Ind2: Atom indices
    '''
    ClusterNum = len(ClusterLst);
    SumLst = [sum(POS['AtomNum'][0:i]) for i in range(POS['EleNum'])];

    for LstInd, Lst in enumerate(ClusterLst):
        for Ind, AtomInd in enumerate(Lst):
            if (Ind1 in AtomInd) | (Ind2 in AtomInd):
                ClusterLst[LstInd].remove(AtomInd);

    for CInd, Cluster in enumerate(ClusterDes):
        CSize = len(Cluster[0]);
        if (CSize == 1) & (Atom1 == Cluster[0][0]):
            ClusterLst[ClusterInd].append([Ind1]);
        elif (CSize == 1) & (Atom2 == Cluster[0][0]):
            ClusterLst[ClusterInd].append([Ind2]);
        else:
            for AtomI, Atom in enumerate([Atom1, Atom2]):
                if Atom in Cluster[0]:
                    AtomInd = [Ind1,Ind2][AtomI];
                    AtomLoc = Cluster[0].index(Atom);
                    IndLst = [0]*(CSize-1); IndLstMax = [];
                    GIndLst = [0]*(CSize-1);
                    GrpLst = MathKit.findGrpLst(Cluster[0]);
                    PermuteLst = MathKit.grporderPermute(Cluster[0],GrpLst);
                    DistRef = MathKit.grpSort(Cluster[1],PermuteLst);
                    ClusterTmp = copy.deepcopy(Cluster[0]);
                    ClusterTmp.remove(Atom);
                    for Ele in ClusterTmp:
                        IndLstMax.append(POS['AtomNum'][Ele] - 1);
                    while (IndLst[-1] <= IndLstMax[-1]):
                        for i, Ind in enumerate(IndLst):
                            GIndLst[i] = Ind + SumLst[Cluster[0][i]];
                        GIndLst.insert(AtomLoc,AtomInd);
                        GIndLst = MathKit.grpSort(GIndLst,GrpLst)
                        Dist = [];
                        for i in range(CSize):
                            for j in range(i+1,CSize):
                                Dist.append(POS['dismat'][GIndLst[i]][GIndLst[j]]);
                        flag = 1;
                        Dist = MathKit.grpSort(Dist,PermuteLst);
                        for Dind, D in enumerate(Dist):
                            if abs (D - DistRef[ind]) > TS:
                                flag = 0;
                                break;
                        if flag:
                            ClusterLst[CInd].append(list(GIndLst));
                        MathKit.lstGrpAdd(IndLst,IndLstMax,GrpLst);
    return ClusterLst;
#######################

#######################
def clusterswap(ClusterDes,POS,ClusterLst,Atom1,Atom2,Ind1,Ind2):
    '''update the cluster information after swapping atom positionA'''
    #print Atom1, Atom2, Ind1, Ind2;
    #TS1 = 0.9; TS2 = 1.1; #The range of bond length variation
    TS0 = 0.5; #So the length variation of two atoms is within 0.5 A
    ClusterLstNew = copy.deepcopy(ClusterLst);
    ClusterNum = len(ClusterLst); #Number of clusters considerred

    AtomSumTab = [sum(POS['AtomNum'][0:i]) for i in range(POS['EleNum'])];

    for LstInd, Lst in enumerate(ClusterLst):
        for Ind, AtomInd in enumerate(Lst):
            if (Ind1 in AtomInd) | (Ind2 in AtomInd):
                ClusterLstNew[LstInd].remove(AtomInd);
    ClusterLst = copy.deepcopy(ClusterLstNew);
    for ClusterInd, Cluster in enumerate(ClusterDes):
        ClusterSize = len(Cluster[0]);
        if ((ClusterSize == 1)&(Atom1 == Cluster[0][0])):
            ClusterLst[ClusterInd].append([Ind1]);
        elif ((ClusterSize == 1)&(Atom2 == Cluster[0][0])):
            ClusterLst[ClusterInd].append([Ind2]);
        elif (ClusterSize == 2):
            for AtomI, Atom in enumerate([Atom1, Atom2]):
                if Atom in Cluster[0]:
                    OccLst = [i for i,val in enumerate(Cluster[0]) if val==Atom];
                    IndSite0 = [Ind1,Ind2][AtomI] - AtomSumTab[Atom];
                    SiteOrd = list(Cluster[0]);
                    Sites = list(SiteOrd); Sites.remove(Atom);
                    Site1 = Sites[0];
                    IndOrd = ['IndS1']; IndOrd.insert(OccLst[0],'IndSite0');
                    for IndS1 in range(POS['AtomNum'][Site1]):
                        GIndS0 = eval(IndOrd[0]) + AtomSumTab[SiteOrd[0]];
                        GIndS1 = eval(IndOrd[1]) + AtomSumTab[SiteOrd[1]];
                        Dis_01 = POS['dismat'][GIndS0][GIndS1];
                        if (Dis_01-TS0<Cluster[1][0]) & (Dis_01+TS0>Cluster[1][0]):
                            noexist = 1;
                            for i in range(len(ClusterLst[ClusterInd])):
                                if (sorted([GIndS0,GIndS1])==\
                                        sorted(ClusterLst[ClusterInd][i])):
                                    noexist = 0;
                            if noexist:
                                ClusterLst[ClusterInd].append([GIndS0,GIndS1]);
        elif (ClusterSize == 3):
            IsCount = 0; IsSort = 0;
            for AtomI, Atom in enumerate([Atom1, Atom2]):
                if Atom in Cluster[0]:
                    OccLst = [i for i,val in enumerate(Cluster[0]) if val==Atom];
                    if len(OccLst) > 1:
                        IsSort = '';
                        for Occ in range(3):
                            if Occ in OccLst:
                                IsSort+='1';
                            else:
                                IsSort+='2';
                        IsSort = int(IsSort);
                    IndSite0 = [Ind1,Ind2][AtomI] - AtomSumTab[Atom];
                    SiteOrd = list(Cluster[0]); 
                    Sites = list(SiteOrd); Sites.remove(Atom);
                    Site1, Site2 = Sites[0], Sites[1];
                    IndOrd = ['IndS1','IndS2']; IndOrd.insert(OccLst[0],'IndSite0');
                    IsCount = 1;
                if (IsCount == 1):
                    for IndS1 in range(POS['AtomNum'][Site1]):
                        for IndS2 in range(POS['AtomNum'][Site2]):
                            GIndS0 = eval(IndOrd[0]) + AtomSumTab[SiteOrd[0]];
                            GIndS1 = eval(IndOrd[1]) + AtomSumTab[SiteOrd[1]];
                            GIndS2 = eval(IndOrd[2]) + AtomSumTab[SiteOrd[2]];
                            Dis_01 = POS['dismat'][GIndS0][GIndS1];
                            Dis_02 = POS['dismat'][GIndS0][GIndS2];
                            Dis_12 = POS['dismat'][GIndS1][GIndS2];
                            if (IsSort == 0):
                                Dis=[Dis_01,Dis_02,Dis_12];
                                ClusterDis = [Cluster[1][0],Cluster[1][1],\
                                        Cluster[1][2]];
                            elif (IsSort == 111):
                                Dis=sorted([Dis_01,Dis_02,Dis_12]);
                                ClusterDis = sorted([Cluster[1][0],Cluster[1][1],\
                                        Cluster[1][2]]);
                            elif (IsSort == 112):
                                Dis = sorted([Dis_02,Dis_12]);
                                Dis.insert(0,Dis_01);
                                ClusterDis = sorted([Cluster[1][1],Cluster[1][2]]);
                                ClusterDis.insert(0,Cluster[1][0]);
                            elif (IsSort == 121):
                                Dis = sorted([Dis_01,Dis_12]);
                                Dis.insert(1,Dis_02);
                                ClusterDis = sorted([Cluster[1][0],Cluster[1][2]]);
                                ClusterDis.insert(1,Cluster[1][1]);
                            elif (IsSort == 211):
                                Dis = sorted([Dis_01,Dis_02]);
                                Dis.append(Dis_12);
                                ClusterDis = sorted([Cluster[1][0],Cluster[1][1]]);
                                ClusterDis.append(Cluster[1][2]);
                            flag_01 = (Dis[0]-TS0<ClusterDis[0]) & \
                                    (Dis[0]+TS0>ClusterDis[0]);
                            flag_02 = (Dis[1]-TS0<ClusterDis[1]) & \
                                    (Dis[1]+TS0>ClusterDis[1]);
                            flag_12 = (Dis[2]-TS0<ClusterDis[2]) & \
                                    (Dis[2]+TS0>ClusterDis[2]);
                            if (sum([flag_01,flag_02,flag_12]) == 3):
                                noexist = 1;
                                for i in range(len(ClusterLst[ClusterInd])):
                                    if (sorted([GIndS0,GIndS1,GIndS2])==\
                                            sorted(ClusterLst[ClusterInd][i])):
                                        noexist = 0;
                                if noexist:
                                    ClusterLst[ClusterInd].append([GIndS0,GIndS1,GIndS2]);
        elif (ClusterSize == 4):
            print "We cannot deal with cluster size large than 4!!!!!!"
    return ClusterLst;
#######################

#######################
def clusterinsert(ClusterDes,ClusterLst,Atom1,Atom2,Ind1,Ind2):
    '''
    Update the cluster information after insert atom position
    Not completed yet!
    '''
    if (Ind1 > Ind2):
        for i in range(Ind2,Ind1):
            clusterswap(dismat,i,Ind1);
    elif (Ind1 < Ind2):
        for i in range(Ind2,Ind1,-1):
            clusterswap(dismat,i,Ind1);
    return ClusterLst;
#######################

#######################
#def clustergen(ClusterDes,POSBase,AtomNum):
#Still under development!!!!!!!!!!!!!!!!!!!!
#    '''Generate lattice containing sepcific cluster:
#    1. ClusterDes would be the same type of data as the same in cluster count so I take
#       the same example as:
#       [[[0,1,2],[2.6,2.7,2.8]],[[1,1],[2.5]],[[2]]];
#       it tells the function that it needs two types of cluster, one consist of three
#       atoms type 0, 1, 2 with distance in the d_01=2.6, d_02=2.7, d_12=2.8, the
#       second cluster consist of two atoms of type 1, with the distance d_01=2.5,
#       while the third cluster consit of one atom of type 2.
#    2. The basic conception is to predefine the lattice and allocate cluster needed.
#    3. It should be noted that usually we only predefine the clusters that can be missed
#       during randomizing the structure since they have less probability to appear, too
#       much predefined structures would add too much constraint to the lattice therefore
#       it sometimes will not work.
#    4. The POSBase dictionary should take a different format from regular POS.It should 
#       be read from a POSCAR file for sure, but the poscar file should only distinguish 
#       different sublattice. Therefore the generation will decide the arrangement of >=1 
#       types of atoms in one lattice.
#    5. The AtomNum will give the number of elements and atoms we need for each atomic type
#       For example, [[16,16],[16,16],[31,1]] tells that we need two types of atoms for each
#       sublattice while for the 1st and second sublattice there will be 16 atoms each. The
#       third sublattice will have 31 and 1 of two types of atoms.
#    the dismat tag'''
#
#    TS1 = 0.9; TS2 = 1.1; #The range of bond length variation
#    ClusterNum = len(Clusterdes);
#    ClusterLst = [[] for i in range(ClusterNum)];
#    for ClusterInd, Cluster in enumerate(Clusterdes):
#        # So hopefully the cluster size will not be larger than 4
#        ClusterSize = len(Cluster[0]);
#        if (ClusterSize == 1):
#            print "One body cluster does not need to be predefined";
#        elif (ClusterSize == 2):
#            NoCluster = 1 #The flag telling whether we found the cluster we need
#            Type1 = Cluster[0][0]; Type2 = Cluster[0][1];
#            while (NoCluster):
#                Ind0 = rd.randint(0,POS['AtomNum'][Type1]]);
#                GInd0 = Ind0 + sum(POS['AtomNum'][0:Type1]]);
#                '''Now we need to determine the number of 1st NN of type2 for 
#                   selected atom'''
#                for Ind1 in range(POS['AtomNum'][Type2]):
#
#        elif (ClusterSize ==3):
#        else:
#            print "You know we do not consider clusters larger than 4";
#
#
#    return POS
#######################

#######################
def clusterRandom(POS,AtomNum,AtomName):
    '''The POSBase dictionary should take a different format from regular POS.It should
       be read from a POSCAR file for sure, but the poscar file should only distinguish
       different sublattice. Therefore the generation will decide the arrangement of >=1
       types of atoms in one lattice.
       The AtomNum will give the number of elements and atoms we need for each atomic type
       For example, [[16,16],[16,16],[31,1]] tells that we need two types of atoms for each
       sublattice while for the 1st and second sublattice there will be 16 atoms each. The
       third sublattice will have 31 and 1 of two types of atoms.'''
    NewAtomNum = []; NewAtomInd = []; NewLattPnt = [];
    for i, NumLst in enumerate(AtomNum):
        AtomPre = sum(POS['AtomNum'][0:i]);
        AtomLst = range(AtomPre,AtomPre+POS['AtomNum'][i]);
        for ind,Num in enumerate(NumLst):
            AtomLstPart = rd.sample(AtomLst,Num);
            NewAtomNum.append(Num); NewAtomInd.append(AtomLstPart);
            for Atom in AtomLstPart:
                AtomLst.remove(Atom);
    #print NewAtomNum, NewAtomInd
    for i, AtomLst in enumerate(NewAtomInd):
        for j in range(len(AtomLst)):
            #print i, NewAtomNum, AtomNum
            ind = AtomLst[j];
            #print ind;
            NewLattPnt.append(POS['LattPnt'][ind]);
    POS['AtomNum'] = NewAtomNum[:]; POS['LattPnt'] = NewLattPnt[:];
    POS['EleName'] = AtomName[:]; POS['AtomSum'] = sum(POS['AtomNum']);
    POS['EleNum'] = len(POS['AtomNum']);
    #print POS['AtomNum'], POS['AtomSum']
    # Allocate sublattice just incase
    flag = 0;
    for i in range(len(POS['AtomNum'])):
        for j in range(POS['AtomNum'][i]):
            for k in range(3):
                #print i,j,k,flag
                POS['SubLatt'][i][j].append(POS['LattPnt'][flag][k]);
                #print POS['SubLatt']
            flag = flag + 1;
            if(j!=POS['AtomNum'][i]-1):
                POS['SubLatt'][i].append([])
        if(i!=POS['EleNum']-1):
            POS['SubLatt'].append([[]])
    return POS
#######################

#######################
def clusterE(ClusterLst,ClusterCoef):
    '''
    Calculate total energy

    Args:
        ClusterLst: List of indentified clusters
        ClusterCoef: ECI of each cluster
    '''
    #ClusterCount = [];
    #for i in range(len(ClusterLst)):
    #    ClusterCount.append(len(ClusterLst[i]));
    ClusterCount = countCluster(ClusterLst);
    #print ClusterCount,ClusterCoef, len(ClusterCount), len(ClusterCoef);
    ECE = 0.0;
    ECE = ECE + ClusterCoef[0];
    #print ECE;
    for i in range(len(ClusterCount)):
        ECE = ECE + ClusterCount[i]*ClusterCoef[i+1];
        #print ECE
    return ECE
#######################

#######################
def dismatswap(dismat,Ind1,Ind2):
    '''
    Update the distance matrix

    Args:
        dismat: distance matrix
        Ind1, Ind2: the indexes of two atoms that swap positions
    '''
    lendismat = len(dismat[1]);
    tmp = dismat[Ind1][:];
    dismat[Ind1][:] = dismat[Ind2][:];
    dismat[Ind2][:] = tmp;
    for i in range(len(dismat[1])):
        if (i!=Ind1)&(i!=Ind2):
            dismat[i][Ind1] = dismat[Ind1][i];
            dismat[i][Ind2] = dismat[Ind2][i];
    tmp = dismat[Ind1][Ind2];
    dismat[Ind1][Ind2] = dismat[Ind1][Ind1];
    dismat[Ind1][Ind1] = tmp;
    tmp = dismat[Ind2][Ind1];
    dismat[Ind2][Ind1] = dismat[Ind2][Ind2];
    dismat[Ind2][Ind2] = tmp;

    return dismat;
#######################

#######################
def dismatinsert(dismat,Ind1,Ind2):
    if (Ind1 > Ind2):
        for i in range(Ind2,Ind1):
            dismatswap(dismat,i,Ind1);
    elif (Ind1 < Ind2):
        #print Ind1,Ind2
        for i in range(Ind2,Ind1,-1):
            dismatswap(dismat,i,Ind1);
    return dismat;
#######################


#######################
def ceFit(Energy,ClusterCount):
    '''
    Least square fitting

    Args:
        TotalE: List of total energies
        ClusterCount: Cluster analaysis results, number of each clusters
        EffInd: Effective cluster index
    '''
    #ClusterCount = countCluster(ClusterLst)
    Energy = np.array(Energy);
    ClusterCount = np.array(ClusterCount);
    OneMat = np.ones([len(ClusterCount),1]);
    A = np.hstack([ClusterCount,OneMat]);
    FitResult = np.linalg.lstsq(A,Energy);
    CECoeff = FitResult[0];
    CECoeff = CECoeff.tolist(); CECoeffLen = len(CECoeff);
    Const = CECoeff[CECoeffLen-1];CECoeff = CECoeff[0:CECoeffLen-1]; 
    CECoeff.insert(0,Const);
    return CECoeff
#######################

#######################
def countCluster(ClusterLst):
    ClusterCount = [];
    for i in range(len(ClusterLst)):
        #print('ClusterLst='+str(ClusterLst[i]))
        ClusterCount.append(len(ClusterLst[i]));
    return ClusterCount
#######################

#######################
def rmseCalc(Energy,Energy_Pre):
    Energy = np.array(Energy); Energy_Pre = np.array(Energy_Pre);
    RMSE = np.sqrt(np.mean((Energy_Pre-Energy)**2));

    return RMSE
#######################

#######################
def cvEvaluate(Energy,ClusterLst):
    ClusterCount = countCluster(ClusterLst);
    CVError = np.array([]);
    for ind, E0 in enumerate(Energy):
        Etmp = list(Energy); Lsttmp = list(ClusterLst)
        Etmp.remove(E0); 
        Lst = Lsttmp[ind]; Lsttmp.remove(Count);
        CECoefftmp = CEFit(Etmp,Lsttmp);
        EPretmp = clusterE(Lsttmp,CECoefftmp);
        np.append(CVError,EPretmp-Etmp);
    CVScore = np.sqrt(np.mean(CVError**2));
    CVError = list(CVError);

    return CVError, CVScore;
#######################

#######################
def ceFind(SubLatt,POSRef,NCut=3,Isprint=0,DCut='default'):
    '''
    Method to find the clusters with a given reference lattice

    Args:
        SubLatt: the projection of solid solution into reference lattice
                 something like [[0,1],[1,2],[3,4]];
        POSRef: POSCAR dictionary for reference lattice
        NCut: Cutoff size of clusters (default: 3)
        DCut: Cutoff length of each dimension of the cluster 
              (default: Half of the box size)
    '''
    print('#############################################################');
    if DCut == 'default':
        DCut = 100.0; TS = 0.3;
        for i in range(3):
            DMax = max(POSRef['Base'][i])/2.0 + TS;
            if DMax < DCut:
                DCut = DMax;
    print('Cutoff cluster length is %f A' %DCut);
    NSubLatt = POSRef['EleNum']; ClusterDesLst = [];
    PrimDistLst = []; AllPrimLattLst = [];
    ClusterNum = []; IndMax = max(SubLatt[-1]);

    FreeSubLatt = []; FreePrim = [];
    for i in range(NSubLatt):
        if len(SubLatt[i]) > 1:
            FreePrim.append(i);
        FreeSubLatt.append(SubLatt[i][0:-1]) #Get rid of last one
    NFreePrim = len(FreePrim); NFreeSubLatt = len(FreeSubLatt);
    FreePrim = np.array(FreePrim); FreeSubLatt = np.array(FreeSubLatt);
    #print(NFreePrim,NFreeSubLatt,FreePrim,FreeSubLatt);

    for N in range(2,NCut+1):
        PrimIndLst = [0]*N;
        PrimDistLst.append([]); AllPrimLattLst.append([]);
        while (PrimIndLst[-1]<=NFreePrim-1):
            #print(PrimIndLst);
            PrimLattLst = list(FreePrim[PrimIndLst]);
            AllPrimLattLst[N-2].append(PrimLattLst);
            DistLst = findCluster(POSRef,PrimLattLst,DCut);
            PrimDistLst[N-2].append(DistLst);
            PrimIndLst=MathKit.lstOrderAdd(PrimIndLst,[NFreePrim-1]*N);
    print('The Distance list of primary lattice is '+str(PrimDistLst));
    print('The cluster made from primary lattice is '+str(AllPrimLattLst));

    ClusterDesLst.append([]); ClusterNum.append(0);
    for SubLatt in FreeSubLatt:
        if SubLatt:
            print(SubLatt);
            for Latt in SubLatt:
                ClusterDesLst[0].append([Latt]);
                ClusterNum[0] += 1;
    #print(ClusterDesLst);

    for N in range(2,NCut+1):
        IndLst = [0]*N;
        ClusterDesLst.append([]); ClusterNum.append(0);
        while (IndLst[-1]<=NFreeSubLatt-1):
            LattLst = list(FreeSubLatt[IndLst]);
            if not [] in LattLst:
                #print('LattLst = '+str(LattLst));
                PrimLattLst = [0]*N;
                for LattInd, Latt in enumerate(LattLst):
                    if Latt in list(FreeSubLatt): 
                        SubInd = list(FreeSubLatt).index(Latt);
                        PrimLattLst[LattInd] = SubInd;
                    else:
                        print('Cannot Latt in FreeSubLatt!!!');
                print('PrimLattLst = '+str(PrimLattLst));
                if PrimLattLst in AllPrimLattLst[N-2]:
                    PrimInd = AllPrimLattLst[N-2].index(PrimLattLst);
                    DistLst = PrimDistLst[N-2][PrimInd];
                else:
                    print('Cannot find the relevant PrimLattLst!!!');
                    break;
                for Dist in DistLst:
                    PermuteLattLst = MathKit.listPermute(LattLst);
                    for PermuteLst in PermuteLattLst:
                        Cluster = [PermuteLst,Dist];
                        if not Cluster in ClusterDesLst[N-2]:
                            ClusterDesLst[N-1].append(Cluster);
                            ClusterNum[N-1] += 1;
            IndLst = MathKit.lstOrderAdd(IndLst,[NFreeSubLatt-1]*N); # Next one
    ClusterSum = sum(ClusterNum);
    print('#############################################################');

    if (Isprint):
        print('#############################################################');
        print('%i indepedent Clusters have been found in this structure' %(ClusterSum));
        for N in range(1,NCut+1):
            print('%i Clusters with %i atoms is given below:' %(ClusterNum[N-1],N));
            ClusterStr = '';
            for Cluster in ClusterDesLst[N-1]:
                ClusterStr+=str(Cluster); ClusterStr+='\t';
            print(ClusterStr);
        print('#############################################################');

    return ClusterSum,ClusterNum,ClusterDesLst;
#######################

#######################
def findCluster(POSRef,LattLst,DCut):
    '''
    Find the Distance Lst for a given PrimAtomLst

    Args:
        POSRef: dictionary of POSRef
        PrimAtomLst: atom list in PrimAtomLst
        DCut: Cutoff distance of

    '''
    NLst = len(LattLst); IndLst = [0]*NLst; GIndLst = [0]*NLst;
    TS = 0.05*NLst;
    DistLst = []; IndLstMax = [];
    for i in range(NLst):
        IndLstMax.append(POSRef['AtomNum'][LattLst[i]]-1);

    while (IndLst[-1]<=IndLstMax[-1]):
        for i, Ind in enumerate(IndLst):
            Indtmp = LattLst[i] - 1;
            GIndLst[i] = Ind + sum(POSRef['AtomNum'][0:Indtmp]);
        Dist = [];
        GrpLst = MathKit.findGrpLst(LattLst);
        for i in range(NLst):
            for j in range(i+1,NLst):
                Distmp = POSRef['dismat'][GIndLst[i]][GIndLst[j]];
                Dist.append(Distmp);
        PermuteGrpLst = MathKit.grporderPermute(LattLst,GrpLst);
        Dist = MathKit.grpSort(Dist,PermuteGrpLst);
        flag = 1;
        for Disttmp in DistLst:
            Distmp = MathKit.grpSort(Disttmp,PermuteGrpLst);
            DistDiff = sum(abs(np.array(Dist)-np.array(Disttmp)));
            if (DistDiff < TS):
                flag = 0;
        if (min(Dist) > 0) & (max(Dist) < DCut) & flag:
            DistLst.append(Dist);
        GrpLst = MathKit.findGrpLst(LattLst);
        IndLst = MathKit.lstGrpAdd(IndLst,IndLstMax,GrpLst);

    return DistLst;
#######################

#######################
def ceAtomFind(SubLatt,POSRef,AtomLatt,AtomInd,NCut=2,Isprint=0,DCut='default'):
    '''
    Method to find the clusters with a given reference lattice, we only
    have partial constraint so only sublattice with one element is eliminated

    Args:
        SubLatt: the projection of solid solution into reference lattice
                 something like [[0,1],[1,2],[3,4]];
        POSRef: POSCAR dictionary for reference lattice
        NCut: Cutoff size of clusters (default: 3)
        DCut: Cutoff length of each dimension of the cluster
              (default: Half of the box size)
    '''
    print('#############################################################');
    if DCut == 'default':
        DCut = 100.0; TS = 0.3;
        for i in range(3):
            DMax = max(POSRef['Base'][i])/2.0 + TS;
            if DMax < DCut:
                DCut = DMax;
        print('Cutoff cluster length is %f A' %DCut);
    NSubLatt = POSRef['EleNum']; KClusterDesLst = [];
    PrimDistLst = []; AllPrimLattLst = [];
    KClusterNum = []; IndMax = max(SubLatt[-1]);
    EffPntLst = [];
    for i in range(POSRef['EleNum']):
        EffPntLst.append([]);

    FreeSubLatt = []; FreePrim = [];
    for i in range(NSubLatt):
        if len(SubLatt[i]) > 1:
            FreePrim.append(i);
            FreeSubLatt.append(SubLatt[i])
        else:
            FreeSubLatt.append([]);
    NFreePrim = len(FreePrim); NFreeSubLatt = len(FreeSubLatt);
    FreePrim = np.array(FreePrim); FreeSubLatt = np.array(FreeSubLatt);

    for N in range(2,NCut+1):
        PrimIndLst = [0]*(N-1);
        PrimDistLst.append([]); AllPrimLattLst.append([]);
        while (PrimIndLst[-1]<=NFreePrim-1):
            PrimLattLst = list(FreePrim[PrimIndLst]);
            AllPrimLattLst[N-2].append(PrimLattLst);
            EffLst, DistLst = \
                    findAtomCluster(POSRef,PrimLattLst,AtomLatt,AtomInd,DCut);
            for LstInd, Lst in enumerate(EffLst):
                for EffInd in Lst:
                    if EffInd not in EffPntLst[LstInd]:
                        EffPntLst[LstInd].append(EffInd);
            PrimDistLst[N-2].append(DistLst);
            PrimIndLst=MathKit.lstOrderAdd(PrimIndLst,[NFreePrim-1]*N);
    print('DistLst found as: %s' %(str(DistLst)));
    print('EffPntLst found as: %s' %(EffPntLst))
    
    POSCluster = copy.deepcopy(POSRef);
    del POSCluster['dismat'];
    POSCluster['EleName'] = []; POSCluster['AtomNum']=[]; 
    POSCluster['AtomSum']=[]; POSCluster['LattPnt'] = [];
    for SubInd, SubLst in enumerate(EffPntLst):
        if SubLst != []:
            POSCluster['EleName'].append(POSRef['EleName'][SubInd]);
            POSCluster['AtomNum'].append(len(SubLst));
            for Ind in SubLst:
                POSCluster['LattPnt'].append(POSRef['LattPnt'][Ind]);
    #print(AtomInd)
    POSCluster['LattPnt'].append(POSRef['LattPnt'][AtomInd]);
    POSCluster['AtomNum'][-1] += 1;
    POSCluster['AtomSum'] = sum(POSCluster['AtomNum']);
    POSCluster['EleNum'] = len(POSCluster['EleName']);
    vp.poswriter('./POSCAR_EffCluster',POSCluster);
    print("Effective cluster writen into POSCAE_EffCluster");

    VacLatt = SubLatt[AtomLatt][-1];
    #print('VacInd = %i' %VacInd);
    for N in range(2,NCut+1):
        IndLst = [0]*(N-1);
        KClusterDesLst.append([]); KClusterNum.append(0);
        while (IndLst[-1] <= NFreeSubLatt - 1):
            LattLst = list(FreeSubLatt[IndLst]);
            if not [] in LattLst: #do not need to count for empty cases
                #print('LattLst = '+str(LattLst));
                PrimLattLst = [0]*(N-1);
                for LattInd, Latt in enumerate(LattLst):
                    if Latt in list(FreeSubLatt):
                        SubInd = list(FreeSubLatt).index(Latt);
                        PrimLattLst[LattInd] = SubInd;
                    else:
                        print('Cannot Latt in FreeSubLatt!!!');
                #print('PrimLattLst = '+str(PrimLattLst));
                if PrimLattLst in AllPrimLattLst[N-2]:
                    PrimInd = AllPrimLattLst[N-2].index(PrimLattLst);
                    DistLst = PrimDistLst[N-2][PrimInd];
                else:
                    print('Cannot find the relevant PrimLattLst!!!');
                    break;
                for Dist in DistLst:
                    PermuteLattLst = MathKit.listPermute(LattLst);
                    for PermuteLst in PermuteLattLst:
                        if VacLatt not in PermuteLst:
                            Cluster = [['Sta']+PermuteLst,Dist];
                            if not Cluster in KClusterDesLst[N-2]:
                                KClusterDesLst[N-2].append(Cluster);
                                KClusterNum[N-2] += 1;
            IndLst = MathKit.lstOrderAdd(IndLst,[NFreeSubLatt-1]*N); # Next one
    KClusterSum = sum(KClusterNum);
    print('#############################################################');

    if (Isprint):
        print('#############################################################');
        print('%i indepedent Clusters have been found in this structure' %(KClusterSum));
        for N in range(2,NCut+1):
            print('%i Clusters with %i atoms is given below:' %(KClusterNum[N-2],N));
            ClusterStr = '';
            for Cluster in KClusterDesLst[N-2]:
                ClusterStr+=str(Cluster); ClusterStr+='\t';
            print(ClusterStr);
        print('#############################################################');

    return EffPntLst,KClusterSum,KClusterNum,KClusterDesLst;
#######################

#######################
def findAtomCluster(POSRef,LattLst,AtomLatt,AtomInd,DCut):
    '''
    Find the Distance Lst for a given PrimAtomLst

    Args:
        POSRef: dictionary of POSRef
        PrimAtomLst: atom list in PrimAtomLst
        DCut: Cutoff distance of

    '''
    SumLst = [sum(POSRef['AtomNum'][0:i]) for i in range(POSRef['EleNum'])];
    print(SumLst);
    NLst = len(LattLst); IndLst = [0]*NLst; GIndLst = [0]*NLst;
    TS = 0.05*NLst;
    DistLst = []; IndLstMax = []; EffLst = [];
    for i in range(POSRef['EleNum']):
        EffLst.append([]);
    for i in range(NLst):
        IndLstMax.append(SumLst[LattLst[i]]-1);
    while (IndLst[-1]<=IndLstMax[-1]):
        for i, Ind in enumerate(IndLst):
            GIndLst[i] = Ind + SumLst[LattLst[i]];
        GIndLst.insert(0,AtomInd); #Insert reference atom
        Dist = [];
        GrpLst = MathKit.findGrpLst(LattLst);
        FullGrpLst = list(GrpLst);
        for GrpIn, Grp in enumerate(FullGrpLst):
            FullGrpLst[GrpIn] = [LstInd + 1 for LstInd in Grp];
        FullGrpLst.insert(0,[0]); #Insert reference atom group
        LattLst.insert(0,AtomLatt); #Insert reference atom latt
        for i in range(NLst+1):
            for j in range(i+1,NLst+1):
                Distmp = POSRef['dismat'][GIndLst[i]][GIndLst[j]];
                Dist.append(Distmp);
        #print(LattLst,FullGrpLst);
        PermuteGrpLst = MathKit.grporderPermute(LattLst,FullGrpLst);
        Dist = MathKit.grpSort(Dist,PermuteGrpLst);
        flag = 1;
        if GIndLst[1] == 50:
            print('Atom50: %s' %str(Dist));
        for Disttmp in DistLst:
            Distmp = MathKit.grpSort(Disttmp,PermuteGrpLst);
            DistDiff = sum(abs(np.array(Dist)-np.array(Disttmp)));
            if (DistDiff < TS):
                flag = 0;
        if (min(Dist) > 0) & (max(Dist) < DCut):
            if flag:
                DistLst.append(Dist);
            for i, GInd in enumerate(GIndLst[1:]):
              #print(GIndLst,LattLst,EffLst,i);
              if GInd not in EffLst[LattLst[i+1]]:
                  #print(EffLst,LattLst,Dist,GIndLst);
                  EffLst[LattLst[i+1]].append(GInd);
                  #print(EffLst);
        IndLst = MathKit.lstGrpAdd(IndLst,IndLstMax,GrpLst);
        LattLst.remove(LattLst[0]); GIndLst.remove(GIndLst[0]);

    return EffLst,DistLst
#######################

