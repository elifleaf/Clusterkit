#!/usr/bin/env python
##########################################################################
#Nanokit module for Python
#Last editted by Dr. Bin Ouyang 2017.07.03
##########################################################################
#This is the lib containing all the methods for Kinetic Monter Carlo
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
import CELib, MathKit
import copy

class kmcobj(object):
    '''
    A class for kmc simulator
    '''
    #def __init__(self,POS,KClusterdes,Clusterdes,KCECoeff,\
    #        CECoeff,Radius,AtomSite,VacSite,AtomIndLst,VacIndLst,Fre0,BoltT):
    def __init__(self,POS,KClusterdes):
        '''
        Initialization;

        parameter POS: the dictionary containg position information
        parameter Clusterdes: description of clusters
        One example of Clusterdes is like:
        [['Sta',1],[2.2]], it describe the dimer made from hoping
        atom and one atom from type 1 group. The length of the dimer
        is 2.2 Angstrom. There is another option of keyword other than
        'Sta', which is 'Fin'. So [['Fin',1],[2.2]] describe the dimer
        made from site of hopping destination (usually occupied by a 
        vacancy).

        '''
        self.POS = POS;
        self.KClusterdes = KClusterdes; #Kinetic clusters

    def __str__(self):
        '''
        printable version of the class
        '''
        Notes = '\nKinetic cluster expansion object\n';
        Notes += '###################################\n';
        Notes += 'Base Vector of the simulating lattice:\n';
        for i in range(3):
            Notes += str(self.POS['Base'][i]) + '\n';
        Notes += '###################################\n';
        Notes += 'Number of atomic species: ' + str(self.POS['EleNum']) + '\n';
        Notes += 'Number of atoms: ' + str(self.POS['AtomSum']) + '\n';
        Notes += 'Name of atoms: ';
        for Ele in self.POS['EleName']:
            Notes += Ele + ' ';
        Notes += '\n';
        Notes += '###################################\n';
        for Site in self.AtomSite:
            Notes += 'Diffusion species: '+self.POS['EleName'][Site];
        Notes += '\n###################################\n';
        Notes += 'Energetic Clusters Considered: '
        Notes += str(self.Clusterdes) + '\n';
        Notes += '###################################\n';
        Notes += 'Kinetic Clusters Considered: '
        Notes += str(self.KClusterdes) + '\n';
        Notes += '###################################\n';
        Notes += 'Diffusion time reset as %f s\n' %self.diffusetime;
        Notes += '###################################\n';

        return Notes;

    def kmcSetup(self,Clusterdes,KCECoeff,CECoeff,Radius,AtomSite,VacSite,\
            AtomIndLst,VacIndLst,Fre0,BoltT):
        """
        Set up a KMC simulation
        """
        self.Clusterdes = Clusterdes; #Regular clusters
        self.KCECoeff = KCECoeff; self.CECoeff = CECoeff;
        self.Radius = Radius;
        self.AtomSite = AtomSite; self.VacSite = VacSite;
        self.AtomIndLst = AtomIndLst; self.VacIndLst = VacIndLst;
        self.Fre0 = Fre0; self.BoltT = BoltT;
        self.CELst = CELib.clustercount(self.Clusterdes,self.POS); #Energetic clusterlist
        self.diffusetime = 0.0;


    def kclusterEvaluate1(self,AtomInd):
        '''
        Evaluate the number of effective clusters within certain enviroment
        In this version AtomInd is the moving atomind
        '''
        if not self.POS.has_key('dismat'):
            print('Did not find dismat key in POS,\
                    Please create the distance matrix first');
            return None;
        TS = 0.2; #The variation range of cluster length
        ClusterNum = len(self.KClusterdes);
        ClusterLst = [[] for i in range(ClusterNum)];
        for CInd, Cluster in enumerate(self.KClusterdes):
            #print Cluster
            CSize = len(Cluster[0]);
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
                GIndLst.insert(0,AtomInd);
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

    def kclusterEvaluate(self,StaInd,FinInd):
        '''
        Evaluate the number of effective clusters within certain enviroment
        In the current version, StaInd and FinInd is the same thing
        '''
        if not self.POS.has_key('dismat'):
            print('Did not find dismat key in POS,\
                    Please create the distance matrix first');
            return None;
        TS1 = 0.98; TS2 = 1.02; #The range of length variation
        ClusterNum = len(self.KClusterdes);
        #vp.poswriter('./POSCAR_debug',self.POS);
        #print(self.POS['dismat'][49][68]);
        #print(self.POS['LattPnt'][68]);
        ClusterLst = [[] for i in range(ClusterNum)];
        for ClusterInd, Cluster in enumerate(self.KClusterdes):
            ClusterSize = len(Cluster[0]);
            if Cluster[0][0] == 'Sta':
                GInd0 = StaInd;
                for Ind1 in range(self.POS['AtomNum'][Cluster[0][1]]):
                    GInd1 = Ind1 + sum(self.POS['AtomNum'][0:Cluster[0][1]]);
                    Dis_01 = self.POS['dismat'][GInd0][GInd1];
                    #if (GInd0 == 68) & (GInd1 == 49):
                       #print('Dis_01=',Dis_01);
                       #print('Cluster[1][0]=',Cluster[1][0]);
                    if (Dis_01*TS1<Cluster[1][0]) & (Dis_01*TS2>Cluster[1][0]):
                        #check if the cluster is already in the list
                        #if so, print 'yes'
                        noexist = 1;
                        for i in range(len(ClusterLst[ClusterInd])):
                            if (sorted([GInd0,GInd1])==sorted(ClusterLst[ClusterInd][i])):
                                noexist = 0;
                                #print GInd0, GInd1, ClusterLst[ClusterInd][i]
                        if noexist:
                            ClusterLst[ClusterInd].append([GInd0,GInd1]);
            elif Cluster[0][0] == 'Fin':
                GInd0 = FinInd;
                for Ind1 in range(self.POS['AtomNum'][Cluster[0][1]]):
                    GInd1 = Ind1 + sum(self.POS['AtomNum'][0:Cluster[0][1]]);
                    Dis_01 = self.POS['dismat'][GInd0][GInd1];
                    if (Dis_01*TS1<Cluster[1][0]) & (Dis_01*TS2>Cluster[1][0]):
                        #check if the cluster is already in the list
                        #if so, print 'yes'
                        noexist = 1;
                        for i in range(len(ClusterLst[ClusterInd])):
                            if (sorted([GInd0,GInd1])==sorted(ClusterLst[ClusterInd][i])):
                                noexist = 0;
                                #print GInd0, GInd1, ClusterLst[ClusterInd][i]
                        if noexist:
                            ClusterLst[ClusterInd].append([GInd0,GInd1]);
        return ClusterLst;
    
#    def saddTab(self,Radius,SiteLst):
#        '''
#        Find all the saddle locations
#        SiteLst: sublattice list that can swap with vacancy
#        '''
#        SiteLst.append(len(self.POS['EleName'])); #Append the Vac sublatt
#        AtomIndLst = [];
#        for Site in SiteLst:
#            SubLst = range(POS['AtomNum'][Site]);
#            PntLst += [ind+sum(POS['AtomNum'][0:Site]) for ind in SubLst];
#        self.

    def rateCalc(self):
        '''
        Calculate all the rates in specified reactions
        '''
        CELen = len(self.CECoeff); KCELen = len(self.KCECoeff);
        self.RateLst = {}; #List to store the rates of all reactions
        for Vackey in self.Reaction:
            self.RateLst[Vackey] = [];
            for Atom in self.Reaction[Vackey]:
                #print(Atom,Vackey);
                for i in range(self.POS['EleNum']):
                    if Atom < sum(self.POS['AtomNum'][0:i]):
                        AtomLatt = i; break;
                Eb=self.rateCalcEach(Atom,Vackey,AtomLatt);
                #print(type(self.Fre0),type(Eb),type(self.BoltT));
                Omega0 = self.Fre0*math.exp(-Eb/self.BoltT);
                self.RateLst[Vackey].append([Atom,Vackey,Omega0]);
        print('RateLst = '+str(self.RateLst));
        self.totalRateCalc();

    def rateCalcEach(self,AtomInd,VacInd,AtomLatt):
        self.makeSadd(AtomInd,VacInd);
        KClusterLst = self.kclusterEvaluate(AtomInd,AtomInd);
        #print(str(KClusterLst));
        ####Calculate KRA
        KRA = CELib.clusterE(KClusterLst,self.KCECoeff);
        #print('KRA='+str(KRA))
        self.quitSadd(AtomInd,VacInd);
        ####Calculate E_Fin and E_Sta
        #print(len(self.CELst),len(self.CECoeff))
        E_Sta = CELib.clusterE(self.CELst,self.CECoeff);
        #print('ESta='+str(E_Sta));
        ###Swap Positions, note the POS dict is not changed because no need for it
        self.POS['dismat'] = CELib.dismatswap(self.POS['dismat'],AtomInd,VacInd);
        self.CELst = CELib.clusterswap(self.Clusterdes,self.POS,self.CELst,\
                AtomLatt,self.VacSite,AtomInd,VacInd);
        E_Fin = CELib.clusterE(self.CELst,self.CECoeff);
        #print('E_Fin='+str(E_Fin));
        ###Swap Positions back
        self.POS['dismat'] = CELib.dismatswap(self.POS['dismat'],AtomInd,VacInd);
        self.CELst = CELib.clusterswap(self.Clusterdes,self.POS,self.CELst,\
                AtomLatt,self.VacSite,AtomInd,VacInd);
        ####Calculate Eb=KRA+abs(E_Fin-E_Sta)/2;
        Eb = KRA+(E_Fin-E_Sta)/2;
        #print('Eb='+str(Eb));
        #print('KRA = %f, E_Sta = %f, E_Fin = %f, Eb = %f' %(KRA, E_Sta, E_Fin, Eb));
        return Eb;

    def makeSadd(self,StaInd,FinInd):
        self.POS['LattPnt'][StaInd] = list((np.array(self.POS['LattPnt'][StaInd]) + \
                np.array(self.POS['LattPnt'][FinInd]))/2.0);
        Pnt1 = np.array(self.POS['LattPnt'][StaInd]);
        #Update the dismat key
        for i in range(self.POS['AtomSum']):
            Pnt2 = np.array(self.POS['LattPnt'][i]);
            PntDis = Pnt1 - Pnt2;
            #if i == 49:
            #    print(StaInd,i,PntDis);
            #Adjust according to periodicity
            for j in range(3):
                if (PntDis[j]>0.5):
                    PntDis[j] = 1 - PntDis[j];
                elif (PntDis[j]<-0.5):
                    PntDis[j] = PntDis[j] + 1;
                else:
                    PntDis[j] = abs(PntDis[j]);
            PntDis = np.dot(PntDis,self.POS['Base']);
            self.POS['dismat'][StaInd][i] = math.sqrt(PntDis[0]**2 + \
                    PntDis[1]**2 + PntDis[2]**2);
            #if i == 49:
            #    print('distance:'+str(math.sqrt(PntDis[0]**2 + PntDis[1]**2 + PntDis[2]**2)));
            self.POS['dismat'][i][StaInd] = self.POS['dismat'][StaInd][i];

    def quitSadd(self,StaInd,FinInd):
        self.POS['LattPnt'][StaInd] = list(2.0*np.array(self.POS['LattPnt'][StaInd]) -\
                np.array(self.POS['LattPnt'][FinInd]));
        Pnt1 = np.array(self.POS['LattPnt'][StaInd]);
        for i in range(self.POS['AtomSum']):
            Pnt2 = np.array(self.POS['LattPnt'][i]);
            PntDis = Pnt1 - Pnt2;
            for j in range(3):
                if (PntDis[j]>0.5):
                    PntDis[j] = 1 - PntDis[j];
                elif (PntDis[j]<-0.5):
                    PntDis[j] = PntDis[j] + 1;
                else:
                    PntDis[j] = abs(PntDis[j]);
            PntDis = np.dot(PntDis,self.POS['Base']);
            self.POS['dismat'][StaInd][i] = math.sqrt(PntDis[0]**2 + \
                    PntDis[1]**2 + PntDis[2]**2);
            self.POS['dismat'][i][StaInd] = self.POS['dismat'][StaInd][i];

    def rateUpdate(self,AtomInd,VacInd):
        '''
        Update the rates for all reactions
        '''
        CELen = len(self.CECoeff); KCELen = len(self.KCECoeff);
        self.reactionUpdate(AtomInd,VacInd); #Identified all the reactions
        #Now we need to recalculate the reaction rates
        for Vackey in self.Reaction:
            if (Vackey == VacInd): #update the rates related to VacInd
                #Since vacancy is moved, we need to calculate new rates
                self.RateLst[Vackey] = [];
                for Atom in self.Reaction[Vackey]:
                    for i in range(self.POS['EleNum']):
                        if Atom < sum(self.POS['AtomNum'][0:i]):
                            AtomLatt = i; break;
                    Eb = self.rateCalcEach(Atom,Vackey,AtomLatt);
                    Omega0 = self.Fre0*math.exp(-Eb/self.BoltT);
                    self.RateLst[Vackey].append([Atom,Vackey,Omega0]);
            else:
                if AtomInd in self.Reaction:
                    Ind = self.Reaction.index(AtomInd);
                    for i in range(self.POS['EleNum']):
                        if Atom < sum(self.POS['AtomNum'][0:i]):
                            AtomLatt = i; break;
                    Eb = self.rateCalcEach(AtomInd,Vackey,AtomLatt);
                    Omega0 = self.Fre0*math.exp(-Eb/self.BoltT);
                    self.RateLst[Vackey][Ind] = \
                            [AtomInd,Vackey,Omega0];
        self.totalRateCalc();


    def reactionUpdate(self,AtomInd,VacInd):
        '''
        Update the reaction available
        '''
        #Update POS dictionary
        tmp=list(self.POS['LattPnt'][VacInd]);
        self.POS['LattPnt'][VacInd] = list(self.POS['LattPnt'][AtomInd]);
        self.POS['LattPnt'][AtomInd] = list(tmp);
        self.POS['dismat'] = CELib.dismatswap(self.POS['dismat'],AtomInd,VacInd);

        #Update Reaction dictionary
        #***Hint: When swap position, only the lattice point is swapped, not the indexes
        #1. Search the list of reaction for the new Vac position
        self.Reaction[VacInd] = [];
        for Atom in self.AtomIndLst: #Iterate the atom sites
            if (self.POS['dismat'][VacInd][Atom] <= self.Radius):
                self.Reaction[VacInd].append(Atom);
        #2. For other vacancy adjacent to AtomInd, it would be a vacancy at that location
        #   instead. So just remove those reactions
        for VacKey in self.Reaction:
            if AtomInd in self.Reaction[VacKey]:
                if (self.POS['dismat'][VacKey][AtomInd] > self.Radius):
                    self.Reaction[VacKey].remove(AtomInd);
            elif (self.POS['dismat'][VacKey][AtomInd] <= self.Radius):
                self.Reaction[VacKey].append(AtomInd);

        return AtomInd, self.POS['dismat'][AtomInd][VacInd];

    def reactionFind(self,IsP=False):
        '''
        Identify the possible reaction pathway
        This version would consider not only 1NN
        '''
        #TS0 = 0.2 #0.2 Angstrom of error tolerance
        self.Reaction = {};
        for Ind in self.VacIndLst:
            self.Reaction[Ind] = [];
            for AtomInd in self.AtomIndLst: #Iterate the atom sites
                if (self.POS['dismat'][Ind][AtomInd] <= self.Radius):
                    self.Reaction[Ind].append(AtomInd);
        if IsP:
            print('#############################');
            print('The identified reactions is given below:');
            print(str(self.Reaction)); #Print the find reactions
            print('#############################');

    def totalRateCalc(self):
        '''
        Caclulate the time length of a certain event:
        '''
        self.totalRates = 0.0;
        for Vackey in self.Reaction:
            Rates = np.array(self.RateLst[Vackey]);
            self.totalRates += sum(Rates[:,2]);

    def findEvent(self,rdnum):
        '''
        Find the jumping event according to random number
        '''
        Ratetmp = 0.0
        for Vackey in self.Reaction:
            for Ind in range(len(self.RateLst[Vackey])):
                flag1 = (rdnum>=Ratetmp);
                flag2 = (rdnum < Ratetmp+self.RateLst[Vackey][Ind][2]/\
                        self.totalRates);
                if flag1 & flag2 :
                    return self.RateLst[Vackey][Ind][0],self.RateLst[Vackey][Ind][1];
                Ratetmp += self.RateLst[Vackey][Ind][2]/self.totalRates;
        print "something must go wrong!!!!";
        print "Ratetmp = %f, rdnum = %f" %(Ratetmp,rdnum);


