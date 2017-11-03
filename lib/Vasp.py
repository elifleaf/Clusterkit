#!/usr/bin/env python
#######################################################################
#Libarary containing all the methods for VASP file processing
#######################################################################
#Format of definition
#functions are all non-capitalized
#dictionaries are all capitalized
#Other variables are first letter capitalized
#######################################################################

__author__ = 'Bin Ouyang'
import os
import math
import numpy as np
import sys
import re

def posreader(PosName='POSCAR'):
    """
    Read the atomic configuration from POSCAR

    Args:
        PosName (str): the name of the POSCAR File, (default: 'POSCAR')
    """
    POS = {} #Initialize the dictionary for POSCAR information
    Fid = open(PosName,'r');
    Line = Fid.readline(); POS['CellName'] = Line.split('\n')[0]; #Comment line
    Line = Fid.readline(); Sline = Line.split(); POS['LattConst'] = float(Sline[0]); #Lattice constant
    POS['Base'] = [[0.0]*3 for i in range(3)] #Initilize the base list
    for i in range(3):
        Line = Fid.readline(); Sline = Line.split(); 
        POS['Base'][i] = [float(Sline[i]) for i in range(3)];
    Line = Fid.readline(); Sline = Line.split(); POS['EleName']=Sline #The name of each element
    POS['EleNum']= len(POS['EleName']);#number of elements involved
    Line = Fid.readline(); Sline = Line.split(); POS['AtomNum'] = [0]*POS['EleNum'];
    POS['AtomSum'] = 0;
    for ind, Num in enumerate(Sline):
       POS['AtomNum'][ind] = int(Num);
       POS['AtomSum'] += int(Num);
    Line = Fid.readline(); Sline = Line.split(); FL=Sline[0][0]; #Check the first letter
    if (FL=='S'):
        POS['IsSel'] = 1; POS['SelMat'] = [['X']*3 for i in range(POS['AtomSum'])]
        Line = Fid.readline(); Sline = Line.split(); FL=Sline[0][0]; #Check the first letter for coord
    else:
        POS['IsSel'] = 0
    # Set up the lattice type
    if (FL=='D') | (FL=='d'):
        POS['LatType'] = 'Direct'
    elif (FL=='C') | (FL=='c'):
        POS['LatType'] = 'Cartesian'
    else:
        print ("Please check the POSCAR file, the lattice type is not direct or cartesian");
    POS['LattPnt'] = [[0.0]*3 for i in range(POS['AtomSum'])] #Initialize lattice points
    for i in range(POS['AtomSum']):
        Line = Fid.readline(); Sline = Line.split();
        POS['LattPnt'][i] = [float(Sline[i]) for i in range(3)];
        if(POS['IsSel']):
            POS['SelMat'][i] = [Sline[i+3] for i in range(3)];
    Fid.close();
    #The current version does not support reading the POSCAR with velocity information!!!!!!!!!!!!!!!!
    return POS

def poswriter(PosName = 'POSCAR',POS):
    """
    Write out the POS into a POSCAR file

    Args:
        PosName: the name of the POSCAR file
        POS: the POS dictionary
    """
    Fid = open(PosName,'w');
    Fid.write('%s ' %POS['CellName']);
    Fid.write('\n')
    Fid.write('%f \n' %POS['LattConst']);
    for i in range(3):
        Fid.write('%f %f %f \n' %(POS['Base'][i][0],POS['Base'][i][1],POS['Base'][i][2]))
    for i in range(POS['EleNum']):
        Fid.write('%s ' %POS['EleName'][i]);
    Fid.write('\n');
    for i in range(POS['EleNum']):
        Fid.write('%i ' %POS['AtomNum'][i]);
    Fid.write('\n');
    if (POS['IsSel']):
        Fid.write('Selective Dynamics \n');
    Fid.write('%s \n' %POS['LatType']);
    for i in range(POS['AtomSum']):
        Fid.write('%f %f %f ' %(POS['LattPnt'][i][0],POS['LattPnt'][i][1],POS['LattPnt'][i][2]))
        if (POS['IsSel']):
            Fid.write('%s %s %s ' %(POS['SelMat'][i][0],POS['SelMat'][i][1],POS['SelMat'][i][2]))
        Fid.write('\n')
    Fid.close();

def isduplicate(LattPnt,RefLatt,TS = 0.001):
    """
    check if LattPnt is in RefLatt

    Args:
        LattPnt: the point to be checked, list
        RefLatt: the reference lattice point, list of list
        TS: the allowered errors of atomic position (default: 0.01)
    """
    IsDup = 0
    for RefInd, Ref in enumerate(RefLatt):
        Ref=np.array(Ref); LattPnt=np.array(LattPnt);
        PntDis = Ref-LattPnt;
        flag0 = abs(PntDis[0])<=0.01;
        flag1 = abs(PntDis[1])<=0.01;
        flag2 = abs(PntDis[2])<=0.01;
        if (flag0 & flag1 & flag2):
            IsDup = 1;
            break;
    return IsDup;

def nebreader(NEBFile='neb.dat'):
    """
    Read the ned.dat generated from nebbarrier.pl script

    Args:
        NEBFile: the neb.dat file
    """
    E_Path = [];
    Fid = open(NEBFile,'r');
    while 1:
        Line = Fid.readline(); 
        if not Line:
            break;
        Sline = Line.split();
        E_Path.append(float(Sline[2]));
    Fid.close();
    E_Sta = E_Path[0]; E_Fin = E_Path[-1]; E_NEB = max(E_Path[1:-1]);

    return E_Sta,E_Fin,E_NEB;

def dismatcreate(POS):
    """
    Create the distance matrix for a given POS

    Args:
        POS: the POS dictionary
    """
    POS['dismat'] = [[0.0]*POS['AtomSum'] for i in range(POS['AtomSum'])];
    for AtomInd1, Pnt1 in enumerate(POS['LattPnt']):
        for AtomInd2, Pnt2 in enumerate(POS['LattPnt']):
            Pnt1=np.array(Pnt1); Pnt2=np.array(Pnt2);
            PntDis = Pnt1-Pnt2;
            for i in range(3):
                if (PntDis[i]>0.5):
                    PntDis[i] = 1-PntDis[i];
                elif (PntDis[i]<-0.5):
                    PntDis[i] = PntDis[i]+1;
                elif (PntDis[i]>=-0.5) & (PntDis[i]<=0.5):
                    PntDis[i] = abs(PntDis[i]);
                else:
                    print "Something is wrong when calculating dist matrix";
            PntDis=np.dot(PntDis,POS['Base']);
            POS['dismat'][AtomInd1][AtomInd2] = math.sqrt(PntDis[0]**2 + \
                    PntDis[1]**2 + PntDis[2]**2);
    return POS;

def findVac(POS,POSRef,Site,RefSite):
    """
    Find the vacancy lattice site, the Vac sublattice will be added in the end

    Args:
        POS: the POS dictionary
        POSRef: the reference POS dictionary without a vacancy
        Site: the sub lattice index in POS
        RefSite: the sub lattice index in RefSite
    """
    #The Vac sublattice should be added in to end
    if POS['AtomSum'] == POSRef['AtomSum']:
        print "There would be no vacant site to add!";
        return;
    if POS['EleName'][-1] == 'Vac':
        print('The vacancy sites have already been specified');
        return;
    POSRef['LattPnt'] = np.array(POS['LattPnt']);
    RefRange = list(np.array(range(POSRef['AtomNum'][RefSite])) + \
            sum(POSRef['AtomNum'][0:RefSite]));
    RefPnts = list(POSRef['LattPnt'][RefRange]);
    CurrrentRange = list(np.array(range(POS['AtomNum'][Site])) + \
            sum(POS['AtomNum'][0:Site]));
    Pnts = list(POS['LattPnt'][CurrentRange]);
    VacPnts = []; 
    POS['EleName'].append('Vac');
    POS['AtomNum'].append(0); POS['EleNum'] += 1;

    for Pnt in RefPnts:
        if not isduplicate(Pnt,Pnts):
            POS['LattPnt'].append(Pnt);
            POS['AtomNum'][-1]+=1;
            POS['AtomSum'] += 1;
    return;

def findMobile(POS1,POS2):
    '''
    Find the mobile atom index after the diffusion step

    Args:
        POS1 POS2: the POS dictionary before and after diffusion step
    '''
    MobileInd = [];
    for Ind1,Pnt1 in enumerate(POS1['LattPnt']):
        if Pnt1 not in POS2['LattPnt']:
            MobileInd.append(Ind1)
    #print len(MobileInd);
    if len(MobileInd) == 1:
        StaInd = MobileInd; FinInd = MobileInd;
        #print StaInd,FinInd;
        return StaInd[0],FinInd[0];
    else:
        print 'Only one moving atom is accepted, go check pls!!!';

