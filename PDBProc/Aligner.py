'''
Created on Oct 10, 2021

Super simple aligner using scipy
@author: Chris Forman
'''
import numpy as np
import copy as cp
from scipy.spatial.transform import Rotation as rot

class Aligner(object):
    '''
        A minimal stapler for stapling together XYZ lists of numpty vectors
        based on indices.  Can use this to use Stapler more easily within
        other programs without all the fluff about file handling etc. 
    '''
    def __init__(self):
        '''
        Constructor
        '''
        pass
    
    def convertNumpy3VecToCoordList(self, vecList):
        # converts a list of numpy vectors to a flat coord list
        outlist = []
        for vec in vecList:
            outlist.append(vec[0])
            outlist.append(vec[1])
            outlist.append(vec[2])
        return np.array(outlist)  

    def convertCoordListToNumpy3Vec(self, coordList):
        # converts a flat coord list to a list of 3 vectors
        return np.reshape(np.array(coordList), (int(float(len(coordList)) / 3.0), 3) )  

    def prepVecs(self, inpVecs, indices):
        # slice the numpy array to extract the subset of vectors
        _SubSetVecs = inpVecs[indices]  

        # compute the Centre of mass of the *sub set*.
        COM = np.sum(_SubSetVecs, 0)/len(_SubSetVecs)

        # move the subset to its centre of gravity
        _VecsCOM = [ vec - COM for vec in _SubSetVecs]
        
        # return the centre of mass and the shifted sub set 
        return COM, _VecsCOM
        
    # takes a list of static and mobile XYZ points and aligns the 
    # the mobile ones with the static ones.  The alignment rotation
    # matrix is calculated by minimising the alignment between
    # a subset of the static and mobile sets.
    # The rotation matrix is then applied to the full set of 
    # mobile XYZ coords and the new XYZ coords are returned.

    def align(self, staticXYZ, mobileXYZ, staticJoinIndex, mobileJoinIndex, alignPoint='COM'):

        # copy and numpify the input vectors
        _statVecs = np.copy(cp.copy(staticXYZ))
        _mobVecs = np.copy(cp.copy(mobileXYZ))
        
        # get the relevant subsets of the static and mobile xyz set both translated to their own COM 
        staticCOM, _statVecsCOM = self.prepVecs( _statVecs, staticJoinIndex )
        mobileCOM, _mobVecsCOM = self.prepVecs( _mobVecs, mobileJoinIndex )
        
        # find the rot matrix that best aligns the points
        rotBest, rmsd = rot.align_vectors(_statVecsCOM, _mobVecsCOM)
        
        # move the entire mobile block to the SUBSET centre of Mass. 
        _mobVecsCOM_ALL = [ mv - mobileCOM for mv in _mobVecs ]
        
        # apply the rotation matrix
        _newMobVecsCom_ALL = rotBest.apply(np.copy(cp.copy(_mobVecsCOM_ALL)))

        # we can choose how to set the final position
        if alignPoint=='COM':
            # translate the set to the STATIC subset's centre of mass.
            outVecs = [ mv + staticCOM for mv in _newMobVecsCom_ALL ]
        elif alignPoint in mobileJoinIndex:
            # The alignPoint refers to the index in the joining sets of the points that should be superimposed.
            StaticPoint = staticXYZ[staticJoinIndex[alignPoint]]
            MobilePoint = _newMobVecsCom_ALL[mobileJoinIndex[alignPoint]]
            # translate the rotated mobile set so one of it's members perfectly coincides with 
            # one of the static members, retaining the rotated character of the mobile set. 
            outVecs = [ mv - MobilePoint + StaticPoint for mv in _newMobVecsCom_ALL ]
        else:
            # translate the set to the STATIC subset's centre of mass.
            outVecs = [ mv + staticCOM for mv in _newMobVecsCom_ALL ]
            
        return outVecs 


    def isEqual(self, p1, p2, epsilon):
        # checks to see if two vectors are equal. 
        # Zero vector defined as having a magnitude less than epsilon. 
        equal = True # assume vectors are equal
        n = p1 - p2
        nMag =  np.linalg.norm(n)
        
        if ( abs(nMag - 0.0) > epsilon):
            equal = False
            
        return equal

    def rotPAboutAxis(self, p, n, angle):
        # angle in radians
        # rotates the vector P about an axis n by an angle.
    
        nHat = n/np.linalg.norm(n) # ensure we are using a normalised vector
        try:
            outVec = p*np.cos(angle) + np.cross(nHat, p)*np.sin(angle) + nHat*np.dot(nHat, p)*(1- np.cos(angle))
        except ValueError:
            print("valueError")
            outVec = p
    
        return outVec

    def rotPAboutAxisBetweenPoints(self, p, P0, P1, angle):
        # angle in radians
        # rotates the vector p about the axis defined by the point P0 and P1
        if not self.isEqual(p, P0, 1e-10) and not self.isEqual(P0, P1, 1e-10) and not self.isEqual(p, P1, 1e-10):
            r = p - P0
            n = P1 - P0
            nHat = n/np.linalg.norm(n) # ensure we are using a normalised vector
            try:
                r_Rot = r*np.cos(angle) + np.cross(nHat, r)*np.sin(angle) + nHat*np.dot(nHat, r)*(1- np.cos(angle))
                outVec = r_Rot + P0
            except ValueError:
                print("valueError")
                outVec = p
        else:
            outVec = p
        return outVec

    def rotPAboutAxisAtPoint(self, p, p0, n, angle):
        # angle in radians
        # rotates the vector P about an axis n through the point p0 by an angle.
        if not self.isEqual(p, p0, 1e-10):
            r = p - p0
            nHat = n/np.linalg.norm(n) # ensure we are using a normalised vector
            try:
                r_Rot = r*np.cos(angle) + np.cross(nHat, r)*np.sin(angle) + nHat*np.dot(nHat, r)*(1- np.cos(angle))
                outVec = r_Rot + p0
            except ValueError:
                print("valueError")
                outVec = p
        else:
            outVec = p
        return outVec

        
if __name__=="__main__":
    vecs1 = np.random.rand(10,3)
    vecs2 = vecs1 
    vecs2[0] = vecs2[0] + 0.1 * np.random.rand(1,3)
    
    testVecs1 = np.array([[0.0, 0.0, 0.0],
                         [0.0, 0.0, 1.0],
                         [0.0, 1.0, 0.0]])
    
    testVecs2 = np.array([[1.0, 0.0, 2.0],
                         [1.0, -1.0, 2.0],
                         [1.0, 0.0, 3.0]])
    
    AL=Aligner()
    newVecs1 = AL.align(vecs1, vecs2, [0, 1, 2], [0, 1, 2])
    
    import Utilities.fileIO as fIO
    
    fIO.saveXYZ(vecs1, 'c', 'vecs1.xyz')
    fIO.saveXYZ(vecs2, 'c', 'vecs2.xyz')
    fIO.saveXYZ(newVecs1, 'c', 'newVecs1.xyz')
    