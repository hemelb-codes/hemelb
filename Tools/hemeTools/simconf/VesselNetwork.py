import numpy as np
import vtk
import quantities as pq
from ..utils.IterPairs import IterPairs

def IndexOfLastShared(a, b):
    """Given two sequences that start the same, find the index of the
    last shared value (by bisection search).
    """
    n = min(len(a), len(b))
    
    lo = 0
    hi = n - 1
    
    assert a[lo] == b[lo]
    assert a[hi] != b[hi]
    
    while hi > lo + 1:
        mid = (lo + hi) / 2
        if a[mid] == b[mid]:
            lo = mid
        else:
            hi = mid
    return lo


class Unknown(object):
    """Tag class used in constructing the equations.
    """
    pass
class UnknownPressure(Unknown):
    units = pq.pascal
    pass
class UnknownFlowRate(Unknown):
    units = pq.metre**3 / pq.second
    pass

class VesselSegment(object):
    def __init__(self, ptIdList):
        if len(ptIdList) == 1:
            self.children = []
            self.ptIds = ptIdList[0]
            return

        # Find the index of the last common value across all lists
        lastSharedIdx = min(IndexOfLastShared(a,b) for a,b in IterPairs(ptIdList))
        # Set this segment to the common path
        self.ptIds = ptIdList[0][:lastSharedIdx+1]
        # Downstream points are the rest
        downstreamPtIdList = [x[lastSharedIdx+1:] for x in ptIdList]

        # Group the child segments such that the first value is the
        # same in each group
        childGroups = []
        for dsList in downstreamPtIdList:
            matchGrp = False
            for grp in childGroups:
                # Check to see if dsList starts the same as any existing group
                first = grp[0]
                if dsList[0] == grp[0][0]:
                    # It does so add to that group and set a flag
                    grp.append(dsList)
                    matchGrp = True
                    break
                continue
            
            if not matchGrp:
                # It didn't match any existing group
                childGroups.append([dsList])
                pass
            continue

        self.children = [VesselSegment(grp) for grp in childGroups]
        
        return

    def BuildTerms(self, inletPressure=None, flowRate=None, finalOutletPressure=0.):
        unknowns = []
        
        if inletPressure is None:
            self.InletPressure = UnknownPressure()
            unknowns.append(self.InletPressure)
        else:
            self.InletPressure = inletPressure
            pass

        if flowRate is None:
            self.FlowRate = UnknownFlowRate()
            unknowns.append(self.FlowRate)
        else:
            self.FlowRate = flowRate
            pass
        
        if len(self.children):
            self.OutletPressure = UnknownPressure()
            unknowns.append(self.OutletPressure)
            
            for ch in self.children:
                unknowns.extend(ch.BuildTerms(inletPressure=self.OutletPressure,finalOutletPressure=finalOutletPressure))
                continue
        else:
            self.OutletPressure = finalOutletPressure
            pass
        return unknowns
    
    def BuildSystem(self, prof):
        ans = []
        # Always do flow
        resist = prof.ResistanceForIds(self.ptIds)
        # self.InletPressure - self.OutletPressure == Constant(resist) * self.FlowRate
        #    always unkn          const/unkn                               const/unkn
        eq = {'rhs': 0.0 * pq.pascal}
        eq[self.InletPressure] = 1.
        
        if isinstance(self.OutletPressure, Unknown):
            eq[self.OutletPressure] = -1.
        else:
            eq['rhs'] += self.OutletPressure
            pass

        if isinstance(self.FlowRate, Unknown):
            eq[self.FlowRate] = - resist
        else:
            eq['rhs'] += resist * self.FlowRate
            pass
        ans.append(eq)
        
        if len(self.children):
            # Do continuity with children
            # sum(ch.FlowRate for ch in self.children) == self.FlowRate
            #     always unknown                        known if root
            eq = {'rhs': 0.0 * pq.metre**3 / pq.second}
            
            for ch in self.children:
                eq[ch.FlowRate] = 1.
                
            if isinstance(self.FlowRate, Unknown):
                eq[self.FlowRate] = -1.
            else:
                eq['rhs'] = self.FlowRate
                pass
            ans.append(eq)
            
            for ch in self.children:
                ans.extend(ch.BuildSystem(prof))
            pass
        return ans
    
    def MyPeakVelocities(self, prof):
        if isinstance(self.FlowRate, Unknown):
            Q = prof.SolutionForUnknown(self.FlowRate)
        else:
            Q = self.FlowRate
            
        radii = prof.CentreLineRadii
        return 2 * Q / (np.pi * radii[self.ptIds]**2)
    
    def MaxVelocity(self, prof):
        myMax = self.MyPeakVelocities(prof).max()
        if len(self.children):
            chMax = max(ch.MaxVelocity(prof) for ch in self.children)
            return max(myMax, chMax)
        return myMax

    def MinRadius(self, prof):
        r = prof.CentreLineRadii
        myMin = r[self.ptIds].min()
        if len(self.children):
            chMin = min(ch.MinRadius(prof) for ch in self.children)
            return min(myMin, chMin)
        return myMin

    def SamplePtIds(self, spacing, start=0):
        listOfArrays = [self.ptIds[start::spacing]]
        childStart = spacing - (len(self.ptIds)-start) % spacing
        for ch in self.children:
            listOfArrays.append(ch.SamplePtIds(spacing, childStart))
            continue
        return np.concatenate(listOfArrays)
    
    def MaxLength(self, prof):
        x = prof.CentreLinePoints[self.ptIds]
        dx = x[1:] - x[:-1]
        ds = np.sum(dx**2,axis=-1)**0.5
        
        myLen = ds.sum()
        if len(self.children):
            myLen += max(ch.MaxLength(prof) for ch in self.children)
            pass
        
        return myLen
    pass
