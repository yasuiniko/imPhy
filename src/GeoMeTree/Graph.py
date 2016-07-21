#! /usr/bin/env python

from GeoMeTree.utils import *
from copy import *

#================================================================================#
# Classes for the Graph 

class Path: #goes from 1 to 0!
    def __init__(self,edgelist):
        self.edges=edgelist
        self.__inverse=False

    def equations(self,bl1,bl2):
        r_equ=[]
        l_equ=[]
        for e in self.edges:
            act_l,act_r=e.get_equations(bl1,bl2)
            l_equ+=act_l
            r_equ+=act_r
        return (l_equ,r_equ)

    def distance_old(self,bl1,bl2,final=1): #compute distance for path and given branch lengths

        l_equ,r_equ=self.equations(bl1,bl2)
       
        slist=[e.s for e in self.edges]
        slist.reverse() #reverse should not change distance computation
        if slist[-1]!=final:slist.append(final)
        dist=[]
        p0=coordinates(r_equ,l_equ,0)
        for s in slist:
            p1=coordinates(r_equ,l_equ,s)
            dist.append(eucl([p0,p1]))
            p0=p1
            
        return dist

    def distance(self,bl1,bl2): #dont need branch lengths any more!
        dlist=[e.dist for e in self.edges]
        return [norm(dlist)]

    def get_L(self):
        if not self.__inverse:
            llist=[int_to_list(e.L,Orthant.dim1) for e in self.edges]+[[]] 
            llist.reverse()
        else:llist=[int_to_list(e.R,Orthant.dim2) for e in self.edges]+[[]]
        return llist
    
    def get_R(self):
        if not self.__inverse:
            rlist=[int_to_list(e.R,Orthant.dim2) for e in self.edges]+[[]]
            rlist.reverse()
        else:rlist=[int_to_list(e.L,Orthant.dim1) for e in self.edges]+[[]]
        return rlist
    
    def get_s(self):
        if not self.__inverse:
            slist=[e.s for e in self.edges]+[0]
            slist.reverse()
        else:slist=[1-e.s for e in self.edges]+[1]
        return slist

    def __repr__(self):return str(self.edges)

    def inverse(self):
        self.__inverse= not self.__inverse

#================================================================================#

class IEdge: #Representation of an edge
    def __init__(self,Neg,todoPos,L,R,anc=None,s=-1):
        self.anc=anc #ancestor orthant
        self.s=s     #transition time
        self.dist=0  #distance of edge in subtreespace
        self.L=L     #list of Left indices
        self.R=R     #list of Right indices

        Neg=self.__verify(Neg)            #create new Neg-set
        todoPos=self.__adapt(Neg,todoPos) #create new todoPos-set
        self.__compute_id(todoPos)        #compute id of orthant at tip of edge, and save it in the list for future processing, stores id of orthant in self.__id

    def __verify(self,Neg): #adapt L with Coordinates necc to shrink

        def is_right_list_compatible(s,l,adj): # s has to be compatible with every split in l
            for split in l:
                if not adj[s][split]:return False
            return True

        Rs=int_to_list(self.R,Orthant.dim2)
        for n in int_to_list(Neg,Orthant.dim1):
            if not is_right_list_compatible(n,Rs,Orthant.adj):
                self.L+=1<<n
                
        return self.anc.Neg-self.L
        
    def __adapt(self,Neg,todo): #coordinates which can be expanded

        def is_left_list_compatible(s,l,adj): # s has to be compatible with every split in l
            for split in l:
                if not adj[split][s]:return False
            return True
        
        Ns=int_to_list(Neg,Orthant.dim1)
        for p in int_to_list(todo,Orthant.dim2): 
            if is_left_list_compatible(p,Ns,Orthant.adj):
                self.R+=1<<p
        return self.anc.todo-self.R
        
    def __compute_id(self,todo): #compute unique id from todo-pattern
        max=2**(Orthant.dim2+1) - 1
        simple_id=max^todo #xor, is the difference
        if simple_id in Orthant.IDList:self.__id=Orthant.IDList[simple_id]
        else:
            plist=int_to_list(simple_id,Orthant.dim2)
            self.__id=index(plist,Orthant.dim2)
            Orthant.IDList[simple_id]=self.__id

            
    def compute_s(self,bl1,bl2,orth): #Computation of transition times t, need orthant for optimization
        global opts
        if self.s!=-1:return True  #first node has no ancestor, s already 0
        
        nl=norm([bl1[i] for i in int_to_list(self.L,Orthant.dim1)])
        nr=norm([bl2[i] for i in int_to_list(self.R,Orthant.dim2)])

        self.dist=nl+nr
        self.s=nl/self.dist
        
        if self.anc.edges and self.s <= min([e.s for e in self.anc.edges]): return False #no valid joining point
        
        if not Orthant.opt:return True
        
        if orth.todo==0:return True
        nl=norm([bl1[i] for i in int_to_list(orth.Neg,Orthant.dim1)])
        nr=norm([bl2[i] for i in int_to_list(orth.todo,Orthant.dim2)])
        if self.s > nl/(nl+nr):return False
        
        return True


    def get_id(self):return self.__id

    def create_orthant(self):
        Neg=self.anc.Neg-self.L
        todoPos=self.anc.todo-self.R
        return Orthant(Neg,todoPos,self)

        
    def __eq__(self,other): return (self.__id==other.__id) and (self.anc==other.anc) #same edge, if same id and same ancestor
    
    def __repr__(self):return  "(L: "+str(int_to_set(self.L,Orthant.dim1))+" R: "+str(int_to_set(self.R,Orthant.dim2))+" t:  "+str(self.s)+")" #string representation
    
    def __create_equations(self,bl1,bl2): #create equations from transition times, L and R and branch lengths
        l_equ=[]
        for i in int_to_list(self.L,Orthant.dim1):
            l_equ.append([self.s,-bl1[i]/self.s,bl1[i]])
        r_equ=[]
        for i in int_to_list(self.R,Orthant.dim2):
            r_equ.append([self.s,bl2[i]/(1-self.s),bl2[i]-bl2[i]/(1-self.s)])
        return (l_equ,r_equ)

    def get_equations(self,bl1,bl2): #return equations
        return self.__create_equations(bl1,bl2)

#================================================================================#

#Orthant class variables: dim1, dim2, IDList, adj
class Orthant(object):
    def __init__(self,Neg,todoPos,edge=None):
       if edge: 
           self.edges=[edge]        #list of ingoing edges
           self.__id=edge.get_id()  #id already computed
       else:
           self.edges=[]            #first orthant has no ingoing edge
           self.__id=0

       if type(Neg)==set:Neg=set_to_int(Neg)
       self.Neg=Neg                 #set of indices from T1
       if type(todoPos)==set:todoPos=set_to_int(todoPos)
       self.todo=todoPos            #set of indices from T2 not yet processed

    def __repr__(self):return "Splits from T1: "+str(int_to_set(self.Neg,Orthant.dim1))+" Splits from T2 Todo: "+str(int_to_set(self.todo,Orthant.dim2))+" ID: "+str(self.__id) #string representation    

    def __hash__(self):return self.__id #order of orthants

    def __cmp__(self,other):return self.__id-other.__id #an orthant is smaller if it has smaller id

    def get_id(self):return self.__id

    def get_todo(self):return int_to_set(self.todo,Orthant.dim2)
    
    def clone(self,c): #create new edge starting in orthant where c is changed to pos
        cint=set_to_int(c)
        return IEdge(Neg=self.Neg,todoPos=self.todo-cint,anc=self,L=0,R=cint)
        
#================================================================================#
