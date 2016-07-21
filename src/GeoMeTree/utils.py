#! /usr/bin/env python

from GeoMeTree.Node import *
from string import *
import math

#================================================================================#
# Functions for Norm ...

def snorm(c): #squared norm
    n=0
    for i in c:
        n+=i**2        
    return n

def norm(c):return math.sqrt(snorm(c))


#================================================================================#
# Functions for Distance Computation
def eucl(c):
    sum=0
    for i in range(0,len(c[0])):
        sum+=(c[0][i]-c[1][i])**2
    return math.sqrt(sum)

def coordinates(r_equ,l_equ,s): #get coordinates for the equations at point s
    p=[]
    for e in l_equ:
        x=0
        if e[0]==-1 or e[0]>=s:x=e[1]*s+e[2]
        if abs(x)<1e-10:x=0
        p.append(x)
    for e in r_equ:
        x=0
        if e[0]!=-1 and e[0]<=s:x=e[1]*s+e[2]
        if abs(x)<1e-10:x=0
        p.append(x)
    return p


def dist_for_geod(r_equ,l_equ,slist=[]):

    def extract_s_from_equ(r_equ,l_equ):
        s=set([0,1])
        for equ in r_equ,l_equ:
            for e in equ:
                if e[0]==-1:continue
                s.add(e[0])
        return sorted(list(s))
    
    if len(slist)==0:slist=extract_s_from_equ(r_equ,l_equ)
    p0=coordinates(r_equ,l_equ,slist[0])
    dist=[]
    for i in range(1,len(slist)):
        p1=coordinates(r_equ,l_equ,slist[i])
        dist.append(eucl([p0,p1]))
        p0=p1
    return dist

def equations(s,bl1,bl2,dec=""): #create right and left equations for all splits

    if not dec: splits=list(bl1.keys())
    else: splits=[i for i in list(bl1.keys()) if i.startswith("%s/"%dec)]
    
    only1=[]; only2=[]; both=[]
    for sp in splits:
        if bl1[sp]==0:only2.append(sp)
        elif bl2[sp]==0:only1.append(sp)
        else:both.append(sp)
    only1.sort(); only2.sort(); both.sort() #the splits lists
    
    ttimes={}.fromkeys(splits)
    for t in list(s.keys()):
        for l in s[t]: #pair of left and right split lists
            for j in l: #all splits in that split list
                ttimes[str(j)]=t
    
    l_equ=[[ttimes[i],-bl1[i]/ttimes[i],bl1[i]] for i in only1]
    
    r_equ=[[ttimes[i],bl2[i]/(1-ttimes[i]),bl2[i]-bl2[i]/(1-ttimes[i])] for i in only2]
    r_equ+=[[0,bl2[i]-bl1[i],bl1[i]] for i in both]
    
    return [l_equ,r_equ],only1+only2+both
    
   

#================================================================================#
# Functions for Indexing ...

def index(plist,dim):

    def binomial(n, k):
        if k==0:return 1
        if 2*k>n:return binomial(n,n-k)
        r=n
        for i in range(2,k+1):
            r*=(n+1-i)
            r/=i
        return r
    
    def binomialsum(n,k): # \sum_{i=0}^{k} (n choose i)
        if k==-1:return 0
        s=1
        p=1
        for i in range(0,k):
            p *= (n-i)/float(i+1)
            s+=p
        return int(round(s))
    
    def binomialsum2(n,k,a): #\sum_{i=1}^{k} (n-i choose a)
        s=0
        for i in range(1,k+1):
            s+=binomial(n-i,a)
        return s

    def index2(plist,dim):
        k=len(plist)
        if k==0: return 0
        min=plist.pop(0)
        s=binomialsum2(dim,min,k-1)
        return s+index2(list(map(lambda x,z=1+min:x-z,plist)),dim-min-1) #reduce plist by min+1


    k=len(plist)
    start=binomialsum(dim,k-1)
    return start+index2(plist,dim)

#================================================================================#
# Functions for Sets ...

def int_to_list(n,dim):
    s=[]
    for i in range(0,dim):
        if n & (1<<i) : s.append(i)
    return s

def int_to_set(n,dim):
    return set(int_to_list(n,dim))

def set_to_int(s):
    n=0
    for i in s:
        n|=1<<i
    return n

#================================================================================#
# Functions for Input ...

def parse_newick(newick_string):
    root=Node(None)
    root.branch=0
    root.name='root'
    ancnode=None
    actnode=root
    bl,name=False,False
    i=1
    for char in newick_string:
        if bl:
            if char==')' or char==',':
                actnode.branch=float(branch)
                bl=False
            else:
                branch+=char
                continue
        elif name:
            if char==')' or char==',' or char=='(' or char==":":
                if len(namestr)==0:
                    namestr="inner_%u"%i
                    i+=1
                actnode.name=namestr
                name=False
            else:
                namestr+=char
                continue
            
        if char=='(':ancnode=actnode
        if char=='(' or char==',':
            actnode=Node(ancnode) #New node
            ancnode.set_child(actnode)
            namestr=""
            name=True
        elif char==')':
            actnode=ancnode
            ancnode=actnode.ancestor
        elif char==':':
            bl=True
            branch=""
    return root

def splits_for_tree(root,term): #Creates split-list,bl-list, boolean poly (True if multifurcations or zero branch lengths), n (number of spp)
    all_leaves=root.get_leaves()
    all_leaves.sort()
    set_leaves=set(all_leaves)
    n=len(all_leaves) #number of taxa
    splits={} #save splits with branch length in dictionary

    for node in root.traverse():
        act_leaves=node.get_leaves()
        act_l=len(act_leaves)
        if act_l==n:continue
        if node.branch==0:continue #do not collect zero-branch-lengths
        if not term and act_l==1:continue
        act_leaves.sort()
        other_leaves=list(set_leaves.difference(set(act_leaves)))
        other_leaves.sort()
        s1=act_leaves[0]
        for i in range (1,act_l):
            s1+="*%s"%act_leaves[i]
        s2=other_leaves[0]
        for i in range (1,n-act_l):
            s2+="*%s"%other_leaves[i]
        if s1<s2:s="%s|%s"%(s1,s2)
        else:s="%s|%s"%(s2,s1)
        if s in splits:splits[s]+=node.branch #root splits one branch
        else:splits[s]=node.branch
    
    split_list=list(splits.keys())
    split_list.sort()
    bl=[splits[s] for s in split_list]
    return split_list,bl,n

#================================================================================#
# Functions for Output ...

def path_to_str(s,L=[],R=[]): #three Lists or s-dictionary
    if type(s)==dict:
        new_s=list(s.keys());new_s.sort()
        L=[];R=[]
        for i in new_s:
            L+=[s[i][0]]
            R+=[s[i][1]]
        s=new_s
    out="t\tL\tR\n"
    for i in range(0,len(s)):
        if not L[i]:L[i]='-'
        if not R[i]:R[i]='-'
        out+="%1.6f\t%s\t%s\n"%(s[i],"*".join([str(j) for j in L[i]]),"*".join([str(j) for j in R[i]]))
    return out

def ind_trans_points(times,r_equ,l_equ):
    out=""
    for t in times:
        out+="g(%1.4f) = ( %s )\n" %(t,", ".join(["%1.4f"%i for i in coordinates(r_equ,l_equ,t)]))

    return out

def trans_points(s,bl1,bl2,dec=""): #three dictionaries, if dec set, then only results for one decomposition
    equ,spl=equations(s,bl1,bl2,dec)
    out="Splits:     ("
    for sp in spl:
        out+="%s%s,"%(" "*(7-len(sp)),sp)
    out=out[:-1]+" )\n"

    l_equ,r_equ=equ
    times=list(s.keys());times.sort()
    return out+ind_trans_points(times,r_equ,l_equ)

def shorter_split(long_split):
    splt=long_split.split("|")
    k,l=splt[0].count("*"),splt[1].count("*")
    if l<k: return splt[1]
    else: return splt[0]  

#================================================================================#
