#! /usr/bin/env python

#2008-10-03 exception introduced when too many internal splits, but still error when decomposition is applied

# own libraries
from GeoMeTree.utils import *
from GeoMeTree.Graph import *

# external libraries
from string import *
from copy import *
import sys
import math
import time

#### OPTS ####
approx = False
decomp = True
graph = True
normalize = False
opt = True
outfile = None
term = True
#### OPTS ####

def parse_options():

    from optparse import OptionParser, OptionGroup, SUPPRESS_HELP
    parser = OptionParser()
    group=OptionGroup(parser,"Input options")
    group.add_option("-f", "--file", dest="infile", help="Name of input file, all pairs of trees in the file are evaluated (no default!)",default=None)
    parser.add_option_group(group)

    group=OptionGroup(parser,"Tree lengths options")
    group.add_option("-n", "--norm", dest="normalize",help="Normalize the branch length vectors to norm 1 (default: no normalization)",action="store_true",default=False)
    group.add_option("-t", "--term",dest="term",help="Ignore branch lengths of terminal splits (then also ignored in normalization, default: terminal branch lengths considered)",action="store_false",default=True)
    parser.add_option_group(group)

    group=OptionGroup(parser,"Algorithmic options")
    group.add_option("-o","--opt",dest="opt",action="store_false",default=True,help=SUPPRESS_HELP) #help="Turn off optimization (default: with optimization)"
    group.add_option("-g","--graph",dest="graph",help=SUPPRESS_HELP,action="store_false",default=True) #help="Turn off evaluation of complete graph (default: on)"
    group.add_option("-a","--approx",dest="approx",help="Compute only the approximations, not geodesic path",action="store_true",default=False)
    group.add_option("-d","--decomp",dest="decomp",help="Do not decompose the trees (default: decomposition when possible)",action="store_false",default=True)
    parser.add_option_group(group)

    group=OptionGroup(parser,"Output options") 
    group.add_option("-m","--name",dest="outfile",default="None",help=SUPPRESS_HELP)#,help="Name of tabular output file, if 0: {file}.dist, default None"
    group.add_option("-s","--silent",dest="silent",help=SUPPRESS_HELP,action="store_true",default=False) #no pair file
    group.add_option("-v","--header",dest="header",help="Header name of output file(s) (default: 'pair', then files pair_i_j are generated for each pair)",default="pair")
    parser.add_option_group(group)

    (options, args) = parser.parse_args()
    
    if not options.infile:
        print("\nName of infile required (option -f)\n")
        parser.print_help()
        sys.exit(0)

    if options.approx:options.graph=False #suppresses output for graph
        
    return options



#================================================================================#  
# Methods for parsing newick trees
#================================================================================#  

def get_splits(newick_string,term=False): #if term: append also terminal splits and branch lengths    
    return splits_for_tree(parse_newick(newick_string),term)

def get_splits_with_decomposition(newick1,newick2,term=False): #if term: append also terminal splits and branch lengths
    global opts
    
    root1=parse_newick(newick1)
    root2=parse_newick(newick2)

    
    taxad,decomp=root1.decomposition(root2) #returns taxa-dictionary of dummy taxa and decompositions as pairs of root nodes

    all_splits=[] #stores pairs of subtrees to be compared

    for droot1,droot2 in decomp:
        act_splits=[]

        for droot in droot1,droot2:
            act_splits.append(list(splits_for_tree(droot,term)))
        all_splits.append(act_splits)
                
    return taxad,all_splits

def get_split_representation(splits1,splits2): #extract splits only in one tree and compute adjacency matrix
    
    def is_compatible(s1,s2):
        s1r,s1l=s1.split('|') #right and left side
        s2r,s2l=s2.split('|')
        for i in s1r,s1l:
            act1=set(i.split('*'))
            for j in s2r,s2l:
                act2=set(j.split('*')) 
                if len(act1.intersection(act2))==0:return True
        return False

    splits1s=set(splits1).difference(set(splits2)) #splits only in T1
    splits2s=set(splits2).difference(set(splits1)) #splits only in T2

    splits=list(splits1s)+list(splits2s) #now the splits are only those which occur in Exactly one tree, the first dim1 are from T1 and the last dim2 are from T2
    
    dim1=len(splits1s)
    dim2=len(splits2s)
    
    adj=[]
    for i in range(0,dim1):
        adj.append([False]*dim2)
        
    for i in range(0,dim1):
        for j in range(0,dim2):
            c=is_compatible(splits[i],splits[dim1+j])
            adj[i][j]=c
    
    return splits, adj, len(splits1s), len(splits2s)


#================================================================================#  
#Methods for finding the geodesic path
#================================================================================#  

def geodesic(adj,bl1,bl2,neg,todo): #returns the last orthant

    bestp=[None]
    mind=[cone(bl1,bl2)+1]
    countp=[0]
    
    def best_path(edge,past=[]): #finds shortest path through graph
        past.append(edge)
        if not edge.anc.edges:
            countp[0]+=1
            newp=Path(past)
            dist=newp.distance(bl1,bl2)
            if sum(dist)<mind[0]:
                mind[0]=sum(dist)
                bestp[0]=newp
        for e in edge.anc.edges:
            if e.s < past[-1].s:best_path(e,deepcopy(past))
            
    def one_path(edge):
        elist=[edge]
        while edge.anc.edges:
            edge=edge.anc.edges[0]
            elist.append(edge)
        return Path(elist)

    def next_sub(str,i): return [x for (pos,x) in zip(list(range(len(str))), str) if (2**pos) & i] #enumerate all 2**str subsets by binary numbers

    #set class variables
    Orthant.opt=opt
    Orthant.dim1=len(bl1)
    Orthant.dim2=len(bl2)
    Orthant.adj=adj

    length=2**Orthant.dim2  #number of indices generated for the orthants
    Orthant.OrthList=[None]*length

    ###Orthant.IDList=[-1]*length #stores mapping of simple IDs (binary numbers) to complex IDs (binomial coefficients)
    Orthant.IDList={}
    
    start=Orthant(Neg=neg,todoPos=todo)
    startid=start.get_id()
    Orthant.OrthList[startid]=start

    max_i=0 #index of last orthant

    for i in range(startid,length):
        if not Orthant.OrthList[i]:continue #many orthants are not generated bec suborthants of a larger one
        max_i=i
        orth=Orthant.OrthList[i]

        actdim=len(orth.get_todo())
    
        for j in range(1,2**actdim): #0. Subset is empty
            switch=next_sub(orth.get_todo(),j) 
            newedge=orth.clone(switch) #switch is the initial R, the complete edge is generated in the constructor

            oid=newedge.get_id()
            oid = int(oid) if int(oid) == oid else oid
            oldorth=Orthant.OrthList[oid]
            if oldorth and oldorth.edges and (newedge in oldorth.edges):
                del newedge
                continue  #path already computed

            if not oldorth: #generate orthant where the edge points to
                neworth=newedge.create_orthant()
                succ=newedge.compute_s(bl1,bl2,neworth)
                if succ: Orthant.OrthList[oid]=neworth #may be not successfull, because transition times have to fullfill several constraints
                else:
                    del newedge
                    del neworth
                
            else:
                succ=newedge.compute_s(bl1,bl2,oldorth)
                if not succ:continue
                if graph:oldorth.edges.append(newedge)
                else: 
                    oldedge=oldorth.edges[0]
                    maxsk=max(newedge.s,oldedge.s)
                    distdiff=sum(one_path(newedge).distance(bl1,bl2))-sum(one_path(oldedge).distance(bl1,bl2))
                    if distdiff<0: oldorth.edges=[newedge]
                    
    lastedge=IEdge(0,0,0,0,anc=Orthant.OrthList[max_i],s=1)
    best_path(lastedge)

    return bestp[0],countp[0]

    


#================================================================================#
#Methods for computing all the distances
#================================================================================#  

def cone(diff1,diff2,shared1=[],shared2=[]):
    sharednorm=snorm([shared1[i]-shared2[i] for i in range(0,len(shared1))])
    return math.sqrt((norm(diff1)+norm(diff2))**2+sharednorm)


def distance(tree1,tree2):

    def permutation_indices(data): 
        return sorted(list(range(len(data))), key = data.__getitem__)

    def branch_score(diff1,diff2,shared1=[],shared2=[]):
        shared=[shared1[i]-shared2[i] for i in range(0,len(shared1))]
        return norm(diff1+diff2+shared)

    def combine(diff,shared,splits1,bl1,dstart,dend): #combine branch length lists so that indicees correspond to splits
        shared_branch=[bl1[splits1.index(s)] for s in shared]
        diff_branch=[bl1[splits1.index(diff[i])] for i in range (dstart,dend)]
        return diff_branch,shared_branch
        
    def create_shared_equ(branch1,branch2):
        equs=[]
        for i in range(0,len(branch1)):
            equs.append([0,branch2[i]-branch1[i],branch1[i]])
        return equs
    
    def inverse(mat):return [[mat[i][j] for i in range(0,len(mat))] for j in range(0,len(mat[0]))]

    #================================================================================#  

    outstring=""


    if decomp:
        taxad,split_decomp=get_splits_with_decomposition(tree1,tree2,term)
        if len(split_decomp)==1:
            header=False
            spp=split_decomp[0][0][2]
        else:
            header=True
            spp=get_splits(tree1,term)[2]
    else:
        split_decomp=[[list(get_splits(tree1,term)),list(get_splits(tree2,term))]]
        header=False
        spp=split_decomp[0][0][2]
        
    all_equs=[]
    dec=1 #number of decomposition

    if outfile:
        outfile.write("Output of GeoMeTree v1.0\n\n")
        outfile.write("2 Trees of %u taxa given:\n"%spp)
        outfile.write("T1=%s\nT2=%s\n" %(tree1,tree2))
        if normalize:outfile.write("\nTree vectors have been normalized to norm 1 !!!\n")
        if header:
            outfile.write("\nTrees have been decomposed at common splits, the following dummy taxa are used:\n")
            for t in sorted(taxad.keys()):
                outfile.write("%s\t%s\n" %(t,"*".join(taxad[t])))

    if normalize:
        for j in [0,1]: #first and second tree
            bl=[]
            for s in split_decomp:
                bl+=s[j][1]
            n=norm(bl)
            for s in split_decomp:
                s[j][1]=[x/n for x in s[j][1]]

    s_dic={} #saves geodesic path in transition times and Left and Right sets
    bl_dic1={} #saves branch length of splits in case of multiple decompositions
    bl_dic2={}
    len_list=["cone","bs","geod","coneall","bsall","geodall"]
    split_list=["shared","diff1","diff2","dim"]
    len_dic={}.fromkeys(len_list) #save length and splits to compute overall numbers in the end
    for l in list(len_dic.keys()):
        len_dic[l]=[0]
    split_dic={}.fromkeys(split_list,0)
    graph=1
    compl_time=0

    for t1,t2 in split_decomp:
        splits1,bl1,spp1=t1
        splits2,bl2,spp2=t2
        new_s_dic={}
        
        # =============================================================================== #
        # creates set S (diff_splits), compatibility matrix (adj), C (shared_splits) and corresponding numbers
 
        diff_splits,adj,dim1,dim2=get_split_representation(splits1,splits2)
    
        shared_splits=list(set(splits1).intersection(set(splits2)))
        shared_splits.sort()
        c=len(shared_splits)
        overlap=set(splits1).intersection(splits2)

        # =============================================================================== #
        # combine branch length lists so that indices correspond to splits
        
        branch1_diff,branch1_shared=combine(diff_splits,shared_splits,splits1,bl1,0,dim1)
        branch2_diff,branch2_shared=combine(diff_splits,shared_splits,splits2,bl2,dim1,dim1+dim2)
    
        # =============================================================================== #
        # output

        outstring+="%u\t%u\t%u\t%u\t"%(spp1,len(overlap),len(branch1_diff),len(branch2_diff))
        split_dic["shared"]+=len(overlap)
        split_dic["diff1"]+=len(branch1_diff)
        split_dic["diff2"]+=len(branch2_diff)
 
        if outfile:
            outfile.write("\n%s\n"%('-'*80))
            if header:
                outfile.write ("\nResults for Decomposition No. %u:\n" %dec)

            outfile.write("\nSplits only in T1:\n")
            for i in range(0,dim1):
                splt=shorter_split(diff_splits[i])
                if header:outfile.write("%u/"%dec)
                outfile.write("%u\t%s\t%1.6f\n" %(i+1,splt,branch1_diff[i]))

            outfile.write("\nSplits only in T2:\n")
            for i in range(0,dim2):
                splt=shorter_split(diff_splits[dim1+i])
                if header:outfile.write("%u/"%dec)
                outfile.write("%u\t%s\t%1.6f\n" %(i+dim1+1,splt,branch2_diff[i]))    

            outfile.write("\nSplits common to both trees and branch length in T1 and T2:\n")
            for i in range (0,c):
                splt=shorter_split(shared_splits[i])
                if header:outfile.write("%u/"%dec)
                outfile.write("%u\t%s\t%1.6f\t%1.6f\n" %(i+dim1+dim2+1,splt,branch1_shared[i],branch2_shared[i]))
                
        # =============================================================================== #
        # fill bl_dics
        if header:
            for i in range(0,dim1):
                bl_dic1["%u/%u"%(dec,i+1)]=branch1_diff[i]
                bl_dic2["%u/%u"%(dec,i+1)]=0
            for i in range(0,dim2):
                bl_dic1["%u/%u"%(dec,dim1+i+1)]=0
                bl_dic2["%u/%u"%(dec,dim1+i+1)]=branch2_diff[i]
            for i in range (0,c):
                bl_dic1["%u/%u"%(dec,dim1+dim2+i+1)]=branch1_shared[i]
                bl_dic2["%u/%u"%(dec,dim1+dim2+i+1)]=branch2_shared[i]
        else:
            for i in range(0,dim1):
                bl_dic1[str(i+1)]=branch1_diff[i]
                bl_dic2[str(i+1)]=0
            for i in range(0,dim2):
                bl_dic1[str(dim1+i+1)]=0
                bl_dic2[str(dim1+i+1)]=branch2_diff[i]
            for i in range (0,c):
                bl_dic1[str(dim1+dim2+i+1)]=branch1_shared[i]
                bl_dic2[str(dim1+dim2+i+1)]=branch2_shared[i]
            
                
        # =============================================================================== #
        # there are different splits
    
        if dim1 or dim2:
            # =============================================================================== #
            # compute bounds          
            
            bs=branch_score(branch1_diff,branch2_diff) #lower bound
            ub=cone(branch1_diff,branch2_diff) #upper bound
            len_dic["cone"].append(ub)
            len_dic["bs"].append(bs)
            if approx: #do not compute geodesic path
                if outfile:
                    outfile.write("\nResults for the different splits:\n")
                    outfile.write("Branch score %1.6f\nCone distance %1.6f\n" %(bs,ub))
                outstring += "0\t%1.6f\t%1.6f\t0\t" %(bs,ub)
                count=0
                timeg=0
            else:
          
                # =============================================================================== #
                # find splits that are compatible to all others

                l_ind=[]
                #There may be some splits in the first tree that are compatible to all splits in the second tree -> they have to end at 1
                for i in range(0,dim1):
                    if adj[i]==[True]*dim2:l_ind.append(i)
                comp_equ_l=[[1,-branch1_diff[i],branch1_diff[i]] for i in l_ind]
                neg=set(range(0,dim1)).difference(set(l_ind))  #neg is the starting neg for the geodesic algorithm, in case of full compatibilities, some are excluded

                r_ind=[]
                #vice versa
                for j in range(0,dim2):
                    comp=True
                    for i in range(0,dim1):
                        if not adj[i][j]:
                            comp=False
                            break
                    if comp:r_ind.append(j)
                comp_equ_r=[[0,branch2_diff[i],0] for i in r_ind]
                todo=set(range(0,dim2)).difference(set(r_ind))        
        
                # =============================================================================== #
                # geodesic distance algorithm if still something todo

                t0=time.clock()
                if len(todo)>0: #it may be that all were compatible bec of polytomies

                    # =============================================================================== #
                    # todo should be the smaller set since algorithm is exponential in len(todo)
                    swap=(len(neg)<len(todo))
                    if swap:
                        adj=inverse(adj)
                        branch1_diff,branch2_diff=branch2_diff,branch1_diff
                        neg,todo=todo,neg
                    
                    outstring+="%u\t" %len(todo)
                    split_dic["dim"]+=len(todo)

                    try:
                        path,count=geodesic(adj,branch1_diff,branch2_diff,neg,todo)
                    except OverflowError:
                        print("To many splits in actual decomposition to compute the geodesic distance exactly:", Orthant.dim2)
                        print("Start the computation with option '-a' to compute the approximations")
                        continue
                    
                    t1=time.clock()
               
                    equ_l,equ_r= path.equations(branch1_diff,branch2_diff) #the assignment to first and second tree is not correct, but changes nothing for distance computation!!!
                    if swap:
                        branch1_diff,branch2_diff=branch2_diff,branch1_diff
                        path.inverse()
            
                else:
                    equ_l=[]
                    equ_r=[]
                    count=1
                    path=None
                    t1=t0
                    outstring+="0\t"

                # =============================================================================== #
                # Get the actual L, R and s and save in s_dic
                
                if path:
                    L = path.get_L();L=[[i+1 for i in l] for l in L]
                    R = path.get_R();R=[[i+dim1+1 for i in r] for r in R]
                    s = path.get_s()
                else:
                    L=[[],[]]
                    R=[[],[]]
                    s=[0,1]
                if r_ind: R[0]+=[i+dim1+1 for i in r_ind]
                if l_ind: L[-1]+=[i+1 for i in l_ind]    
                if header:
                    L=[["%u/%u"%(dec,i) for i in l] for l in L]
                    R=[["%u/%u"%(dec,i) for i in r] for r in R]
 
                for i in range(0,len(s)):
                   new_s_dic[s[i]]=[L[i],R[i]]


                # =============================================================================== #
                # Compute the distances for different splits and output them

                equ_l+=comp_equ_l
                equ_r+=comp_equ_r

                dist=dist_for_geod(equ_r,equ_l)
                timeg=t1-t0
                outstring+="%1.6f\t%1.6f\t%1.6f\t"%(bs,ub,sum(dist))

                    
                if outfile:
                    outfile.write("\nThe geodesic path:\n")
                    outfile.write(path_to_str(s,L,R))
                    if header:outfile.write("\nTransition points:\n%s\n"%trans_points(new_s_dic,bl_dic1,bl_dic2,dec))
                    else:outfile.write("\nTransition points:\n%s\n"%trans_points(new_s_dic,bl_dic1,bl_dic2))
                    outfile.write("\nResults for the different splits:\n")
                    outfile.write("Branch score %1.6f\nGeodesic distance %1.6f\nCone distance %1.6f\n" %(bs,sum(dist),ub))
                    
                len_dic["geod"].append(sum(dist))
            
        # =================================================================================== #
        # There are no different splits
        
        else:
            timeg=0
            equ_r=[];equ_l=[];count=1 #count is number of geodesic paths
            outstring+="0\t"*4
            new_s_dic[0]=[[],[]]
            new_s_dic[1]=[[],[]]
            if outfile:
                outfile.write("\nThe geodesic path:\nThere are no different splits, it equals the branch score.\n")
                if header:outfile.write("\nTransition points:\n%s\n"%trans_points(new_s_dic,bl_dic1,bl_dic2,dec))
                else:outfile.write("\nTransition points:\n%s\n"%trans_points(new_s_dic,bl_dic1,bl_dic2)) 

        # =============================================================================== #
        # Compute the distances for all splits and output them
        
        bs_compl=branch_score(branch1_diff,branch2_diff,branch1_shared,branch2_shared)
        cone_compl=cone(branch1_diff,branch2_diff,branch1_shared,branch2_shared)
        len_dic["coneall"].append(cone_compl)
        len_dic["bsall"].append(bs_compl)
        
        if approx:
            outstring+="%1.6f\t%1.6f\t0\t0\t"%(bs_compl,cone_compl)
            if outfile:
                outfile.write("\nResults for all splits:\n")
                outfile.write("Branch score %1.6f\nCone distance %1.6f\n" %(bs_compl,cone_compl))

        else:
            dist_compl=dist_for_geod(equ_r+create_shared_equ(branch1_shared,branch2_shared),equ_l)
            outstring+="%1.6f\t%1.6f\t%1.6f\t%1.6f\t"%(bs_compl,cone_compl,sum(dist_compl),timeg)
            if graph: outstring+="%u\t"%count
        
            if outfile:
                outfile.write("\nResults for all splits:\n")
                outfile.write("Branch score %1.6f\nGeodesic distance %1.6f\nCone distance %1.6f\n" %(bs_compl,sum(dist_compl),cone_compl))
                if graph: outfile.write("\nComplete graph enumerated: %u path(s) found!\n" %count)
                outfile.write("\nComputation time for geodesic path: %1.4f s\n"%(timeg))

            len_dic["geodall"].append(sum(dist_compl))


        for s in list(new_s_dic.keys()):
            if s in s_dic:
                s_dic[s][0]+=new_s_dic[s][0]
                s_dic[s][1]+=new_s_dic[s][1]
            else: s_dic[s]=new_s_dic[s]
  
        compl_time+=timeg
        if graph:graph*=count
        dec+=1

    #===============================================================================#
    #output in case of decompositions

    if header:
        for i in list(len_dic.keys()):
            len_dic[i]=norm(len_dic[i])
        outstring2=("%u\t"*5+"%1.6f\t"*7) %(spp,split_dic["shared"],split_dic["diff1"],split_dic["diff2"],split_dic["dim"],len_dic["bs"],len_dic["cone"],len_dic["geod"],len_dic["bsall"],len_dic["coneall"],len_dic["geodall"],compl_time)
        if graph:outstring2+="%u\t"%graph
        outstring=outstring2+outstring
        
        if outfile:
            outfile.write("\n%s\n"%('-'*80))
            if not approx:
                outfile.write("\nThe complete geodesic path:\n%s"%path_to_str(s_dic))
                outfile.write("\nTransition points:\n%s"%trans_points(s_dic,bl_dic1,bl_dic2))
            else:outfile.write("\nComplete results:\n")
            outfile.write("\nResults for the different splits:\n")
            outfile.write("Branch score %1.6f\n"%len_dic["bs"])
            if not approx:outfile.write("Geodesic distance %1.6f\n"%len_dic["geod"])
            outfile.write("Decomposed cone distance %1.6f\n" %len_dic["cone"])
            outfile.write("\nResults for all splits:\n")
            outfile.write("Branch score %1.6f\n"%len_dic["bsall"])
            if not approx:outfile.write("Geodesic distance %1.6f\n"%len_dic["geodall"])
            outfile.write("Decomposed cone distance %1.6f\n" %len_dic["coneall"])
            if graph: outfile.write("\nComplete graphs enumerated: with independent decompositions %u path(s) found!\n" %graph)
            if not approx:outfile.write("\nComputation time for all geodesic paths: %1.4f s\n"%(compl_time))
    
    # return outstring
    return len_dic["geodall"][1] if type(len_dic["geodall"]) is list else len_dic['geodall']

#===============================================================================#

def main():
    global opts
    opts=parse_options()
    trees=open(opts.infile).readlines()

    if opts.outfile=='None': outfile=None
    elif opts.outfile=='0': outfile=open(opts.infile+'.dist',"w")
    else: outfile=open(opts.outfile,'w')

    for i in range(0,len(trees)):
        for j in range(i+1,len(trees)):
            if opts.silent: actpair=None
            else: actpair=open("%s_%u_%u" %(opts.header,i+1,j+1),"w")
            outstr=distance(trees[i][:-1],trees[j][:-1],actpair)
            if outfile: outfile.write("%u\t%u\t%s\n"%(i+1,j+1,outstr))
            if actpair: actpair.close()
    if outfile: outfile.close()
    
if __name__ == "__main__":main()
