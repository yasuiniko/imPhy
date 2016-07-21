class Node(object):
    def __init__(self,anc,bl=1):
        self.ancestor=anc                   # ancestor Node
        self.children=[]                    # children Nodes
        self.branch=bl                      # branch length of branch leading to actual node
        self.name=""                        # name of node, eg. taxa name
        self.depth=len(self.ancestors())-1  # depth of Node, root has depth 0

    def set_child(self,child): self.children.append(child) #append child to list

    def ancestors(self): #return ancestor nodes and self node in order root->bottom
        if self.ancestor==None: return [self.name]
        else: return self.ancestor.ancestors()+[self.name]

    def traverse(self): #returns self node + children traversing
        trav = [self]
        for child in self.children:
            trav+=child.traverse()
        return trav

    def get_leaves(self): #returns list of leaves under node
        if len(self.children)==0:return [self.name]
        l=[]
        for child in self.children:
            l+=child.get_leaves()
        return l       

    def newick(self,intern=False,bl=False): #creates newick string from tree
        if len(self.children)==0:
            if bl: return "%s:%f" %(self.name,self.branch)
            else: return self.name
        s="("+self.children[0].newick(intern,bl)
        for i in range(1,len(self.children)):
            s+=","+self.children[i].newick(intern,bl)
        s+=')'
        if self.ancestor==None:return s
        if intern:s+=self.name
        if bl:s+=":%f" %self.branch
        return s


    def __split_representation(self,all_leaves,act_leaves): #+/- representation which splits are present
        n=len(all_leaves)
        act_string=''
        for i in range (0,n):
            if (all_leaves[i] in act_leaves):act_string+='+'
            else: act_string+='-'
        return act_string
    
    def decomposition(self,other): #create decompositions on common edges of self and other --> returns list of pairs
        all_trees=[[self,other]]
        nodedic={}
        all_leaves=self.get_leaves()
        all_leaves.sort()
        n=len(all_leaves)
        count=1
        taxadic={}
        
        for node in other.traverse():
            act_leaves=set(node.get_leaves())
            nodedic[self.__split_representation(all_leaves,act_leaves)]=node
        for node in self.traverse():
            act_leaves=set(node.get_leaves())
            if len(act_leaves)==n or len(act_leaves)==1:continue
            act_string=self.__split_representation(all_leaves,act_leaves)
            if act_string in nodedic:
                othernode=nodedic[act_string]
                otherchild=set(othernode.children[0].get_leaves())
                distinct=True
                for child in node.children: #are only 2
                    if set(child.get_leaves())==otherchild:distinct=False
                    
                if distinct: #decompose!
                    newroots=[]
                    taxa=node.get_leaves()
                    taxa.sort()
                    taxadic["d%i"%count]=taxa
                    othert=list(set(all_leaves).difference(set(taxa)))
                    othert.sort
                    taxadic["d%i"%(count+1)]=othert
                    
                    for act_node in node,othernode:
                        
                        newnode=Node(act_node.ancestor,act_node.branch)
                        newnode.name="d%i"%count #dummy node
                        ac=act_node.ancestor.children
                        change=0
                        if ac[1]==act_node:change=1
                        ac[change]=newnode

                        act_node.branch=0
                        newroot=Node(None,0)
                        newroot.name='root'
                        newroot.set_child(act_node)
                        newnode=Node(newroot,0)
                        newnode.name="d%i"%(count+1)
                        newroot.set_child(newnode)
                        act_node.ancestor=newroot
                        newroots.append(newroot)
                    count+=2
                    
                    all_trees.append(newroots)
        return taxadic,all_trees

    def get_bl(self,term=True): #return branch lengths list
        bl=[]
        for node in self.traverse():
            if not node.ancestor:continue
            if not term and len(node.children)==0:continue
            bl.append(node.branch)
        return bl

    def get_graph(self): #only topology, returns adjacency matrix
        nodes=self.traverse()
        num=len(nodes)
        
        leaves=nodes[0].get_leaves()
        splits=[nodes[0].__split_representation(leaves,leaves)]
        namedic={nodes[0].name:0}
        
        adj=[]
        for i in range(0,num):
            adj.append([0]*num)
        for i in range(1,num):
            anc=nodes[i].ancestor
            adj[i][namedic[anc.name]]=1;adj[namedic[anc.name]][i]=1;
            for c in anc.children:
                if c.name in namedic: adj[i][namedic[c.name]]=1;adj[namedic[c.name]][i]=1; #else it is set when brother is processed
            splits.append(nodes[i].__split_representation(leaves,nodes[i].get_leaves()))
            namedic[nodes[i].name]=i
        return adj,splits,leaves
            
