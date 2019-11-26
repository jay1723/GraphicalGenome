import cjson
import numpy
import json
import re
import time

class Nodes(dict):

    def __init__(self, inp):
        """Initialize Nodes object
        def __init__(self:<Nodes>, inp:<dict>) -> Void
        """
        self.nodes = inp
    # Override getitem to account for lazy evaluation
    def __getitem__(self,key):
        """Get node annotations tied to specific node
        def __getitem__(self:<Nodes>,key:<str>) -> <dict>
        """
        return self.getNode(key)
    
    def __setitem__(self, key, value):
        """Creates/Updates the value in the Nodes class dictionary specified by the key
        def __setitem__(self:<Nodes>, key:<str>, value:<dict>) -> Void
        """
        self.nodes[key] = value
    
    def __contains__(self, key):
        """Checks to see whether the given key exists in the Nodes class dictionary
        def __contains__(self:<Nodes>, key:<str>) -> <bool>
        """
        return key in self.nodes
            
    def __repr__(self):
        """Defines console representation of object. Defaults to dictionary representation
        def __repr__(self:<Nodes>) -> <dict>
        """
        return self.nodes
    
    def keys(self):
        """Returns an iterable set of the key values for the Nodes class dictionary
        def keys(self:<Nodes>) -> <dictview>
        """
        return self.nodes.keys()
    
    # getNode allows for lazy evaluation
    def getNode(self, nodename):
        """Returns the dictionary entry for the given nodename. If the header has yet to be parsed for the given key then it 
        is evaluated and added to the dictionary.
        def getNode(self:<Nodes>, nodename:<str>) -> <dict>
        """
        if "hdr" in self.nodes[nodename].keys():
            seq = self.nodes[nodename]["seq"]
            hdr = self.nodes[nodename]["hdr"]
            del self.nodes[nodename]
            self.nodes[nodename] = cjson.decode(hdr)
            self.nodes[nodename]["seq"] = seq
            return self.nodes[nodename]
        return self.nodes[nodename]
    
    def removeNode(self, nodename):
        """Remove the given node from the dictionary if it exists. Otherwise return None.
        def removeNode(self:<Nodes>, nodename:<str>) -> <dict or None>
        """
        return self.nodes.pop(nodename, None)
    
class Edges(dict):
    def __init__(self, inp):
        """Initialize the Edges object with the provided dictionary
        def __init__(self:<Edges>, inp:<dict>) -> Void
        """
        self.edges = inp
        
    def __getitem__(self,key):
        """Returns the value in the Edges class dictinoary specified by the key.
        def __getitem__(self:<Edges>, key:<str>) -> <dict>
        """
        return self.getEdge(key)
    
    def __setitem__(self, key, value):
        """Updates/Creates the value for the specified key
        def __setitem__(self:<Edges>, key:<str>, value:<dict>) -> Void
        """
        self.edges[key] = value
    
    def __contains__(self, key):
        """Checks to see if the key exists within the Edges class dictionary
        def __contains__(self:<Edges>, key:<str>) -> <bool>
        """
        return key in self.edges
    
    def __repr__(self):
        """Defines the console representation of the object. Defaults to dictionary console representation
        def __repr__(self:<Edges>) -> <dict>
        """
        return self.edges
            
    def keys(self):
        """Returns an iterable set of the key values for the Edges class dictionary
        def keys(self:<Edges>) -> <dictview>
        """
        return self.edges.keys()            
    
    def getEdge(self, edgename):
        """Returns the dictionary entry for the given edgename. If the header has yet to be parsed for the given key then it 
        is evaluated and added to the dictionary.
        def getEdge(self:<Edges>, edgename:<str>) -> <dict>
        """
        if "hdr" in self.edges[edgename].keys():
            seq = self.edges[edgename]["seq"]
            hdr = self.edges[edgename]["hdr"]
            del self.edges[edgename]
            self.edges[edgename] = cjson.decode(hdr)
            self.edges[edgename]["seq"] = seq
            return self.edges[edgename]
        return self.edges[edgename]    
    
    def removeEdge(self, edgename):
        """Remove the given edge from the dictionary if it exists. Otherwise return None.
        def removeEdge(self:<Edges>, edgename:<str>) -> <dict or None>
        """
        return self.edges.pop(edgename, None)

class GraphicalGenome:
    
    def __init__(self, nodefile, edgefile, verbose=False, enamelen=12, nnamelen=12):
        """Initializes the Graphical Genome data structure. Requires only 2 positional parameters: nodefile and edgefile
        def __init__(self:<GraphicalGenome>, nodefile:<str>, edgefile:<str>, 
                        verbose=False:<bool>, enamelen=12:<int>, nnamelen:<int>) -> Void
        Parameters:
            Nodefile: <str> - Path to the node fasta file
            Edgefile: <str> - Path to the edge fasta file
            verbose: <bool> - True for verbose print statements, False for silent
            enamelen: <int> - Length of the edge key names contained in the edges dictionary
            nnamelen: <int> - Length of the node key names contained in the nodes dictionary
        """
        
         # Global variables
        self.edgenamelength = enamelen
        self.nodenamelength = nnamelen
        
        self.verbose = verbose
        self.nodes, self.edges, self.outgoing, self.incoming, self.max_node_name, self.max_edge_name, self.genes = self.loadChromosomeGraph(nodefile, edgefile)
        self.nodes = Nodes(self.nodes)
        self.edges = Edges(self.edges)
        self.sortednodes = []
        self.subpath = []
        self.subsequence = []
           
    def __repr__(self):
        """Returns the console representation of the Graphical Genome Class. Returns the number of nodes and edges.
        def __repr__(self:<GraphicalGenome>) -> <str>
        """
        return "Number of Nodes %d - Number of Edges %d" % (len(self.nodes), len(self.edges))
        
    def __str__(self):
        """Defines representation of Graphical Genome when using the print statement. Returns the nodefile/edgefile paths
        def __str__(self:<GraphicalGenome>) -> <str>
        """
        return "Origin - Nodefile %s    Edgefile %s" % (nodefile, edgefile)
    
    def encodeASCII(self, data):
        """Encodes the data as a dictionary in ascii format. 
        def encodeASCII(self:<GraphicalGenome>, data:<dict>) -> <dict>
        """
        return dict((k.encode('ascii'), [i.encode('ascii') for i in v] if isinstance(v, list) else v.encode('ascii')) for k, v in data.items())

    

# Inputs:
    # nodefile - Absolute path to the nodefile you wish to open
    # edgefile - Absolute path to the edgefile you wish to open
# Outputs the nodes and edges as dictionaries as well as auxilary dictionaries containing node -> list(edge) pairs 
    def loadChromosomeGraph(self, nodefile, edgefile):
        """Reads in the nodefile/edgefile fasta files and parses them into the GraphicalGenome class objects
        def loadChromosomeGraph(self:<GraphicalGenome>, nodefile:<str>, edgefile:<str>) -> <dict,dict,dict,dict,int,int,dict>
        
        Parameters:
            nodefile: <str> - Path to the node FASTA file
            edgefile: <str> - Path to the edge FASTA file
           
        Output:
            nodes: <dict> - Nodes dictionary
            edges: <dict> - Edges dictionary
            sources: <dict> - Dictionary of the lists of source nodes with the edgenames as keys
            destinations: <dict> - Dictionary of the lists of destination nodes with edgenames as keys
            max_node: <int> - The maximum node number found during parsing
            max_edge: <int> - The maximum edge number found during parsing
            
            Not in use:
            genes: <dict> - Dictionary containing list of nodes on which a given node is found. Gene names are the keys
        """
        # Instantiate returned datastructures
        max_node = 0
        max_edge = 0
        nodes = {}
        edges = {}
        sources = {}
        destinations = {}
        genes = {}
        
        snh = time.time()
        # Node iterations
        with open(nodefile, "r") as fp:
            hdr = fp.readline()[1:]
            seq = fp.readline()
            while hdr and seq:
                node =  hdr[0:self.nodenamelength]
                nodehdr = hdr[self.nodenamelength+1:]
                nodes[node] = {}
                # Strip newlines if they exist at the end
                if nodehdr[-1] == "\n":
                    nodehdr = nodehdr[:-1]
                if seq[-1] == "\n":
                    seq = seq[:-1]
                nodes[node]["hdr"] = nodehdr 
                nodes[node]["seq"] = seq 

                hdr = fp.readline()[1:]
                seq = fp.readline()

#             # Check existence of gene annotation and add to dictionary
#             gene_full = re.search(gene_pattern, hdr)
#             if gene_full:
#                 gene = gene_full.group()
#                 genelist = gene[8:-1].split(",")
#                 for i in genelist:
#                     genekey = i[1:-1]
#                     genes[genekey] = genes.get(genekey, []) + [node]

        if self.verbose:
            print "Node load %5.2f" % (time.time() - snh)
        
        seh = time.time()        
        # Edge iteration
        with open(edgefile, "r") as fp:
            # Calculating the positions of the source and dst nodes
            sstart = self.edgenamelength + 9
            send = sstart + self.nodenamelength
            dstart = send + 9
            dend = dstart + self.nodenamelength
            
            # Read in the data and parse 
            hdr = fp.readline()[1:]
            seq = fp.readline()
            while hdr and seq:
                # Extract the src, dst, key from the header
                src = hdr[sstart:send]
                dst = hdr[dstart:dend]
                key = hdr[0:self.edgenamelength]
                
                # Source and sink nodes are not the same size as other nodes so special case needed
                if "SOURCE" in src:
                    # "SOURCE" is 6 letters long so change the values of the indices to reflect this
                    send = sstart + 6
                    dstart = send + 9
                    dend = dstart + self.nodenamelength
                    dst = hdr[dstart:dend]
                    src = "SOURCE"
                if "SINK" in dst:
                    dst = "SINK"
                
                # Add header and seq to the edges dictionary
                # Strip newlines if they exist at the end
                if hdr[-1] == "\n":
                    hdr = hdr[:-1]
                if seq[-1] == "\n":
                    seq = seq[:-1]
                edges[key] = {}
                edges[key]["hdr"] = hdr[self.edgenamelength+1:]
                edges[key]["seq"] = seq
                
                # Lazy eval the sources and destinations dictionary
                sources[src] = sources.get(src, []) + [key]
                destinations[dst] = destinations.get(dst, []) + [key]

                # Update global counters
                if "F" in key and "S" not in key and "K" not in key:
                    max_edge = max(max_edge, int(key[-7:]))
                if "F" == dst[0]:
                    max_node = max(max_node, int(dst[-8:]))
                if "F" == src[0]:
                    max_node = max(max_node, int(src[-8:]))
                    
                # Load the next line for the next iteration
                hdr = fp.readline()[1:]
                seq = fp.readline()  
                
        if self.verbose:
            print "Edge load %5.2f" % (time.time() - seh)
            print "# Nodes %5.2d" % (len(nodes))
            print "# Edges %5.2d" % (len(edges))
        return nodes, edges, sources, destinations, max_node, max_edge, genes
    
    def writeFasta(self, filename, input_dict, keylist=["src", "dst"]):
        """Write the given node or edge file as a FASTA file. Overwrites existing file. Will create file if it doesn't exist. 
        def writeFasta(self:<GraphicalGenome>, filename:<str>, input_dict:<dict>, keylist=["src","dst"]:<list[str]>) -> Void
        
        Parameters:
            filename: <str> - absolute path for the file you wish to write. 
            input_dict: <dict> - Node or Edge dictionary you wish to write. 
            keylist: <list[str]> - list of strings of dictionary keys that you want to ignore during write. 
                    This will require you to write these values on your own or just discard them. 
        """
        sorted_keys = sorted(input_dict.keys()) 
        with open(filename, "w+") as fastafile:
            # If iterating through the edges, write the edges in the correctly ordered format
            if (sorted_keys[0][0] == "E"):
                for edge in sorted_keys:
                    # If header has not been evaluated, just re-write the header wholesale without any analysis
                    if "hdr" in input_dict[edge].keys():
                        line = ">" + edge + ";" + input_dict[edge]["hdr"] + "\n"
                        line += input_dict[edge]["seq"] + "\n"
                        continue
                    line = ">" + edge + ";{" 
                    # Source
                    line += '"src":"' + input_dict[edge]["src"] + '",'
                    # Destination
                    line += '"dst":"' + input_dict[edge]["dst"] + '"'
                    for key in input_dict[edge].keys():
                        if key == "seq":
                            continue
                        if key in keylist:
                            continue
                        line += ',"' + key + '":' + json.dumps(input_dict[edge][key], separators=(",", ":"))
                    line += "}\n"
                    line += input_dict[edge]["seq"] + "\n"
                    fastafile.write(line)
            # If iterating over nodes, just write the nodes normally
            else:
                for i in sorted_keys:
                    line = ">" + i + ";"
                    obj = {}
                    for j in input_dict[i].keys():
                        if j == 'seq':
                            continue
                        obj[j] = input_dict[i][j]
                    line += json.dumps(obj, separators=(",", ":"))
                    line += "\n" + input_dict[i]['seq'] + "\n"
                    fastafile.write(line)
        

    def returnNodeEdgePath(self, strain, chromo):
        """Returns a list of [node,edge,node,edge...] that traces out a path for a strain on a given chromo
        def returnNodeEdgePath(self:<GraphicalGenome>, strain:<str>, chromo:<str>) -> list[str]
        
        Parameters:
            strain: <str> - founder or CC strain you wish to trace. 
            chromo: <str> - chromosome you are iterating through
        """
        edges = self.edges
        nodes = self.nodes
        outgoing = self.outgoing
        chrToName = {}
        chromoDict = {"Chr1":"01","Chr2":"02","Chr3":"03","Chr4":"04","Chr5":"05", "Chr6":"06", "Chr7":"07", "Chr8":"08",
             "Chr9":"09", "Chr10":"10", "Chr11":"11", "Chr12":"12", "Chr13":"13", "Chr14":"14", 
              "Chr15":"15", "Chr16":"16", "Chr17":"17", "Chr18":"18", "Chr19":"19", "ChrX":"20"}
        path = []
        # Source requires haplotype, initialize empty constructed sequence
        source = "S" + chromoDict[chromo] + "." + strain

        # Iterate through the edges in order based on the destinations and outgoing edges from those destinations
        currentEdge = source
        while True:
            dst = edges[currentEdge]["dst"]
            path.append(currentEdge)
            if dst == "SINK":
                return path
            else:
                path.append(dst)
                for edge in outgoing[dst]:
                    if strain in edges[edge]["strain"]:
                        currentEdge = edge 
        return path 
    

# 
    def findMissingPaths(self, strain, chromo):
        """Given a CC path, look for all node pairs that do not contain an edge between them that has the specified CC strain on them.
        def findMissingPaths(self:<GraphicalGenome>, strain:<str>, chromo:<str>) -> <list[tuple(node,node)]>
        
        Parameters:
            strain: <str> - Strain you wish to trace 
            chromo: <str> - Chromosome you are tracing the strain through 
        
        Notes:
            THIS CURRENTLY WORKS FOR EDGES WITHOUT FLOATING NODES, WHEN FLOATING NODES ARE ADDED THIS NEEDS TO BE UPDATED
        """
        nodePairs = set()
        visited = set()
        for i in self.nodes:
            sawStrain = False
            for j in self.outgoing[i]:
                if strain in self.edges[j]["strain"]:
                    sawStrain = True
            for k in self.incoming[i]:
                if strain in self.edges[k]["strain"]:
                    sawStrain = True
            if sawStrain:
                continue
            else:
                visited.add(i)
        for edge in self.edges:
            if self.edges[edge]["src"] in visited and self.edges[edge]["dst"] in visited:
                nodePairs.add((self.edges[edge]["src"], self.edges[edge]["dst"]))
        return nodePairs

    # Reconstruct an entire sequence for a given strain across the entire genome
    def reconstructSequence(self, strain, path=0):
        """Reconstruct the sequence given a strain (ABCDEFGH or CC) and can return more than one sequence if heterozygosity is known
        def reconstructSequence(self:<GraphicalGenome>, strain:<str>, path:<int>) -> <str>
                
        Parameters:
            strain: <str> - ABCDEFGH or CC strain you wish to trace
            path: <int> - 0 for base path, 1 for first heterozygous path identified, 2 for second heterozygous path identified. 
        """
        nodes = self.nodes
        edges = self.edges
        outgoing = self.outgoing
        # Source requires haplotype, initialize empty constructed sequence
        source = ""
        het = strain
        if path == 0:
            het = strain + 'a'
        elif path == 1:
            het += "a"
        else:
            het += "b"

        for src in outgoing["SOURCE"]:
            for edge in edges[src]["strain"]:
                if het in edge or strain in edge:
                    source = src
        if source == "":
            print het + "strain not found on any source path"
            return ""
        conseq = ""

        # Iterate through the edges in order based on the destinations and outgoing edges from those destinations
        currentEdge = source
        firstNode = edges[source]["dst"]
        founder_on_that_Edge = list(set(edges[source]['strain']) & set('ABCDEFGH'))[0]
        numberN = int(nodes[firstNode][founder_on_that_Edge]) - 1 - len(edges[source]["seq"])
        conseq += "N" * numberN
        while True:
            dst = edges[currentEdge]["dst"]
            if dst == "SINK":
                conseq += edges[currentEdge]["seq"]
                conseq += ("N" * int(edges[currentEdge]["addNs"]))
                return conseq.upper()
            else:
                if "F" in dst or 'B' in dst:
                    conseq += edges[currentEdge]["seq"]
                else:
                    conseq += edges[currentEdge]["seq"] + nodes[dst]["seq"]
                
                elist = outgoing[dst]
                if len(elist) == 1:
                    currentEdge = elist[0]
                else:
                    for edge in outgoing[dst]:
                        if strain in edges[edge]["strain"] or het in edges[edge]["strain"]:
                            currentEdge = edge


    # Used in the CLI, writes a sequence to a specified file
    def writeSequence(self, strain, filedest, path=0):
        """CLI function that writes a specified sequence to a file. 
        def writeSequence(self:<GraphicalGenome>, strain:<str>, path:<int>):
        
        Parameters:
            strain: <str> - ABCDEFGH or CC strain you wish to trace
            path: <int> - 0 for base path, 1 for first heterozygous path identified, 2 for second heterozygous path identified. 
        """
        reconstructed = self.reconstructSequence(strain, path)
        with open(filedest, "w+") as seqfile:
            seqfile.write(reconstructed)
            
    # Given an anchor, provide the next anchor
    def nextAnchor(self, node):
        """Return the next monotonically incrementing anchor given a nodename in the graph
        def nextAnchor(self:<GraphicalGenome>, node:<str>) -> <str>
        """
        foundAnchor = False
        currentEdge = self.outgoing[node][0]
        currentNode = self.edges[currentEdge]["dst"]
        
        while not foundAnchor:
            if currentNode[0] == "A" or currentNode == "SINK":
                foundAnchor = True
            else:
                currentEdge = self.outgoing[currentNode][0]
                currentNode = self.edges[currentEdge]["dst"]
        return currentNode

    def prevAnchor(self, node):
        """Return the prev monotonically incrementing anchor given a nodename in the graph
        def prevAnchor(self:<GraphicalGenome>, node:<str>) -> <str>
        """
        foundAnchor = False
        currentEdge = self.incoming[node][0]
        currentNode = self.edges[currentEdge]["src"]
        
        while not foundAnchor:
            if currentNode[0] == "A" or currentNode == "SOURCE":
                foundAnchor = True
            else:
                currentEdge = self.incoming[currentNode][0]
                currentNode = self.edges[currentEdge]["src"]
        return currentNode
    # 
    def tracePath(self, strain, start, end):
        """Return a [node, edge, node...] path tracing a strain between two anchors
        def tracePath(self:<GraphicalGenome>, strain:<str>, start:<str>, end:<str>):
        
        Parameters:
            strain: <str> - Founder or CC path that you want to follow. 
            start: <str> - Start anchor name
            end: <str> - End anchor name. Must occur later in the genome than the start anchor.
        """
        startAnchor = start
        endAnchor = end
        path = []
        current = startAnchor
        while current != endAnchor:
            for elem in self.outgoing[current]:
                if strain in self.edges[elem]["strain"]:
                    path.append(current)
                    path.append(elem)
                    current = self.edges[elem]["dst"]
        path.append(endAnchor)
        return path

    def entityPath(self, strain="B"):
        genome = self
        path = []
        for i in genome.outgoing["SOURCE"]:
            edge = genome.edges[i]
            # Get the first B edge
            if strain in edge["strain"]:
                current = i 
                path.append(current)
        while current != "SINK":
            # Floating node
            if current[0] == "F" or current[0] == "B":
                for i in genome.outgoing[current]:
                    edge = genome.edges[i]
                    if strain in edge["strain"]:
                        current = i
                        path.append(current)
            # Anchor node
            elif current[0] == "A":
                for edge in genome.outgoing[current]:
                    cEdge = genome.edges[edge]
                    if strain in cEdge["strain"]:
                        current = edge
                        path.append(current)
            # Edge
            else:
                current = genome.edges[current]["dst"]    
                path.append(current)
        # Add the sink node into the path
        for i in genome.incoming["SINK"]:
            edge = genome.edges[i]
            if strain in edge["strain"]:
                path.append(i)
                return path
        return path

    def findEdge(self, node, strain, ori = 'src'):#given a node and a strain, return the incoming or outgoing edge that contains the strain 
        het_a = strain + 'a'
        het_b = strain + 'b'
        edgelist = []
        if ori == 'src':
            for edge in self.outgoing[node]:
                strainlist = self.edges[edge]['strain']
                if strain in strainlist or het_a in strainlist or het_b in strainlist or len(self.outgoing[node])==1:
                    edgelist.append(edge)
            return edgelist
        elif ori == 'dst':
            for edge in self.incoming[node]:
                strainlist = self.edges[edge]['strain']
                if strain in strainlist or het_a in strainlist or het_b in strainlist or len(self.incoming[node])==1:
                    edgelist.append(edge)
            return edgelist
        
    def anchorPosList(self):   #populates global field: self.sortednodes but also returns it
        #will return a dictionary of anchors as keys with 8 [founder] values that show their coordinates
        tmp = []
        for i in sorted(self.nodes.keys(), key = lambda n: int(self.nodes[n]["B"])):
            tup = (i, int(self.nodes[i]["A"]), int(self.nodes[i]["B"]), int(self.nodes[i]["C"]),
              int(self.nodes[i]["D"]), int(self.nodes[i]["E"]), int(self.nodes[i]["F"]), 
              int(self.nodes[i]["G"]), int(self.nodes[i]["H"]))
            tmp.append(tup)
        dt = [("node", "U20"), ("A", "i8"), ("B", "i8"), ("C", "i8"), ("D", "i8"), 
          ("E", "i8"), ("F", "i8"), ("G", "i8"), ("H", "i8")]
        Sortednodes = numpy.array(tmp, dtype=dt)
        return Sortednodes    
            
    def boundingAnchors(self, strain, startPos, endPos): #
        """returns the bouding anchors of the given strain and coordinate positions, takes time for the first query
        def boundingAnchors(self:<GraphicalGenome>, strain:<str>, startPos:<int>, endPos:<int>):
        
        Parameters:
            strain: <str> - Founder or CC path that you want to follow. 
            start: <int> - Start linear coordinate
            end: <int> - End linear coordinate. Could be the same or larger than start.
        """
        assert endPos >= startPos, "endPos must be larger than the startPos"
        
        if len(self.sortednodes) == 0:
            self.sortednodes = self.anchorPosList()           
        Sortednodes = self.sortednodes
        i1 = numpy.searchsorted(Sortednodes[strain], startPos) - 1
        i2 = numpy.searchsorted(Sortednodes[strain], endPos) 
        if i1 < 0 and i2 != len(Sortednodes):
            start = 'SOURCE'
            end = Sortednodes["node"][i2]
        elif i1 < 0 and i2 == len(Sortednodes):
            start = 'SOURCE'
            end = 'SINK'
        elif i2 == len(Sortednodes):
            start = Sortednodes["node"][i1]
            end = 'SINK'
        else:
            start = Sortednodes["node"][i1]
            end = Sortednodes["node"][i2]
        return start, end
    
    def findNumOfFloatingNodes(self,startingAnchor, destAnchor):
        """returns the number of FloatingNodes between anchor pairs
        def findNumOfFloatingNodes(self (Graph),startingAnchor, destAnchor)
        
        Parameters:
            starting Anchor: <int> - Start anchor
            ending Anchor: <int> - End anchor
        """
        founders = ['A','B','C','D','E','F','G','H']
        dictionary = {}
        for founder in founders:
            #get the founder's path including floating nodes
            path = self.tracePath(founder,startingAnchor,destAnchor)
            #find all the floating nodes in the path and add them to a dictionary.
            for item in path:
                if 'B' in item:
                    try:
                        dictionary[item]+=1
                    except KeyError:
                        dictionary[item]=1
        return len(dictionary) 
    

    def Enumerate_anchor_by_anchor(self, start, end, sofar = [], strain = set('ABCDEFGH')):
        """helper function for getSubPath function, implement Path entity to self.subpath
        def Enumerate_anchor_by_anchor(self, start, end, sofar = [], strain = set('ABCDEFGH')):
        
        Parameters:
            start: <int> - Start anchor
            end: <int> - the anchor next to the start anchor
        """ 
        # termination by node
        if start == end:   
            sofar1 = sofar
            if len(strain)>0:
                self.subpath.append((sofar1, strain))
            return 

        if len(strain) <1:
            return  

        for edge in self.outgoing[start]:               
            self.Enumerate_anchor_by_anchor(self.edges[edge]['dst'], end, sofar = sofar + [start, edge], strain = strain & set(self.edges[edge]['strain']))

    
    def findsequence(self, pathlist, countinganchor = False):
        """Given a list for path entity, construct the sequence on the path
        def findsequence(self, pathlist, countinganchor = False):
        
        Parameters:
            pathlist: <list> - entity list that construct the path, including nodes and edges
            countinganchor: <logic> default False: do not counting the first and the last anchor, else, counting anchor sequence
        """
        seq = ''
        for item in pathlist:
            if item.startswith('A'):
                if countinganchor == True:
                    seq += self.nodes[item]['seq']
                else:
                    seq += '' # do not count anchor length
            elif item.startswith('L') or item.startswith('E')or item.startswith('K'):
                seq += self.edges[item]['seq']
            elif item.startswith('S') and item != "SOURCE" and item != 'SINK':
                seq += self.edges[item]['seq']
            else:
                seq += ''
        return seq

# Updated Aug 20, 2019, construct subsequences between each anchor pair then merge them together
    def getSubPath(self, startanchor, endanchor, init_strain = set("ABCDEFGH"), seq = False, counting_anchor = True): 
        """Get all subpath between two anchor, return a list with tuple: if seq = "False" return (itemlist, shared strain) , else return (sequence, shared strain), the path includes the start and end anchors
        def getSubPath(self, startanchor, endanchor, seq = False

        Parameters:
            startanchor: <str> - starting anchor
            endanchor: <str> ending anchor
            init_strain: <set> the set for all possible strains; default: the 8 founder strain set
            seq: <logic> default False: return (path entity, strain);if True: return (sequence, strain)
            counting_anchor <logic> default True, including the bounding anchor sequences
        """
        Subgraph = []
        anchor = startanchor
        while anchor != endanchor:
            self.subpath = []
            self.Enumerate_anchor_by_anchor(anchor, self.nextAnchor(anchor), strain = init_strain)

            if len(Subgraph) == 0:
                Subgraph = self.subpath
            else:
                New_Subgraph = []
                for prev_edgelist, prev_strain in Subgraph:
                    for new_edge, new_strain in self.subpath:                    
                        strain = prev_strain & new_strain
                        if len(strain)>0:
                            edgelist = prev_edgelist + new_edge
                            New_Subgraph.append((edgelist, strain))
                Subgraph = New_Subgraph        
            anchor = self.nextAnchor(anchor)

        # add the last anchor
        for item in Subgraph:
            item[0].append(endanchor)
        
        # return seq or edgelist
        if seq == True:
            Subsequence = []
            for itemlist, strain in Subgraph:
                Sequence = self.findsequence(itemlist, countinganchor = counting_anchor) # revise Sep 30, 19. When returning sequence, default include the start and end anchor sequences.
                Subsequence.append((Sequence, strain))
            return Subsequence
        else:
            return Subgraph

##
    # demote anchors as suggest in FigureS2
    def DemoteAnchor(self, anchor):
        """Demote Anchor sequence to floating nodes, edge, floating nodes
        def DemoteAnchor(self, anchor)
        
        Parameters:
            anchor: <str> - anchor name
        """
        chromo = int(anchor.split('.')[0][1:])
        #edgelist = Graph.edges.keys()
        # incoming edge
        self.incoming['F%02d.%08d' % (chromo, int(anchor.split('.')[1]))] = []
        for edge in self.incoming[anchor]:
            self.edges[edge]['dst'] = 'F%02d.%08d' % (chromo, int(anchor.split('.')[1]))
            self.incoming['F%02d.%08d' % (chromo, int(anchor.split('.')[1]))].append(edge)
        #outgoing edge
        index = 0
        self.outgoing['F%02d.%08d' % (chromo, int(anchor.split('.')[1])+1)] = []
        for edge in self.outgoing[anchor]:
            self.edges[edge]['src'] = 'F%02d.%08d' % (chromo, int(anchor.split('.')[1])+1)
            index += 1
            if edge.startswith('E'):
                name = 'E%02d.%07d%01d' % (chromo, int(anchor.split('.')[1])+1, index)
                #assert name not in edgelist
                self.edges[name] = self.edges[edge]
                self.outgoing['F%02d.%08d' % (chromo, int(anchor.split('.')[1])+1)].append(name)
                dst = self.edges[edge]['dst']
                self.incoming[dst].remove(edge)
                self.incoming[dst].append(name)
                self.edges.removeEdge(edge)
            else:
                self.outgoing['F%02d.%08d' % (chromo, int(anchor.split('.')[1])+1)].append(edge)
        # create new edge
        newname = 'E%02d.%07d%01d' % (chromo, int(anchor.split('.')[1]), 1)
        New_Dict = {}
        New_Dict['seq'] = self.nodes[anchor]['seq']
        New_Dict['strain'] = list('ABCDEFGH')
        New_Dict['src'] = 'F%02d.%08d' % (chromo, int(anchor.split('.')[1]))
        New_Dict['dst'] = 'F%02d.%08d' % (chromo, int(anchor.split('.')[1])+1)
        New_Dict['singleton'] = 'Conserved'
        New_Dict['variants'] = '45='
        for key in ['gene', 'exon', 'repeatclass']:
            if key in self.nodes[anchor].keys():
                New_Dict[key] = self.nodes[anchor][key]

        self.edges[newname] = New_Dict
        self.outgoing['F%02d.%08d' % (chromo, int(anchor.split('.')[1]))]= [newname]
        self.incoming['F%02d.%08d' % (chromo, int(anchor.split('.')[1])+1)]= [newname]
        # demote anchor
        self.nodes.removeNode(anchor)
        
## Graph-coordinate System
## linear-graph entity transformation
## Vertical Projection
## Homologous Sequence abstraction
## Under Construction
    def processCigar(self, cigar):
        """Helper Function, may not be used directly, expand Cigar string
        
        Parameters:
            cigar: <str> - compressed cigar
        """
        out = ''
        N = 0
        for symbol in cigar:
            if symbol in '0123456789':
                N = 10*N + int(symbol)
            else:
                #if (symbol != 'D'):
                if (N == 0):
                    out += symbol
                else:
                    out += N*symbol
                N = 0
        return out

    def combineCigar(self, cigar):
        """Helper Function, may not be used directly, compress Cigar string
        
        Parameters:
            cigar: <str> - expanded cigar
        """
        cigar = cigar +'$'
        out = ''
        N = 0
        start = 0
        for i in range(1,len(cigar)):
            if cigar[i-1] != cigar[i]:
                out += str(i-start) + cigar[i-1]
                start = i
        return out  
    
## helper function for get floating node coordinate: start from floating nodes, traverse forward until the first anchor
    def Traverseforward(self, fnode, sofar = '', strain = set('ABCDEFGH')):
        """Helper Function, may not be used directly, traverse path from floating node to find the proximal anchor
        subpath are stored in self.subpath
        Parameters:
            fnode: <str> - floating node name
        """
        # termination by anchor node
        if fnode.startswith('A'):   
            sofar1 = self.nodes[fnode]['seq'] + sofar
            if len(strain)>0:
                self.subpath.append((sofar1, strain))
            return 

        if len(strain) <1:
            return  

        for edge in self.incoming[fnode]:               
            self.Traverseforward(self.edges[edge]['src'], sofar =  self.edges[edge]['seq'] + sofar , strain = strain & set(self.edges[edge]['strain']))

# given node name, return dictionary for node coordinates: for floating node, return the coordinates for the base right before the floating node (1 base coordinates)
    def Nodecoordinates(self, node):
        """Given node name, return dictionary for linear coordinates for each node
        for floating node, by coordinates it refer to the base right after the floating node (1 base coordinates)
        Parameters:
            fnode: <str> - node name, compatible for both anchor and floating nodes
        """
        node_coord = {}
        if node.startswith('A'):
            for s in 'ABCDEFGH':
                node_coord[s] = int(self.nodes[node][s])
        elif node == 'SOURCE':
            for s in 'ABCDEFGH':
                node_coord[s] = 0        
        else:
            self.subpath = []
            self.Traverseforward(node)
            prev_anchor = self.prevAnchor(node)
            for seq, strain in self.subpath:
                for s in strain:
                    node_coord[s] = int(self.nodes[prev_anchor][s]) + len(seq) 
        return node_coord
    
    ## given graph entity(anchor or edge), offset and strain, return linear coordinates
    def entity_to_linearcoord(self, entity, offset, strain):
        """Given Graph entity (anchor or edge), offset and strain attributes, return linear coordinates for the base position
        
        Parameters:
            entity: <str> - anchor or edge that contain sequences
            offset: <int> - base position on the entity
            strain: <str> - specify strain attributes for the linear coords, given strain attributes to the entity, this could be multiple choices.
        """
        if entity.startswith('A'):
            node_coord = self.Nodecoordinates(entity)
            return node_coord[strain] + offset
        else:
            node = self.edges[entity]['src']
            node_coord = self.Nodecoordinates(node)
            if node.startswith('A'):
                return node_coord[strain] + 45 + offset
            else:
                return node_coord[strain] + offset
            
    ## given linear coordinates in strain, return graph entity(anchor or edge) plus offset 
    def linear_to_entityoffset(self, coord, strain):
        """Given linear coordinates, return graph entity + offset
        
        Parameters:
            coord: <int> - linear coordinate 
            strain: <str> - strain attributes for the linear coords 
        """
        sanchor, eanchor = self.boundingAnchors(strain, coord, coord)
        offset = coord - int(self.nodes[sanchor][strain])
        itemlist = self.tracePath(strain, sanchor, eanchor,)
        seq = 0
        for item in itemlist:
            if item.startswith('A'):
                newseq = seq + len(self.nodes[item]['seq'])
            elif item.startswith('B') or item.startswith('F') or item == "SOURCE" or item == 'SINK':
                newseq = seq
            else:
                newseq = seq + len(self.edges[item]['seq'])
            if newseq > offset:
                break
            else:
                seq = newseq
        current_offset = offset - seq
        entity = item
        return (entity, current_offset)
    
    # Graph Entity plus offset to anchors, return dictionary
    def entity_to_anchor(self, entity, offset, ref_node = False):
        """Given Graph entity (anchor or edge), offset and strain attributes, return linear coordinates for the base position

        Parameters:
            entity: <str> - anchor or edge that contain sequences
            offset: <int> - base position on the entity
            strain: <str> - specify strain attributes for the linear coords, given strain attributes to the entity, this could be multiple choices.
            ref_node: <str> - default proximal anchor, or 
                            - "SOURCE" return linear coordinates, specifying the anchor that refer to
        """
        Strain_Offset = {}
        if entity.startswith('A'):
            if ref_node == False: # default proximal anchor
                for s in 'ABCDEFGH':
                    Strain_Offset[s] = offset
                return (entity, Strain_Offset)
            elif ref_node == "SOURCE": # linear coordinates
                node_coord = self.Nodecoordinates(entity)
                for s in "ABCDEFGH":
                    Strain_Offset[s] = node_coord[s] + offset
                return ("SOURCE", Strain_Offset) 
            else: # other anchor
                node_coord = self.Nodecoordinates(entity)
                for s in "ABCDEFGH":
                    ref_coord = int(self.nodes[ref_node][s])
                    Strain_Offset[s] = node_coord[s] - ref_coord + offset         
                return (ref_node, Strain_Offset)
        else:
            node = self.edges[entity]['src']
            node_coord = self.Nodecoordinates(node)
            strainlist = node_coord.keys()
            if ref_node == False: # default
                if node.startswith('A'):
                    for s in strainlist:
                        Strain_Offset[s] = node_coord[s] + 45 + offset
                    return (node, Strain_Offset)
                else:
                    proximal_anchor = self.prevAnchor(node)
                    for s in strainlist:
                        pos = int(self.nodes[proximal_anchor][s])
                        Strain_Offset[s] = node_coord[s] + offset - pos                
                    return (proximal_anchor, Strain_Offset)

            elif ref_node == 'SOURCE': # Linear            
                if node.startswith('A'):
                    for s in strainlist:
                        Strain_Offset[s] = node_coord[s] + 45 + offset
                else:
                    for s in strainlist:
                        Strain_Offset[s] = node_coord[s] + offset

                return (ref_node, Strain_Offset)
            else:
                if node.startswith('A'):
                    for s in strainlist:
                        ref_pos = int(self.nodes[ref_node][s])
                        Strain_Offset[s] = node_coord[s] + 45 + offset - ref_pos
                else:
                    for s in strainlist:
                        ref_pos = int(self.nodes[ref_node][s])
                        Strain_Offset[s] = node_coord[s] + offset - ref_pos
                return (ref_node, StrainOffset)

    def linear_to_anchoroffset(self, coord, strain, ref_node = False ):
        """Given linear coordinates, return graph entity + offset
        Parameters:
            coord: <int> - linear coordinate 
            strain: <str> - strain attributes for the linear coords 
            ref_node: <str> - default proximal anchor, or 
                            - "SOURCE" return linear coordinates, specifying the anchor that refer to
        """
        sanchor, eanchor = self.boundingAnchors(strain, coord, coord)
        offset = coord - int(self.nodes[sanchor][strain])
        if ref_node == False:        
            return (sanchor, offset)
        elif ref_node == "SOURCE":
            return (ref_node, coord)
        else:
            current_offset = coord - int(self.nodes[ref_node][strain])
            return(ref_node, current_offset)
        
        
    def Parallel_edge_offset_exchange(self, entity, offset, alter_strain, strain = 'B'):
        """Alignment-based entity offset exchange
        (only apply to parallel edge, edge share the same src and dst, offset exchange in Path scale are not finished yet)
        Given entity+offset on B path, return the position on the alternative strain
        If alignment results are not applicable, return None
        
        Parameters:
            entity: <str> - Graph entity on B path
            offset: <int> - offset on the entity
            alter_strain: <str> - strain attributes for the target alternative path position
            strain: <default> "B" path, when multi-alignment cigar are added, this could be further implemented
        """
        if entity.startswith('A'):
            alter_item = entity
            return (alter_item, offset)
        else:
            s_anchor = self.edges[entity]['src']        
            for edge in self.outgoing[s_anchor]:
                if alter_strain in self.edges[edge]['strain']:
                    alter_item = edge
                    assert self.edges[entity]['dst'] == self.edges[alter_item]['dst'], "Not parallel edge"
                    break
            else:
                return 

            if 'variants' in self.edges[alter_item].keys():
                cigar = self.processCigar(self.edges[alter_item]['variants'])
                ref_index = 0
                alt_index = 0
                if ref_index == offset:
                    if cigar[ref_index] == 'D':
                        return (alter_item, alt_index, 'D')
                    else:
                        return (alter_item, alt_index)

                for i in range(len(cigar)): # index of all cigar string
                    if cigar[i] in 'MDS=X':
                        ref_index += 1 # coordinates in ref
                    if cigar[i] in 'MIS=X':
                        alt_index += 1 # coordinate in alt
                    if ref_index == offset:
                        if cigar[i] == 'D':
                            return (alter_item, alt_index, 'D')
                        else:
                            return (alter_item, alt_index)
            else:
                return 
            
    def Vertical_linearpos_exchange(self, coord, alter_strain, strain = 'B'):
        """Alignment-based entity offset exchange
        Given linear coord on B path, return corresponding linear coordinates on the alternative genome
        If alignment results are not applicable, return None
        
        Parameters:
            coord: <int> - linear coordinates on B path
            alter_strain: <str> - strain attributes for the target alternative genome position
            strain: <default> "B" path, when multi-alignment cigar are added, this could be further implemented
        """
        entity, offset = self.linear_to_entityoffset(coord, 'B')
        if self.Parallel_edge_offset_exchange(entity, offset, alter_strain) != None:
            Item = self.Parallel_edge_offset_exchange(entity, offset, alter_strain)
            alt_entity = Item[0]
            alt_offset = Item[1]
            return self.entity_to_linearcoord(alt_entity, alt_offset, alter_strain)
        else:
            return