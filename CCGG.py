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
            het = strain
        elif path == 1:
            het += "a"
        else:
            het += "b"

        for src in outgoing["SOURCE"]:
            for edge in edges[src]["strain"]:
                if het in edge or strain in edge:
                    source = src
        if source == "":
            print "strain not found on any source path"
            return ""
        conseq = ""

        # Iterate through the edges in order based on the destinations and outgoing edges from those destinations
        currentEdge = source
        firstNode = edges[source]["dst"]
        numberN = int(nodes[firstNode][source[-1]]) - 1 - len(edges[source]["seq"])
        conseq += "N" * numberN
        while True:
            dst = edges[currentEdge]["dst"]
            if dst == "SINK":
                conseq += edges[currentEdge]["seq"]
                conseq += ("N" * (int(edges[currentEdge]["seqlen"]) - len(conseq)))
                return conseq.upper()
            else:
                if "F" in dst:
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
        while True:
            edge = self.outgoing[node][0]
            if self.edges[edge]["dst"][0] == "A" or self.edges[edge]["dst"] == "SINK":
                return self.edges[edge]["dst"]
            else:
                node = self.edges[edge]["dst"]
                
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
    