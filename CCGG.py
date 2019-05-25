import cjson
import numpy
import json
import re
import time

class Nodes(dict):
    def __init__(self, inp):
        self.nodes = inp
    # Override getitem to account for lazy evaluation
    def __getitem__(self,key):
        return self.getNode(key)
    
    def __setitem__(self, key, value):
        self.nodes[key] = value
    
    def __contains__(self, key):
        return key in self.nodes
            
    def __repr__(self):
        return self.nodes
    
    def keys(self):
        return self.nodes.keys()
    
    # getNode allows for lazy evaluation
    def getNode(self, nodename):
        if "hdr" in self.nodes[nodename].keys():
            seq = self.nodes[nodename]["seq"]
            hdr = self.nodes[nodename]["hdr"]
            del self.nodes[nodename]
            self.nodes[nodename] = cjson.decode(hdr)
            self.nodes[nodename]["seq"] = seq
            return self.nodes[nodename]
        return self.nodes[nodename]
    
    def returnDict(self):
        return self.nodes
    
    def removeNode(self, nodename):
        output = self.nodes.pop(nodename, None)
        if output == None:
            output = "Key does not exist in dictionary"
        return output
    
class Edges(dict):
    def __init__(self, inp):
        self.edges = inp
        
    def __getitem__(self,key):
        return self.getEdge(key)
    
    def __setitem__(self, key, value):
        self.edges[key] = value
    
    def __contains__(self, key):
        return key in self.edges
    
    def __repr__(self):
        return self.edges
            
    def keys(self):
        return self.edges.keys()            
    
    def getEdge(self, edgename):
        if "hdr" in self.edges[edgename].keys():
            seq = self.edges[edgename]["seq"]
            hdr = self.edges[edgename]["hdr"]
            del self.edges[edgename]
            self.edges[edgename] = cjson.decode(hdr)
            self.edges[edgename]["seq"] = seq
            return self.edges[edgename]
        return self.edges[edgename]    
       
    def returnDict(self):
        return self.edges
    
    def removeEdge(self, edgename):
        output = self.edges.pop(edgename, None)
        if output == None:
            output = "Key does not exist in dictionary"
        return output

class GraphicalGenome:
    
    def __init__(self, nodefile, edgefile, verbose=False):
        self.verbose = verbose
        self.nodes, self.edges, self.outgoing, self.incoming, self.max_node_name, self.max_edge_name, self.genes = self.loadChromosomeGraph(nodefile, edgefile)
        self.nodes = Nodes(self.nodes)
        self.edges = Edges(self.edges)
        
    def __repr__(self):
        return "Number of Nodes %d - Number of Edges %d" % (len(self.nodes), len(self.edges))
        
    def __str__(self):
        return "Origin - Nodefile %s    Edgefile %s" % (nodefile, edgefile)
    
    def encodeASCII(self, data):
        return dict((k.encode('ascii'), [i.encode('ascii') for i in v] if isinstance(v, list) else v.encode('ascii')) for k, v in data.items())

    

# Inputs:
    # nodefile - Absolute path to the nodefile you wish to open
    # edgefile - Absolute path to the edgefile you wish to open
# Outputs the nodes and edges as dictionaries as well as auxilary dictionaries containing node -> list(edge) pairs 
    def loadChromosomeGraph(self, nodefile, edgefile):
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
                node =  hdr[0:12]
                nodehdr = hdr[13:]
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
            hdr = fp.readline()[1:]
            seq = fp.readline()
            while hdr and seq:
                # Extract the src, dst, key from the header
                src = hdr[21:33]
                dst = hdr[42:54]
                key = hdr[0:12]
                
                # Source and sink nodes are not the same size as other nodes so special case needed
                if "SOURCE" in src:
                    dst = hdr[36:48]
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
                edges[key]["hdr"] = hdr[13:]
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
    
# Given a CC path, look for all node pairs that do not contain an edge between them that has the specified CC strain on them. 
# * Parameters: 
#     * (str) strain - Strain you wish to trace 
# * (str) chromo - Chromosome you are tracing the strain through 
# * Output: 
#     * list(tuple(node, node)) - List of node tuples that do not have any edges between them with the specified CC strain on them.
# THIS CURRENTLY WORKS FOR EDGES WITHOUT FLOATING NODES, WHEN WE REINTRODUCE FLOATING NODES THIS NEEDS TO BE FIXED
    def findMissingPaths(self, strain, chromo):
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
        reconstructed = self.reconstructSequence(strain, path)
        with open(filedest, "w+") as seqfile:
            seqfile.write(reconstructed)
            
    # Given an anchor, provide the next anchor
    def nextAnchor(self, node):
        while True:
            edge = self.outgoing[node][0]
            if self.edges[edge]["dst"][0] == "A" or self.edges[edge]["dst"] == "SINK":
                return self.edges[edge]["dst"]
            else:
                node = self.edges[edge]["dst"]
                
    # Return a [node, edge, node...] path tracing a strain between two anchors
    def tracePath(self, strain, start, end):
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
    
    
    # Bounding Anchors
        # Given 