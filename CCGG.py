import cjson
import numpy
import json
import time

# Overridden dict class to allow for lazy evaluation of nodes
class Nodes(dict):
    def __init__(self, inp):
        self.nodes = inp
        
    # Override getitem to account for lazy evaluation
    def __getitem__(self,key):
        return self.getNode(key)

    def getNode(self, nodename):
        if "hdr" in self.nodes[nodename].keys():
            seq = self.nodes[nodename]["seq"]
            hdr = self.nodes[nodename]["hdr"]
            del self.nodes[nodename]
            self.nodes[nodename] = cjson.decode(hdr)
            self.nodes[nodename]["seq"] = seq
            return self.nodes[nodename]
        return self.nodes[nodename]

# Overridden dict class to allow for lazy evaluation of edges
class Edges(dict):
    def __init__(self, inp):
        self.edges = inp
    
    # Override getitem to account for lazy evaluation
    def __getitem__(self,key):
        return self.getEdge(key)
    
    def getEdge(self, edgename):
        if "hdr" in self.edges[edgename].keys():
            seq = self.edges[edgename]["seq"]
            hdr = self.edges[edgename]["hdr"]
            del self.edges[edgename]
            self.edges[edgename] = cjson.decode(hdr)
            self.edges[edgename]["seq"] = seq
            return self.edges[edgename]
        return self.edges[edgename]    
    
class GraphicalGenome:
    def __init__(self, nodefile, edgefile):
        self.nodes, self.edges, self.outgoing, self.incoming, self.max_node_name, self.max_edge_name, self.genes = self.loadChromosomeGraph(nodefile, edgefile)
        self.nodes = Nodes(self.nodes)
        self.edges = Edges(self.edges)
#         tmp = []
#         for i in sorted(self.nodes.keys(), key = lambda n: int(self.nodes[n]["B"])):
#             tup = (i, int(self.nodes[i]["A"]), int(self.nodes[i]["B"]), int(self.nodes[i]["C"]),
#                   int(self.nodes[i]["D"]), int(self.nodes[i]["E"]), int(self.nodes[i]["F"]), 
#                   int(self.nodes[i]["G"]), int(self.nodes[i]["H"]))
#             tmp.append(tup)
#         dt = [("node", "U20"), ("A", "i8"), ("B", "i8"), ("C", "i8"), ("D", "i8"), 
#               ("E", "i8"), ("F", "i8"), ("G", "i8"), ("H", "i8")]
#         self.sortednodes = numpy.array(tmp, dtype=dt)

    def __repr__(self):
        return "Number of Nodes %d - Number of Edges %d" % (len(self.nodes), len(self.edges))
        
    def __str__(self):
        return "Origin - Nodefile %s    Edgefile %s" % (nodefile, edgefile)
    
    def encodeASCII(self, data):
        return dict((k.encode('ascii'), [i.encode('ascii') for i in v] if isinstance(v, list) else v.encode('ascii')) for k, v in data.items())

    def loadChromosomeGraph(self, nodefile, edgefile):
        """
        Inputs:
            nodefile - Absolute path to the nodefile you wish to open
            edgefile - Absolute path to the edgefile you wish to open
        Outputs the nodes and edges as dictionaries as well as auxilary dictionaries containing node -> list(edge) 
        pairs 
        """
        # Instantiate returned datastructures
        max_node = 0
        max_edge = 0
        nodes = {}
        edges = {}
        sources = {}
        destinations = {}
        genes = {}

        # Load Nodefile
        snh = time.time()
        nodehdr, nodeseq = self.loadFasta(nodefile)
        print "Node FASTA load %5.2f" % (time.time() - snh)
        # Populate base node dictionary without evaluating any JSON
        sn = time.time()
        for i in xrange(len(nodehdr)):
            part = nodehdr[i].split(';',1)
            node = part.pop(0)
            nodes[node] = {}
            nodes[node]["hdr"] = part.pop(0)
            nodes[node]["seq"] = nodeseq[i]
            
#             # Check existence of gene annotation and add to dictionary
#             gene_full = re.search(gene_pattern, nodehdr[i])
#             if gene_full:
#                 gene = gene_full.group()
#                 genelist = gene[8:-1].split(",")
#                 for i in genelist:
#                     genekey = i[1:-1]
#                     genes[genekey] = genes.get(genekey, []) + [node]
                
        print "Node loop %5.2f" % (time.time() - sn)
#         del nodehdr
#         del nodeseq
        # Load Edgefile
        
        seh = time.time()
        edgehdr, edgeseq = self.loadFasta(edgefile)
        print "Edge FASTA load %5.2f" % (time.time()-seh)
        se = time.time()
        for i in xrange(len(edgehdr)):
            part = edgehdr[i].split(";",1)
#             src = part[1].split('"src":',1)[1].split('"', 2)
#             dst = src[2].split('"dst":', 1)[1].split('"',2)[1]
#             src = src[1]
#             key = part[0]
            src = edgehdr[i][21:33]
            dst = edgehdr[i][42:54]
            key = edgehdr[i][0:12]

            # Populate base edge dictionary without evaluating any JSON
            edges[key] = {}
            edges[key]["hdr"] = edgehdr[i][13:]
            edges[key]["seq"] = edgeseq[i]

            # Populate sources and destinations dictionaries
            sources[src] = sources.get(src, []) + [key]
            destinations[dst] = destinations.get(dst, []) + [key]
            
            # Check existence of gene annotation and add to dictionary
#             gene_full = re.search(gene_pattern, edgehdr[i])
#             if gene_full:
#                 gene = gene_full.group()
#                 genelist = gene[8:-1].split(",")
#                 for i in genelist:
#                     genekey = i[1:-1]
#                     genes[genekey] = genes.get(genekey, []) + [key]
                
#             # Update edge global counters
            if "F" in key and "S" not in key and "K" not in key:
                max_edge = max(max_edge, int(key[-7:]))
#             # Update node global counters
            if "F" == dst[0]:
                max_node = max(max_node, int(dst[-8:]))
            if "F" == src[0]:
                max_node = max(max_node, int(src[-8:]))
        print "Edge loop %5.2f" % (time.time() - se)
        del edgehdr
        del edgeseq
        return nodes, edges, sources, destinations, max_node, max_edge, genes

    
    ## FIXZ THIS TO NOT EVALUATE JSON THAT HAS NOT BEEN EVALUATED YET

    def writeFasta(input_dict, filename, keylist=["src", "dst"]):
        """ Writes FASTA formatted file in the correct order to allow for lazy evaluation in init
        """
        sorted_keys = sorted(input_dict.keys()) 
        with open(filename, "w+") as fastafile:
            for edge in sorted_keys:
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

    def loadFasta(self, filename):
        """ Parses a classically formatted and possibly 
            compressed FASTA file into a list of headers 
            and fragment sequences for each sequence contained.
            The resulting sequences are 0-indexed! """
        if (filename.endswith(".gz")):
            fp = gzip.open(filename, 'rb')
        else:
            fp = open(filename, 'rb')
        # split at headers
        data = fp.read().split('>')
        fp.close()
        # ignore whatever appears before the 1st header
        data.pop(0)     
        headers = []
        sequences = []
        for sequence in data:
            lines = sequence.split('\n')
            headers.append(lines.pop(0))
            sequences.append(''.join(lines))
        return (headers, sequences)

    def showParallelEdges(self, edge):
        """  
        """
        src = self.edges[edge]["src"]
        dst = self.edges[edge]["dst"]
        parallelEdges = []
        found = True
        outgoing = set(self.outgoing[src])
        incoming = set(self.incoming[dst])
        return list(outgoing and incoming)
    
        for i in self.outgoing[src]:
            currentEdge = i
            while(self.edges[currentEdge]["dst"] != dst):
                currentEdge = dst
                if self.edges[currendEdge] == "SINK":
                    found = False
                    break
            if found:
                parallelEdges.append(i)
        return parallelEdges

    def deleteNode(self, nodename):
        inEdges = self.incoming[nodename]
        prevAnchor = ""
        for i in inEdges:
            if "A" in self.edges[i]["src"]:
                prevAnchor = self.edges[i]["src"]
                break
        
        outEdges = self.outgoing[nodename]
        nextAnchor = ""
        for j in outEdges:
            if "A" in self.edges[j]["dst"]:
                nextAnchor = self.edges[j]["dst"]
                break
        
        for inedge in inEdges:
            self.edges[inedge]["dst"] = nextAnchor
        for outedge in outEdges:
            self.edges[outedge]["src"] = prevAnchor
        del self.nodes[nodename]
        

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

    # Bounding anchors reqiures strains to be in founder sequences (ABCDEFGH)
    def boundingAnchors(self, strain, startPos, endPos):
        i1 = numpy.searchsorted(self.sortednodes[strain], startPos) - 1
        i2 = numpy.searchsorted(self.sortednodes[strain], endPos) - 1
        if i1 < 0:
            i1 = 0
        if i2 < 0:
            i2 = 0
        start = self.sortednodes["node"][i1]
        end = self.sortednodes["node"][i2]
        return start, end

    def writeSequence(self, strain, filedest, path=0):
        reconstructed = self.reconstructSequence(strain, path)
        with open(filedest, "w+") as seqfile:
            seqfile.write(reconstructed)
            
        