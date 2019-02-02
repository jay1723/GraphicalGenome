import cjson
import numpy
import time
import json
from collections import defaultdict
class GraphicalGenome:
    
    def __init__(self, nodefile, edgefile):
        self.nodes, self.edges, self.outgoing, self.incoming = self.loadChromosomeGraph(nodefile, edgefile)
        tmp = []
        for i in sorted(self.nodes.keys(), key = lambda n: int(self.nodes[n]["B"])):
            tup = (i, int(self.nodes[i]["A"]), int(self.nodes[i]["B"]), int(self.nodes[i]["C"]),
                  int(self.nodes[i]["D"]), int(self.nodes[i]["E"]), int(self.nodes[i]["F"]), 
                  int(self.nodes[i]["G"]), int(self.nodes[i]["H"]))
            tmp.append(tup)
        dt = [("node", "U20"), ("A", "i8"), ("B", "i8"), ("C", "i8"), ("D", "i8"), 
              ("E", "i8"), ("F", "i8"), ("G", "i8"), ("H", "i8")]
        self.sortednodes = numpy.array(tmp, dtype=dt)

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
        node = {}
        edge = {}
        sources = {}
        destinations = {}
        hdr, seq = self.loadFasta(nodefile)
        # Nodes
        for i in xrange(len(hdr)):
            part = hdr[i].split(';')
            key = part.pop(0)

            node[key] = cjson.decode(part.pop(0)) # popping the second item always pops the JSON
            node[key]['seq'] = seq[i]
        hdr, seq = self.loadFasta(edgefile)
        # Edges
        for i in xrange(len(hdr)):
            part = hdr[i].split(';')
            key = part.pop(0)
            edge[key] = cjson.decode(part.pop(0)) # popping the second item always pops the JSON
            # Make these their own data structures
            sources[edge[key]['src']] = sources.get(edge[key]['src'], []) + [key]
            destinations[edge[key]['dst']] = destinations.get(edge[key]['dst'], []) + [key]
            edge[key]['seq'] = seq[i]
        return node, edge, sources, destinations
    

    def writeFasta(self, filename, file_dict):
        """
        Inputs:
            filename - Absolute path to where you want the file to be written including the file name
            file_dict - The updated node/edge dictionary that you wish to write
        Outputs FASTA formatted file at path specified in argument. Assumes that 
        """
        sorted_keys = sorted(file_dict.keys())
        with open(filename, "w+") as fastafile: 
            for i in sorted_keys:
                line = ">" + i + ";"
                obj = {}
                for j in file_dict[i].keys():
                    if j == 'seq':
                        continue
                    obj[j] = file_dict[i][j]
                line += json.dumps(obj, separators=(",", ":"))
                line += "\n" + file_dict[i]['seq'] + "\n"
                fastafile.write(line)
        return

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
            

                
            

