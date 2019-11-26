# Graphical Genome Class
Found on github [here](https://github.com/jay1723/GraphicalGenome)

## What is the Graphical Genome
The Graphical Genome is an implementation of a pangenome that uses a directed graph to represent sequence data from disparate animals in one representation. 

## Graphical Genome basic usage
The Graphical Genome, when initialized, returns an object with four instance data dictionaries, nodes, edges, incoming, outgoing. 

`GraphicalGenome.nodes` : Dictionary of nodes. Referenced with string node names
* Each Node currently contains the following keys:
	* `A-H` = (str) : Each key A-H contains the relative position the original founder sequence that the node 45-mer can be found
	* `repeatclass` = list(str) : Contains the list of all repeat classes for the given 45-mer
	* `seq` = (str) : Sequence
`GraphicalGenome.edges` : Dictionary of edges. Referenced with string edge names
* Each Edge in the dictionary currently contains the following keys:
	* `src` = (str) : Source node for the edge
	* `dst` = (str) : Destination node for the edge
	* `strains` = list(str) : List of strains that lie on the edge
	* `seq` = (str) : Sequence
`GraphicalGenome.outgoing` : Dictionary of incoming edges. Referenced with string node names
`GraphicalGenome.incoming` : Dictionary of outgoing edges. Referenced with string edge names


### Source and Sink Edges
The 8 founder strain sequences start and end with a "Source" and "Sink" edge. 

The following is an example of a source and sink edge on chromosome 1 for strain A. While normal edges begin with "E", source nodes begin with "S" and sink nodes begin with "K". The remainder of their name includes the relevant chromosome information and the strain for which they are a source of sink. 

```python
GraphicalGenome.edges["S01.00000001"] = {"src" : "SOURCE", "seq" : <sequence>, "dst" : <node>, "strain" : ["A"]}
```

```python
GraphicalGenome.edges["K01.00000001"] = {"src" : <node>, "seq" : <sequence>, "dst" : "SINK", "strain" : ["A"]}
```
## API

### Initialization
Initalize a new graphical genome instance
* Parameters:
	* (str) `nodefile` - Absolute path to the nodefile
    * (str) `edgefile` - Absolute path to the edgefile
* Optional Parameters
	* (bool) `verbose=False` - Prints verbose output as Graphical genome is being generated.
	* (int) `enamelen=12` - Defines the edge key name length.
	* (int) `nnamelen=12` - Defines the node key name length.
* Output: new instance of Graphical Genome

```python
import CCGG
genome = GraphicalGenome(nodefile, edgefile) # Initialize Graphical Genome instance with the specified node and edge file
genome.edges[<edgename>]["strains"] # Example of looking at the list of strains given an edgename
genome.nodes[<nodename>]["repeatclass"] # Example of returning the list of repeat classes of a given node
```

### `writeFasta`
Write updated node and edge dictionaries in a FASTA format.
* Parameters:
	* (str) `filename` - Absolute path for the file that you would like to write. Files that do not exist will be created in the specified location with the specified name. 
	* (dict) `file_dict` - Node or Edge dictionary that you would like to write. 
	* (list) `keylist` - List of strings of dictionary keys tha tyou want to ignore during the file write. This will require you to handle these keys on your own. 
* Output: FASTA file created in the specified location
```python
import CCGG
genome = GraphicalGenome(nodefile, edgefile) # Initialize Graphical Genome instance with the specified node and edge file
genome.nodes[<nodename>]["new_annotation"] = <new_data> # Add a new annotation to the specified node
genome.writeFasta("FASTAFolder/newfile.fa", newinstance.nodes) # Writes updated node dictionary to the specified FASTA file
```

### `returnNodeEdgePath`
Returns a list of [node,edge,node,edge...] that traces out a path for a strain on a given chromosome.
* Parameters:
	* (str) `strain` - Founder or CC strain you wish to trace.
	* (str) `chromo` - Chromosome you are iterating through.

```python
import CCGG
genome = GraphicalGenome(nodefile, edgefile) # Initialize Graphical Genome instance with the specified node and edge file
genome.returnNodeEdgePath(<strain>, <chromo>) # Get the list of nodes and edges for a strain
```

### `findMissingPaths`
Given a Collaborative Cross path, look for all node pairs that do not contain an edge between them that has the specified Collaborative Cross strain on them.
* Parameters:
	* (str) `strain` - Founder or CC strain you wish to trace.
	* (str) `chromo` - Chromosome you are iterating through. 
* Deprecation warning:
	*  THIS CURRENTLY WORKS FOR EDGES WITHOUT FLOATING NODES, WHEN FLOATING NODES ARE ADDED THIS NEEDS TO BE UPDATED
* Output:
	* (list(tuples(node, node))) - Retuns a list of node tuples that represent the missing paths
```python
import CCGG
genome = GraphicalGenome(nodefile, edgefile) # Initialize Graphical Genome instance with the specified node and edge file
paths = genome.findMissingPaths(<strain>, <chromo>) # Get the missing paths as a list
```

### `reconstructSequence` (May need to be updated)
Reconstruct the sequence given a strain (ABCDEFGH or CC) and can return more than one sequence if heterozygosity is known.
* Parameters:
	* (str) `strain` - The strain that you want to trace
* Optional parameters:
	* (int) `path` - 0 for base path, 1 for first heterozygous path identified, 2 for second heterozygous path identified. 
* Output:
	The entire genomic sequence, from start to finish, for the given strain
```python
import CCGG
genome = GraphicalGenome(nodefile, edgefile) # Initialize Graphical Genome instance with the specified node and edge file
sequence = genome.reconstructSequence(<strain>)
```


### `nextAnchor`
Return the next monotonically incrementing anchor given a nodename in the graph
* Parameters:
	* (str) `anchor` - The current anchor node
* Output:
	* (str) - The next monotonically increasing anchor node
```python
import CCGG
# Basic usage
genome = CCGG.GraphicalGenome(nodefile, edgefile)
nextAnchor = genome.nextAnchor(<anchor>)

# Iterating through the genome anchor to anchor
currentAnchor = <anchor>
while currentAnchor != "SINK":
	# Add any analysis here #
	currentAnchor = genome.nextAnchor(currentAnchor)

```
### `prevAnchor`
Return the previously monotonically decreasing anchor given a nodename in the graph
* Parameters:
	* (str) `anchor` - The current anchor node
* Output:
	* (str) - The next monotonically decreasing anchor node
```python
import CCGG
# Basic usage
genome = CCGG.GraphicalGenome(nodefile, edgefile)
nextAnchor = genome.prevAnchor(<anchor>)

# Iterating through the genome backwards anchor to anchor
currentAnchor = <anchor>
while currentAnchor != "SOURCE":
	# Add any analysis here #
	currentAnchor = genome.prevAnchor(currentAnchor)
```
### `tracePath`
Return a [node, edge, node...] path tracing a strain between two anchors
* Parameters
	* (str) `strain` - the strain that you want to trace through the genome
	* (str) `start` - the anchor where you would like to start the trace
	* (str) `end` - the anchor where you would like to end the trace
* Output
	* (list(node,edge,node,...)) - Output of alternating node, edge pairs that define a path through the graph from the start node to the end node. 

```python
import CCGG
genome = CCGG.GraphicalGenome(nodefile, edgefile) # Instantiate the genome object
path = genome.tracePath(<strain>, <start>, <end>)
```

### `entityPath`
Return the path through the entire genome for a given strain as a node,edge,node... list.
* Optional Parameters:
	* (str) `strain` - The strain you would like to trace through the genome. (ABCDEFGH or CC)
		* Defaults to "B"
* Output:
	* The [node, edge, node...] path from the beginning to the end of the genome that traces the specified strain. 
```python
import CCGG
genome = CCGG.GraphicalGenome(nodefile, edgefile) # Instantiate the path
path = genome.entityPath("A") # Returns the node, edge, node... list that defines the path through the entire genome for the "A" strain
```
### `findEdge`
Given a node and a strain, return the incoming or outgoing edge that contains the strain. In the case of CC strains this will also look for cases of heterozygosity and potentially return more than one edge that contains each heterozygous strain. 
* Parameters:
	* (str) `node` - The node who's edge you want to examine
	* (str) `strain` - the specified strain you wish to find on the edge
* Optional Parameters:
	* (str) `origin` - If `origin` == "src" then return outgoing edges with the specified strain on them. If `origin` == "dst" then return incoming edges with the specified strain on them. 
		* Defaults to "src".
		* Must be an option either in "src" or "dst".
* Output:	
	* Return an edge or multiple edges that contains the specified strain. 
```python
import CCGG
genome = CCGG.GraphicalGenome(nodefile, edgefile)
edges = genome.findEdge(<node>, <strain>, src="dst") # return the incoming edge/edges that contain the desired strain
```

### `anchorPosList` (Potentially Internal Method?)
populates global field: self.sortednodes. It will return a dictionary of anchors as keys with 8 [founder] values that show their.

### `boundingAnchor`
Returns the bouding anchors of the given strain and coordinate positions. 
* Parameters:
	* (str) `strain` - Founder or CC path that you want to follow
	* (int) `start` - Linear coordinate representing desired left bound
	* (int) `end` - Linear coordinate representing desired right bound
* Output:
	* Returns the two bounding anchors that most closely contain the start and end linear coordinates provided. 
```python
import CCGG
genome = CCGG.GraphicalGenome(nodefile, edgefile)
leftBoundingAnchor, rightBoundingAnchor = genome.boudingAnchor("B", <start>, <end>)
```

### `findNumOfFloatingNodes`
Returns the number of FloatingNodes between anchor pairs
* Parameters:
	* (str) `start` - The starting anchor
	* (str) `end` - The ending anchor
* Output
	Returns the number of floating nodes between two anchors
```python
import CCGG
genome = CCGG.GraphicalGenome(nodefile, edgefile)
numOfFloatingNodes = genome.findNumOfFloatingNodes(<startAnchor>, <endAnchor>)
```

### `findsequence`
Given a [node, edge, node, edge ....] list defining a path we construct the sequence on the provided path
* Parameters:
	* (list[str]) `pathlist` - The list of [node, edge, node...] names that define a path through the genome
* Optional Parameters:
	* (bool) `countinganchor` - Flag that states whether to count the first and last element in the `pathlist`
* Output:
	* Returns the sequence for the input `pathlist`
```python
import CCGG
genome = CCGG.GraphicalGenome(nodefile, edgefile)
pathlist = [<node1>, <edge1>, <node2>, <edge2>, <node3>, ...]
sequence = genome.findsequence(pathlist)
```

### returnNodeEdgePath
Given a strain, return an alternating list of edges and nodes that traces the strain through the graphical genome. 
* Parameters:
	* (str) `strain` - Strain you wish to trace
	* (str) `chromo` - Chromosome you are tracing the strain through
* Output:
	* list(str) - List of edges and nodes starting with a "sink" edge and ending with a "sink" edge
```python
import CCGG
genome = CCGG.GraphicalGenome(nodefile, edgefile)
tracepath = genome.returnNodeEdgePath(<strain>, <chromo>)
```

### DemoteAnchor
Given an anchor node, convert it into a sequence of floating node -> edge -> floating node
* Parameters:
	* (str) `anchor` - The anchor being demoted
* Output:
	* Modifies internal dictionary. Does not return a value
```python
import CCGG
genome = CCGG.GraphicalGenome(nodefile, edgefile)
genome.DemoteAnchor(<anchor>)
```

## Example Usage
### Constant values we will use
```python
nodefile = "EdgeFile_Chr17.fa"
edgefile = "NodeFile_Chr17.fa"
```