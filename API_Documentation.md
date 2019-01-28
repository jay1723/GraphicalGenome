# Graphical Genome Class

## Graphical Genome basic usage
The Graphical Genome, when initialized, returns an object with four instance data dictionaries. 

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
GraphicalGenome.edges["S01.A"] = {"src" : "SOURCE", "seq" : <sequence>, "dst" : <node>, "strain" : ["A"]}
```

```python
GraphicalGenome.edges["K01.A"] = {"src" : <node>, "seq" : <sequence>, "dst" : "SINK", "strain" : ["A"]}
```
## API

### Initialization
Initalize a new graphical genome instance
	* Parameters:
		* (str) `nodefile` - Absolute path to the nodefile
    	* (str) `edgefile` - Absolute path to the edgefile
    * Output: new instance of Graphical Genome
```python
newinstance = GraphicalGenome(nodefile, edgefile) # Initialize Graphical Genome instance with the specified node and edge file
newinstance.edges[<edgename>]["strains"] # Example of looking at the list of strains given an edgename
newinstance.nodes[<nodename>]["repeatclass"] # Example of returning the list of repeat classes of a given node
```

### writeFasta
Write updated node and edge dictionaries in a FASTA format.
	* Parameters:
		* (str) `filename` - Absolute path for the file that you would like to write. Files that do not exist will be created in the specified location with the specified name. 
		* (dict) `file_dict` - Node or Edge dictionary that you would like to write. 
	* Output: FASTA file created in the specified location
```python
newinstance = GraphicalGenome(nodefile, edgefile) # Initialize Graphical Genome instance with the specified node and edge file
newinstance.nodes[<nodename>]["new_annotation"] = <new_data> # Add a new annotation to the specified node
newinstance.writeFasta("FASTAFolder/newfile.fa", newinstance.nodes) # Writes updated node dictionary to the specified FASTA file
```

### showParallelEdges
Show all edges that share source and destination nodes
	* Parameters: 
		* (str) `edge` - Name of the edge
	* Output: 
		* list(str) - list of edges that share source and destination nodes
```python
newinstance = GraphicalGenome(nodefile, edgefile) # Initialize Graphical Genome instance with the specified node and edge file
parallelledges = newinstance.showParallelEdges("E01.00000001") # Returns a list of all edges that share the same source and destination as "E01.00000001"
```
### deleteNode
Delete the given node while also resolving the edges connected to or originating from that node
	* Parameters:
		* (str) `node` - Name of the node
	* Output:
		* Updated node and edge dictionaries without the specified node.
```python
newinstance = GraphicalGenome(nodefile, edgefile) # Initialize Graphical Genome instance with the specified node and edge file
newinstance.deleteNode("A01.00000001") # Deletes node "A01.00000001" from the nodes dictionary and redirects the edges that start and end in the node. 
```
### returnNodeEdgePath
Given a strain, return an alternating list of edges and nodes that traces the strain through the graphical genome. 
	* Parameters:
		* (str) `strain` - Strain you wish to trace
		* (str) `chromo` - Chromosome you are tracing the strain through
	* Output:
		* list(str) - List of edges and nodes starting with a "sink" edge and ending with a "sink" edge
```python
newinstance = GraphicalGenome(nodefile, edgefile) # Initialize Graphical Genome instance with the specified node and edge file
tracepath = newinstance.returnNodeEdgePath("A", "Chr1") # Returns the nodes and edges that make up the A strain on chromosome 1
```
### findMissingPath
Given a CC path, look for all node pairs that do not contain an edge between them that has the specified CC strain on them. 
	* Parameters:
		* (str) `strain` - Strain you wish to trace
		* (str) `chromo` - Chromosome you are tracing the strain through
	* Output:
		* list(tuple(node, node)) - List of node tuples that do not have any edges between them with the specified CC strain on them. 
```python
newinstance = GraphicalGenome(nodefile, edgefile) # Initialize Graphical Genome instance with the specified node and edge file
missingnodes = newinstance.findMissingPath("CC027", "Chr1") # Returns node tuples that do not have any edge that has a portion of "CC027" on it
```
### reconstructSequence
Reconstruct the sequence of a given strain
	* Parameters:
		* (str) `strain` - Strain you wish to trace
		* (str) `chromo` - Chromosome you are tracing the strain through
	* Output:
		* (str) - Reconstructed sequence for the given strain
```python
newinstance = GraphicalGenome(nodefile, edgefile) # Initialize Graphical Genome instance with the specified node and edge file
seq = newinstance.reconstructSequence("A", "Chr1") # Reconstruct strain A on chromosome 1
```

## Example Usage

### Reading and writing a file

### Reconstructing a sequence and writing the full sequence to a FASTA file