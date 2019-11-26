#! /usr/bin/python

import argparse
import CCGG_SH as CCGG

def Main():
    parser = argparse.ArgumentParser(description="Graphical Genome Command line CLI. Type -h for help options")

    # Define the subparsers
    subparsers = parser.add_subparsers(help="defining subparsers for different commands", dest="command")



    """
    Dump sequence to file
    """
    parser_sequence = subparsers.add_parser("extract", help="Extract the desired sequence and dump to a file")
    parser_sequence.add_argument("--nodefile", help="Absolute path to the nodefile", type=str, required=True, dest="nodefile")
    parser_sequence.add_argument("--edgefile", help="Absolute path to the edgefile", type=str, required=True, dest="edgefile")
    parser_sequence.add_argument("--outfile", type=str, required=True, dest="outfile") 
    parser_sequence.add_argument("--sequence", help="Extract the strain for the sequence. 2 required arguments: Strain, filepath", type=str, nargs=2)
    


    """
    Get Sub Path
    Parameters: 
        startAnchor <str> - Starting anchor for the subpath
        endAnchor <str> - Ending anchor for the subpath
        strain <str> - The strain whose path should be followed between the two anchors
        outFile <str> - File that output will be dumped to
    Optional Parameters:
        anchor <bool> - Whether or not to count the first and last anchor in the subpath
        sequence <bool> - Dump the sequence between the two anchors instead of dumping the anchor,edge,anchor... vals
    Output:
        If sequence == true: dumps the sequence to the specified file
        If sequence == false: dumps the [anchor,edge,anchor,edge...] list to the specified file
    """
    parser_subpath = subparsers.add_parser("subpath", help="construct subsequences between each anchor pair then merge them together")
    parser_subpath.add_argument("-nf", "--nodefile", help="Absolute path to the nodefile", type=str, required=True, dest="nodefile")
    parser_subpath.add_argument("-ef", "--edgefile", help="Absolute path to the edgefile", type=str, required=True, dest="edgefile")
    parser_subpath.add_argument("-s", "--start", help="The anchor that you want to start on", type=str, required=True, dest="start")
    parser_subpath.add_argument("-e", "--end", help="The anchor that you want to end on", type=str, required=True, dest="end")
    parser_subpath.add_argument("-o", "--outfile", type=str, help="The file you want to save the output to", required=True, dest="outfile")
    parser_subpath.add_argument("-seq", "--sequence", help="Whether you want to dump the sequence or the anchor,edge,anchor... list", action="store_true", dest="sequence")
    parser_subpath.add_argument("-a", "--anchors", help="Whether or not the first and last anchors are included in the subpath", action="store_false", dest="anchor")


# Parse the arguments
    args = parser.parse_args()
    command = args.command
    if (command == "extract"):
        print("hello")
    elif (command == "subpath"):
        nf = args.nodefile
        ef = args.edgefile
        s = args.start 
        e = args.end
        out = args.outfile
        seq = args.sequence
        anc = args.anchor
        print(nf, ef, s, e, out, seq, anc)
        g = CCGG.GraphicalGenome(nf, ef, verbose=True).getSubPath(s, e, seq=seq, counting_anchor=anc)
        print(g)
        
    else:
        print("Command not found, please try again with appropriate command. Choices are ['extract', 'subpath']")

if __name__ == '__main__':
    Main()
    

    # graph = CCGG.GraphicalGenome(args.nodefile, args.edgefile)
    # strain = args.extract[0]
    # filepath = args.extract[1]
    # chromo = args.chromosome
    # graph.writeSequence(strain, chromo, filepath) 