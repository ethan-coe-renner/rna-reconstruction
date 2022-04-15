# RNA Reconstruction
## Description
For my part of the project, I implemented the RNA reconstruction algorithm in Python. The program will ask if the user wants to see an example. If the example is chosen, the program will run through a hardcoded set of digests from problem 1 of Homework 7. Otherwise, the user will be prompted to enter the digests manually and then the program will find the reconstruction.
	The program runs through the algorithm as discussed in class, first by finding the start and end. It does this by finding abnormal fragments, singletons and extended bases in the middle of fragments. Then, it finds all the vertices in a similar way as it found the extended bases, along with the edge strings. Using these, it constructs a graph using the graph class I defined. This graph is modeled as a list of vertices and a dictionary. The dictionary has vertices as keys, and values which are lists of edges. Each edge is a tuple of a number, representing the destination vertex, and a string, representing the strings on that edge. Once the graph has been constructed, the program tries to find an eulerian path through the graph, and prints it out as the reconstructed string. The program uses Fleury’s Algorithm to find the eulerian path. This algorithm works by making sure there are either 0 or 2 odd vertices (vertices with odd degree). Then, choose a vertex to start at. This vertex must either be one of the two odd vertices if they exist, or any vertex if they don’t. Then, just follow edges one at a time, and choose the non-bridge if there is a choice. The program prints out the edges as it follows them.
	The draft implementation currently has at least one flaw, it seems to not properly handle graphs with multiple edges with the same start and end. This is a problem that I intend to fix by the final draft.

## How to Run
In order to run this project, simply run `python rna.py` in this directory
