#!/usr/bin/python

from pylab import inf

class SparseGraph(object):
    
    def __init__(self):
        self._data = {}

    def __len__(self):
        return len(self._data)
        
    def __contains__(self, v):
        if (v in self._data):
            return True
        for (v1, neigh) in self._data.iteritems():
            if (v in neigh):
                return True
        return False

    def __getitem__(self, v):
        if (not v in self._data):
            self._data[v] = {}
        return self._data[v]

    def __eq__(self, other):
        if not isinstance(other, SparseGraph):
            return False
        return self._data == other._data
    
    def __ne__(self, other):
        if not isinstance(other, SparseGraph):
            return True
        return self._data != other._data

    def __hash__(self):
        return hash(self._data)

    def itervertices(self):
        vertices = set()
        for (v1, v2, w) in self:
            vertices.add(v1)
            vertices.add(v2)
        
        for v in vertices:
            yield v

    def __iter__(self):
        for (v1, neigh) in self._data.iteritems():
            for (v2, weight) in neigh.iteritems():
                yield (v1, v2, weight)
                
    def __repr__(self):
        return "\n".join(["(%s, %s) - %s" % (str(v1), str(v2), str(w)) for (v1, v2, w) in self])
    
    def __delitem__(self, v):
        if (v in self._data):
            del self._data[v]
        for (v1, neigh) in self._data.iteritems():
            if (v in neigh):
                del neigh[v]
    
    # Dijkstra's algorithm for shortest paths
    # David Eppstein, UC Irvine, 4 April 2002
    # http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/117228
    def Dijkstra(self, start, end=None):
        """
        Find shortest paths from the start vertex to all
        vertices nearer than or equal to the end.
    
        The input graph G is assumed to have the following
        representation: A vertex can be any object that can
        be used as an index into a dictionary.  G is a
        dictionary, indexed by vertices.  For any vertex v,
        G[v] is itself a dictionary, indexed by the neighbors
        of v.  For any edge v->w, G[v][w] is the length of
        the edge.  This is related to the representation in
        <http://www.python.org/doc/essays/graphs.html>
        where Guido van Rossum suggests representing graphs
        as dictionaries mapping vertices to lists of neighbors,
        however dictionaries of edges have many advantages
        over lists: they can store extra information (here,
        the lengths), they support fast existence tests,
        and they allow easy modification of the graph by edge
        insertion and removal.  Such modifications are not
        needed here but are important in other graph algorithms.
        Since dictionaries obey iterator protocol, a graph
        represented as described here could be handed without
        modification to an algorithm using Guido's representation.
    
        Of course, G and G[v] need not be Python dict objects;
        they can be any other object that obeys dict protocol,
        for instance a wrapper in which vertices are URLs
        and a call to G[v] loads the web page and finds its links.
        
        The output is a pair (D,P) where D[v] is the distance
        from start to v and P[v] is the predecessor of v along
        the shortest path from s to v.
        
        Dijkstra's algorithm is only guaranteed to work correctly
        when all edge lengths are positive. This code does not
        verify this property for all edges (only the edges seen
         before the end vertex is reached), but will correctly
        compute shortest paths even for some graphs with negative
        edges, and will raise an exception if it discovers that
        a negative edge has caused it to make a mistake.
        """
        
        from priodict import priorityDictionary
    
        D = {}    # dictionary of final distances
        P = {}    # dictionary of predecessors
        Q = priorityDictionary()   # est.dist. of non-final vert.
        Q[start] = 0
        
        for v in Q:
            D[v] = Q[v]
            if v == end: break
            
            for w in self[v]:
                vwLength = D[v] + self[v][w]
                if w in D:
                    if vwLength < D[w]:
                        raise ValueError, "Dijkstra: found better path to already-final vertex"
                elif w not in Q or vwLength < Q[w]:
                    Q[w] = vwLength
                    P[w] = v
        
        return (D, P)
                
    def shortest_path(self, start, end):
        """
        Find a single shortest path from the given start vertex
        to the given end vertex.
        The input has the same conventions as Dijkstra().
        The output is a list of the vertices in order along
        the shortest path.
        """
    
        (D, P) = self.Dijkstra(start, end)
        Path = []
        while 1:
            Path.append(end)
            if end == start: break
            if (not end in P):
                return None
            end = P[end]
        Path.reverse()
        return Path
    
    def shortest_distance(self, start, end):
        (D, P) = self.Dijkstra(start, end)
        if (end in D):
            return D[end]
        else:
            return inf
    
    def find_all_paths(self, s, t, max_weight=inf, max_length=inf, forbidden_vertices=set([])):
        """loop-avoiding DFS that stops at a certain depth and/or total-weight,
           and returns all possible paths from 's' to 't'.
           G is in a sparse graph representation.
           Uses recursion for the DFS.
        """
        if (max_weight <= 0 or max_length <= 0):
            return []
        if (s == t):
            return [[s]]
        
        paths = []
        for (i, weight) in self[s].iteritems():
            if (not i in forbidden_vertices):
                for path in self.find_all_paths(i, t, max_weight - weight, max_length - 1, forbidden_vertices | set([s])):
                    paths.append([s] + path)
        return paths
    
    def to_string_dfs(self, r, width=4, prefix="", connector=" ", predecessors=set()):
        """Returns a tree representation of the graph, while doing a DFS
           Assumes no loops in the graph.
        """

        next_predecessors = predecessors | set([r])
        children = list(set(self[r]) - next_predecessors) # avoid loops!

        
        if (children != []):    
            string = prefix + '+' + '-' * width + '+' + ' ' + str(r) + "\n"
            string += prefix + connector + " " * width + "|" + "\n"
            while (len(children) > 1):
                string += self.to_string_dfs(children.pop(), width, prefix + connector + " " * width, "|", next_predecessors)
            string += self.to_string_dfs(children.pop(), width, prefix + connector + " " * width, " ", next_predecessors)
        else:
            string = prefix + '+' + '-' * width + '-' + ' ' + str(r) + "\n"
            string += prefix + connector + "\n"
        
        return string

def test():
    G = SparseGraph()
    G[0][1] = 1
    G[0][3] = 0
    G[1][2] = 1
    G[1][1] = 4
    G[2][5] = 0
    G[3][4] = 1
    G[3][5] = 1
    G[4][5] = 10
    G[5][1] = 1
    G[5][2] = 1

    print list(G.itervertices())
    print "G = \n" + str(G)
    
    
    for s in range(0, 6):
        (D, P) = G.Dijkstra(s)
        print "********** s = %d, D = %s, P = %s" % (s, str(D), str(P))
        for t in range(0, 6):
            if (not t in D):
                print "t = %d : cannot be reached" % t
            else:
                print "t = %d : weight(%s) = %s" % (t, str(G.shortest_path(s, t)), D[t])
  
    print G.to_string_dfs(0)
    
    print "[DONE]"

    
    
if __name__ == '__main__': test()
