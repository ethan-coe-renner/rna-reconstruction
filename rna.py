# RNA reconstruction
# Ethan Coe-Renner

from collections import defaultdict


# frag is a fragment, search is a list of the characters to break at, returns list of fragments
def break_frag(frag, search):
    fragments = [""]
    c = 0
    for singleteton in frag:
        fragments[c] += singleteton
        if singleteton in search:
            c += 1
            fragments.insert(c, "")
    if fragments[-1] == "":
        return fragments[:-1]

    return fragments


# takes list of fragments and breaks each by enzyme (a list of singletons to break at), returning a list of fragments
def break_frags(frags, enzyme):
    total = []
    for frag in frags:
        total += break_frag(frag, enzyme)
    return total


#takes the cu digest and g digest, returns the one or maybe two abnormal fragments
def find_abnormals(cu_frags, g_frags):
    abnormals = []
    for frag in cu_frags:
        if frag[-1] in ['g', 'a']:
            abnormals.append(frag)

    for frag in g_frags:
        if frag[-1] in ['c', 'u', 'a']:
            abnormals.append(frag)

    return abnormals


# return except of two lists
def exceptList(lst1, lst2):
    for i in lst2:
        for ii in lst1:
            if i == ii:
                lst1.remove(ii)
    return lst1


# return an extended base from a fragment, if it exists
def get_midextended_base(fragment, enzyme):
    bases = break_frag(fragment, enzyme)
    if len(bases) < 3:
        return None
    return bases[1:-1]


# return extended bases from digests
def find_mid_extended_bases(cu_digest, g_digest):
    extended_bases = []
    for frag in cu_digest:
        meb = get_midextended_base(frag, ['g'])
        if meb:
            extended_bases += meb

    for frag in g_digest:
        meb = get_midextended_base(frag, ['c', 'u'])
        if meb:
            extended_bases += meb

    return extended_bases


# get singleton if it is a singleton
def is_singleton(fragment, enzyme):
    bases = break_frag(fragment, enzyme)
    if len(bases) == 1:
        return bases[0]
    return None


def get_singletons(cu_digest, g_digest):
    singletons = []

    for frag in cu_digest:
        singleton = is_singleton(frag, ['g'])
        if singleton:
            singletons.append(singleton)

    for frag in g_digest:
        singleton = is_singleton(frag, ['c', 'u'])
        if singleton:
            singletons.append(singleton)

    return singletons


# start_and_end is a list of two bases, one is the start
def get_end(start_and_end, abnormals):
    end = None
    for base in start_and_end:
        for abnormal in abnormals:
            if base in abnormal:
                end = abnormal
                break
        if not end:
            break
    return end


def get_start_and_end(uc_digest, g_digest):
    singletons = get_singletons(uc_digest, g_digest)

    midextended_bases = find_mid_extended_bases(uc_digest, g_digest)

    print("singletons: ", singletons)
    print("midextended_bases: ", midextended_bases)

    start_and_end = exceptList(singletons, midextended_bases)

    abnormals = find_abnormals(uc_digest, g_digest)

    print("abnormals: ", abnormals)

    end = get_end(start_and_end, abnormals)
    print("end: ", end)
    start = start_and_end[0] if len(
        start_and_end) < 2 or end in start_and_end[1] else start_and_end[1]
    print("start: ", start)

    return (start, end)


def get_interior_vertices(uc_digest, g_digest):

    vertices = []

    for frag in uc_digest:
        bases = break_frag(frag, ['g'])
        if len(bases) >= 3:
            vertices.append(bases[0])
            vertices.append(bases[-1])

    for frag in g_digest:
        bases = break_frag(frag, ['u', 'c'])
        if len(bases) >= 3:
            vertices.append(bases[0])
            vertices.append(bases[-1])

    return list(set(vertices))


def add_end(vertices, end):
    for i in vertices:
        if i in end:
            vertices.remove(i)
            vertices.append(i)
            return

    vertices.append(end)


def get_vertices(uc_digest, g_digest):
    vertices = get_interior_vertices(uc_digest, g_digest)

    start, end = get_start_and_end(uc_digest, g_digest)

    vertices.insert(0, start)
    vertices.insert(0, "*")

    add_end(vertices, end)
    return vertices


def get_edge(fragment, enzyme):
    bases = break_frag(fragment, enzyme)
    if len(bases) == 2:
        return (bases[0], bases[1])
    elif len(bases) >= 3:
        return (bases[0], "".join(bases[1:-1]), bases[-1])  # (v1, edge, v2)
    else:
        return None


class Graph:

    def __init__(self, vertices):
        print(vertices)
        self.vertices = vertices
        self.V = len(vertices)  #No. of vertices
        self.graph = defaultdict(list)  # default dictionary to store graph

    def create_graph(self, uc_digest, g_digest):
        self.add_edge(0, 1, "")  # add start edge

        for frag in uc_digest:
            edge = get_edge(frag, ['g'])
            if not edge:
                continue
            if len(edge) == 2:
                self.add_edge(self.vertices.index(edge[0]),
                              self.vertices.index(edge[1]), "")
            elif len(edge) == 3:
                self.add_edge(self.vertices.index(edge[0]),
                              self.vertices.index(edge[2]), edge[1])

        for frag in g_digest:
            edge = get_edge(frag, ['u', 'c'])
            if not edge:
                continue
            if len(edge) == 2:
                self.add_edge(self.vertices.index(edge[0]),
                              self.vertices.index(edge[1]), "")
            elif len(edge) == 3:
                self.add_edge(self.vertices.index(edge[0]),
                              self.vertices.index(edge[2]), edge[1])

        self.add_edge(self.V - 1, 0, "")  # add ending edge

    # function to add an edge to graph
    def add_edge(self, u, v, base):
        self.graph[u].append((v, base))

    # This function removes edge u-v from graph
    def remove_edge(self, u, v):
        for index, key in enumerate(self.graph[u]):
            if key[0] == v:
                x = self.graph[u][index][1]
                self.graph[u].pop(index)
                return x
        print("returning none from rmv")
        return ""

    def DFSCount(self, v, visited):
        count = 1
        visited[v] = True
        for i in self.graph[v]:
            # print(i)
            if visited[i[0]] == False:
                count = count + self.DFSCount(i[0], visited)
        return count

    # The function to check if edge u-v can be considered as next edge in
    # Euler Tour
    def is_valid_next_edge(self, u, v):
        if len(self.graph[u]) == 1:
            return True
        else:
            visited = [False] * (self.V)
            count1 = self.DFSCount(u, visited)

            base = self.remove_edge(u, v)
            visited = [False] * (self.V)
            count2 = self.DFSCount(u, visited)

            self.add_edge(u, v, base)

            return False if count1 > count2 else True

    # Print Euler tour starting from vertex u
    def print_euler_util(self, u):
        for v in self.graph[u]:
            if self.is_valid_next_edge(u, v[0]):
                print(f"{self.vertices[u]}{v[1]}", end=""),
                self.remove_edge(u, v[0])
                self.print_euler_util(v[0])

    def print_euler_cycle(self):
        #Find a vertex with odd degree
        u = 0
        for i in range(self.V):
            if len(self.graph[i]) % 2 != 0:
                u = i
                break
        # Print tour starting from odd vertex
        print("\n")
        self.print_euler_util(u)
        print("*")


def main():
    print("RNA Reconstruction")
    response = input("Would you like to see an example (y/n)?")
    if response != 'n':
        example()
        return
    print("Enter uc fragments, just enter <return> to finish")
    uc_digest = []
    fragment = input("fragment: ")
    while fragment != "":
        uc_digest.append(fragment)
        print("fragments so far in uc_digest: ", uc_digest)
        fragment = input("fragment: ")

    print("Enter g fragments, just enter <return> to finish")
    g_digest = []
    fragment = input("fragment: ")
    while fragment != "":
        g_digest.append(fragment)
        print("fragments so far in g_digest: ", g_digest)
        fragment = input("fragment: ")

    v = get_vertices(uc_digest, g_digest)
    g = Graph(v)
    g.create_graph(uc_digest, g_digest)

    print("Reconstruction: ")
    g.print_euler_cycle()


def example():
    # HW 7 prob 1
    g_digest = ["ag", "aag", "ucucag"]
    uc_digest = ["c", "c", "u", "aagu", "agag"]

    # uc_digest = ['c', 'c', 'agu', 'gagu', 'ggau', 'agu']
    # g_digest = ['g', 'u', 'ag', 'aug', 'uag', 'uccg']

    print("Example:")
    print("UC-Digest:", uc_digest)
    print("G-Digest:", g_digest)

    v = get_vertices(uc_digest, g_digest)
    g = Graph(v)
    g.create_graph(uc_digest, g_digest)

    g.print_euler_cycle()


main()
