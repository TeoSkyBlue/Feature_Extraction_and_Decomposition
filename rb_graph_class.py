class ReebGraph:
    def __init__(self, index = None, parents = None):
        self.index = index
        self.left_bound = []
        self.right_bound = []
        self.parents = parents
        self.Tsets = []

class Rnode:
    def __init__(self, index = None, parent = None, a = None, l = None):
        self.index = index
        self.parent = parent
        self.Tsets = []
        self.left_bound = []
        self.right_bound = []
        self.attribute = attributeElement(a, l)
        self.children = []
        self.MLIST = []


class attributeElement:
    def __init__(self, a = None, l = None):
        self.a = a
        self.l = l 



class MPairElement:
    def __init__(self, node1 = None, node2 = None, vector_nodes1 = [], vector_nodes2 = []):
        self.node1 = node1 #node from MRG1
        self.node2 = node2 #node from MRG2
        self.vector_nodes1 = vector_nodes1
        self.vector_nodes2 = vector_nodes2

