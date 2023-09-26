class ReebGraph:
    def __init__(self, index = None, parents = None, Tsets = None):
        self.index = index
        self.left_bound = []
        self.right_bound = []
        self.parents = parents
        self.Tsets = Tsets

