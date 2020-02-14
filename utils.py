# Some utility classes to represent the aa sequence

class AASeq:
    """
    A simple class for an amino acid sequence
    """

    def __init__(self, name):
        self.name = name
        self.partialsequence = []

    # Overload the __repr__ operator to make printing simpler.
    def __repr__(self):
        return self.name
