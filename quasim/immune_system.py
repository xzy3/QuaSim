"""
Code related to immune system functionality
AntiBody class representing B-cell either targeted to specific
epitope or naive.
EpitopePos class representing positions of specific epitope
EpitopeVariant class for epitope with specific nucleotides on epitope positions
"""
from quasim.virus import EpitopeVariant


class AntiBody:

    def __init__(self, epitope_var, delay=0):
        assert isinstance(epitope_var, EpitopeVariant)
        self.epitope_var = epitope_var
        self.delay = delay

    def tick(self):
        self.delay -= 1
        if self.is_active():
            self.epitope_var.is_targeted = True

    def is_active(self):
        return self.delay <= 0