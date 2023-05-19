'''
Created on 27 Mar 2013

@author: dunbar

These are selection classes for the save method of the AbPDB entity
They are based on the Bio.PDB.PDBIO Selection class

'''

#from ABDB.AB_Utils import regions_tuples,is_interface_region

class select_all(object):
    """
    Default selection (everything) during writing - can be used as base class
    to implement selective output. This selects which entities will be written out.
    """

    def __repr__(self):
        return "<Select all>"

    def accept(self, ob):
        if ob.level == "A":
            return self.accept_atom(ob)
        elif ob.level == "R":
            return self.accept_residue(ob)
        elif ob.level == "C":
            return self.accept_chain(ob)
        elif ob.level == "F":
            return self.accept_fragment(ob)
        elif ob.level == "H":
            return self.accept_holder(ob)
        elif ob.level == "M":
            return self.accept_model(ob)


    def accept_model(self, model):
        """
        Overload this to reject models for output.
        """
        return 1
    
    def accept_holder(self, model):
        """
        Overload this to reject holders for output. (fabs, abchains-holder,agchains-holder)
        """
        return 1


    def accept_chain(self, chain):
        """
        Overload this to reject chains for output.
        """
        return 1

    def accept_fragment(self, fragment):
        """
        Overload this to reject residues for output.
        """
        return 1

    def accept_residue(self, residue):
        """
        Overload this to reject residues for output.
        """
        return 1

    def accept_atom(self, atom):
        """
        Overload this to reject atoms for output.
        """
        return 1
    

