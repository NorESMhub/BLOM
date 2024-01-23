class CaseDict(dict):
    """Class to handle using a case as a dictionary rather 
    than a CIME case object"""

    def get_value(self, key):
        if key in self:
            return self[key]
        else:
            return None
