"""
Chromosome and sex chromosome information for common species used in GWAS.
"""

class Chromosomes:
    """
    Class to manage chromosome numbers and sex chromosomes for common species.
    Uses self.chromosomes to store all chromosome identifiers.
    """
    
    # Species chromosome data: (num_autosomes, sex_chromosomes, has_mitochondrial)
    # sex_chromosomes can be empty list, ["X", "Y"], or ["Z", "W"]
    _SPECIES_DATA = {
        "homo sapiens": (22, ["X", "Y"], True),
        "human": (22, ["X", "Y"], True),
        "mus musculus": (19, ["X", "Y"], True),
        "mouse": (19, ["X", "Y"], True),
        "rattus norvegicus": (20, ["X", "Y"], True),
        "rat": (20, ["X", "Y"], True),
        "gallus gallus": (28, ["Z", "W"], True),
        "chicken": (28, ["Z", "W"], True),
        "danio rerio": (25, [], True),
        "zebrafish": (25, [], True),
        "drosophila melanogaster": (None, ["X", "Y"], True),  # Special case: autosomes are ["2", "3", "4"]
        "fruit fly": (None, ["X", "Y"], True),
        "sus scrofa": (18, ["X", "Y"], True),
        "pig": (18, ["X", "Y"], True),
        "bos taurus": (29, ["X", "Y"], True),
        "cattle": (29, ["X", "Y"], True),
        "cow": (29, ["X", "Y"], True),
        "canis lupus familiaris": (38, ["X", "Y"], True),
        "dog": (38, ["X", "Y"], True),
        "equus caballus": (31, ["X", "Y"], True),
        "horse": (31, ["X", "Y"], True),
        "oryza sativa": (12, [], True),
        "rice": (12, [], True),
        "arabidopsis thaliana": (None, [], True),  # Special case: autosomes are ["1", "2", "3", "4", "5"]
        "arabidopsis": (None, [], True),
    }
    
    def __init__(self, species="homo sapiens"):
        """
        Initialize chromosome information for a given species.
        
        Parameters:
        -----------
        species : str, default="homo sapiens"
            Species name (case-insensitive)
        """
        self.species = species.lower()
        self.chromosomes = list()
        self.autosomes = list()
        self.sex_chromosomes = list()
        self.mitochondrial = None
        
        self._initialize_species_data()
    
    def _initialize_species_data(self):
        """Initialize chromosome data based on species."""
        species_key = self.species
        
        # Get species data or use default (human)
        if species_key not in self._SPECIES_DATA:
            species_key = "homo sapiens"
        
        num_autosomes, sex_chr_list, has_mt = self._SPECIES_DATA[species_key]
        
        # Set autosomes
        if num_autosomes is None:
            # Special cases with non-sequential autosomes
            if species_key in ["drosophila melanogaster", "fruit fly"]:
                self.autosomes = ["2", "3", "4"]
            elif species_key in ["arabidopsis thaliana", "arabidopsis"]:
                self.autosomes = ["1", "2", "3", "4", "5"]
            else:
                self.autosomes = []
        else:
            self.autosomes = [str(i) for i in range(1, num_autosomes + 1)]
        
        # Set sex chromosomes
        self.sex_chromosomes = sex_chr_list.copy()
        
        # Set mitochondrial
        self.mitochondrial = "MT" if has_mt else None
        
        # Build complete chromosome list
        self.chromosomes = self.autosomes.copy()
        if self.sex_chromosomes:
            self.chromosomes.extend(self.sex_chromosomes)
        if self.mitochondrial:
            self.chromosomes.append(self.mitochondrial)
    
    def get_sex_chromosomes_numeric(self, xymt_num=[23, 24, 25]):
        """
        Get sex chromosomes as numeric values (for compatibility with existing code).
        
        Parameters:
        -----------
        xymt_num : list, default=[23, 24, 25]
            Numeric values for X, Y, MT chromosomes (human convention)
            
        Returns:
        --------
        list
            List of numeric sex chromosome identifiers
        """
        if len(self.sex_chromosomes) >= 2:
            return [xymt_num[0], xymt_num[1]]  # First two sex chromosomes
        elif len(self.sex_chromosomes) == 1:
            return [xymt_num[0]]  # Single sex chromosome
        else:
            return []
    
    def get_all_sex_chromosomes_numeric(self, xymt_num=[23, 24, 25]):
        """
        Get all non-autosomal chromosomes (sex + mitochondrial) as numeric values.
        
        Parameters:
        -----------
        xymt_num : list, default=[23, 24, 25]
            Numeric values for X, Y, MT chromosomes (human convention)
            
        Returns:
        --------
        list
            List of numeric non-autosomal chromosome identifiers
        """
        sex_chr_numeric = self.get_sex_chromosomes_numeric(xymt_num)
        if self.mitochondrial:
            return sex_chr_numeric + [xymt_num[2]]  # Add MT
        return sex_chr_numeric
    
    def get_chromosome_mappings(self, xymt_num=[23, 24, 25]):
        """
        Get chromosome mappings as tuples (label, numeric_value) for x, y, mt.
        
        Parameters:
        -----------
        xymt_num : list, default=[23, 24, 25]
            Numeric values for X, Y, MT chromosomes (human convention)
            
        Returns:
        --------
        tuple
            (x_tuple, y_tuple, mt_tuple) where each tuple is (label, numeric_value)
        """
        sex_chr_numeric = self.get_sex_chromosomes_numeric(xymt_num)
        all_sex_chr_numeric = self.get_all_sex_chromosomes_numeric(xymt_num)
        
        # Get x mapping (first sex chromosome)
        if len(self.sex_chromosomes) >= 1:
            x = (self.sex_chromosomes[0], sex_chr_numeric[0] if len(sex_chr_numeric) > 0 else 23)
        else:
            x = ("X", 23)  # Default, won't be used if no sex chromosomes
        
        # Get y mapping (second sex chromosome)
        if len(self.sex_chromosomes) >= 2:
            y = (self.sex_chromosomes[1], sex_chr_numeric[1] if len(sex_chr_numeric) > 1 else 24)
        else:
            y = ("Y", 24)  # Default, won't be used if only one or no sex chromosomes
        
        # Get mt mapping (mitochondrial)
        if self.mitochondrial:
            mt_num = all_sex_chr_numeric[-1] if len(all_sex_chr_numeric) > len(sex_chr_numeric) else 25
            mt = (self.mitochondrial, mt_num)
        else:
            mt = ("MT", 25)  # Default, won't be used if no mitochondrial
        
        return x, y, mt
    
    def get_min_chromosome(self):
        """
        Get the minimum autosome number.
        
        Returns:
        --------
        int
            Minimum autosome number, defaults to 1 if no autosomes
        """
        if not self.autosomes:
            return 1
        
        try:
            numeric_autosomes = [int(a) for a in self.autosomes if a.isdigit()]
            return min(numeric_autosomes) if numeric_autosomes else 1
        except:
            return 1
    
    def is_sex_chromosome(self, chrom):
        """
        Check if a chromosome identifier is a sex chromosome.
        
        Parameters:
        -----------
        chrom : str or int
            Chromosome identifier
            
        Returns:
        --------
        bool
            True if chromosome is a sex chromosome
        """
        chrom_str = str(chrom).upper()
        return chrom_str in [c.upper() for c in self.sex_chromosomes]
    
    def is_autosome(self, chrom):
        """
        Check if a chromosome identifier is an autosome.
        
        Parameters:
        -----------
        chrom : str or int
            Chromosome identifier
            
        Returns:
        --------
        bool
            True if chromosome is an autosome
        """
        chrom_str = str(chrom)
        return chrom_str in self.autosomes
    
    def is_mitochondrial(self, chrom):
        """
        Check if a chromosome identifier is mitochondrial.
        
        Parameters:
        -----------
        chrom : str or int
            Chromosome identifier
            
        Returns:
        --------
        bool
            True if chromosome is mitochondrial
        """
        chrom_str = str(chrom).upper()
        return chrom_str == self.mitochondrial.upper() if self.mitochondrial else False
    
    def get_chr_to_number_dict(self, out_chr=False, xymt_num=[23, 24, 25], max_chr=200):
        """
        Create a dictionary mapping chromosome identifiers to numeric representations.
        
        Parameters:
        -----------
        out_chr : bool, default=False
            If True, returns dictionary with string keys and values
        xymt_num : list, default=[23, 24, 25]
            Numeric values for X, Y, MT chromosomes (human convention)
        max_chr : int, default=200
            Maximum chromosome number to include in dictionary
            
        Returns:
        --------
        dict
            Dictionary mapping chromosome identifiers to numeric values or strings
        """
        sex_chr_numeric = self.get_sex_chromosomes_numeric(xymt_num)
        all_sex_chr_numeric = self.get_all_sex_chromosomes_numeric(xymt_num)
        
        if out_chr:
            dic = {str(i): str(i) for i in range(1, max_chr + 1)}
            # Map sex chromosomes
            if len(self.sex_chromosomes) >= 1:
                dic[self.sex_chromosomes[0]] = str(sex_chr_numeric[0]) if len(sex_chr_numeric) > 0 else "23"
            if len(self.sex_chromosomes) >= 2:
                dic[self.sex_chromosomes[1]] = str(sex_chr_numeric[1]) if len(sex_chr_numeric) > 1 else "24"
            # Map mitochondrial
            if self.mitochondrial:
                mt_num = all_sex_chr_numeric[-1] if len(all_sex_chr_numeric) > len(sex_chr_numeric) else 25
                dic[self.mitochondrial] = str(mt_num)
                dic["M"] = str(mt_num)  # Also support "M" as alias
        else:
            dic = {str(i): i for i in range(1, max_chr + 1)}
            # Map sex chromosomes
            if len(self.sex_chromosomes) >= 1:
                dic[self.sex_chromosomes[0]] = sex_chr_numeric[0] if len(sex_chr_numeric) > 0 else 23
            if len(self.sex_chromosomes) >= 2:
                dic[self.sex_chromosomes[1]] = sex_chr_numeric[1] if len(sex_chr_numeric) > 1 else 24
            # Map mitochondrial
            if self.mitochondrial:
                mt_num = all_sex_chr_numeric[-1] if len(all_sex_chr_numeric) > len(sex_chr_numeric) else 25
                dic[self.mitochondrial] = mt_num
                dic["M"] = mt_num  # Also support "M" as alias
        
        return dic
    
    def get_number_to_chr_dict(self, in_chr=False, xymt_num=[23, 24, 25], prefix="", max_chr=200):
        """
        Create a dictionary mapping chromosome numbers to string representations.
        
        Parameters:
        -----------
        in_chr : bool, default=False
            If True, returns dictionary with string keys and values
        xymt_num : list, default=[23, 24, 25]
            Numeric values for X, Y, MT chromosomes (human convention)
        prefix : str, default=""
            Optional prefix for chromosome identifiers
        max_chr : int, default=200
            Maximum chromosome number to include in dictionary
            
        Returns:
        --------
        dict
            Dictionary mapping chromosome numbers to string representations
        """
        sex_chr_numeric = self.get_sex_chromosomes_numeric(xymt_num)
        all_sex_chr_numeric = self.get_all_sex_chromosomes_numeric(xymt_num)
        
        if in_chr:
            dic = {str(i): prefix + str(i) for i in range(1, max_chr + 1)}
            # Map sex chromosomes
            if len(self.sex_chromosomes) >= 1:
                dic[str(sex_chr_numeric[0]) if len(sex_chr_numeric) > 0 else "23"] = prefix + self.sex_chromosomes[0]
            if len(self.sex_chromosomes) >= 2:
                dic[str(sex_chr_numeric[1]) if len(sex_chr_numeric) > 1 else "24"] = prefix + self.sex_chromosomes[1]
            # Map mitochondrial
            if self.mitochondrial:
                mt_num = all_sex_chr_numeric[-1] if len(all_sex_chr_numeric) > len(sex_chr_numeric) else 25
                dic[str(mt_num)] = prefix + self.mitochondrial
        else:
            dic = {i: prefix + str(i) for i in range(1, max_chr + 1)}
            # Map sex chromosomes
            if len(self.sex_chromosomes) >= 1:
                dic[sex_chr_numeric[0] if len(sex_chr_numeric) > 0 else 23] = prefix + self.sex_chromosomes[0]
            if len(self.sex_chromosomes) >= 2:
                dic[sex_chr_numeric[1] if len(sex_chr_numeric) > 1 else 24] = prefix + self.sex_chromosomes[1]
            # Map mitochondrial
            if self.mitochondrial:
                mt_num = all_sex_chr_numeric[-1] if len(all_sex_chr_numeric) > len(sex_chr_numeric) else 25
                dic[mt_num] = prefix + self.mitochondrial
        
        return dic
    
    def get_chr_list(self, add_number=False, only_number=False):
        """
        Generate a list of chromosome identifiers for this species.
        
        Parameters:
        -----------
        add_number : bool, default=False
            If True, include both string and numeric representations
        only_number : bool, default=False
            If True, return only numeric chromosome numbers
            
        Returns:
        --------
        list
            List of chromosome identifiers
        """
        if only_number:
            # Return numeric autosomes only
            numeric_autosomes = [int(a) for a in self.autosomes if a.isdigit()]
            return numeric_autosomes if numeric_autosomes else []
        
        chrom_list = self.chromosomes.copy()
        
        if add_number:
            # Add numeric representations
            numeric_autosomes = [int(a) for a in self.autosomes if a.isdigit()]
            chrom_list = self.chromosomes.copy() + numeric_autosomes
        
        return chrom_list


from typing import Optional
def get_sex_chromosomes(species: str = "homo sapiens") -> 'Chromosomes':
    """
    Convenience function to get chromosome information for a species.
    
    Parameters:
    -----------
    species : str, default="homo sapiens"
        Species name
        
    Returns:
    --------
    Chromosomes
        Chromosomes instance for the specified species
    """
    return Chromosomes(species)
