import os 
import re
import io 
import sys
import urllib
import tempfile
import requests
import collections
import numpy as np
import matplotlib
from Bio import SeqIO
import Bio
from typing import List, Dict, Tuple
from pycircos import color_dict

matplotlib.rcParams['ps.fonttype']       = 42
matplotlib.rcParams['pdf.fonttype']      = 42
matplotlib.rcParams['font.sans-serif']   = ["Arial","Lucida Sans","DejaVu Sans","Lucida Grande","Verdana"]
matplotlib.rcParams['font.family']       = 'sans-serif'
matplotlib.rcParams['font.size']         = 10.0
matplotlib.rcParams["axes.labelcolor"]   = "#000000"
matplotlib.rcParams["axes.linewidth"]    = 1.0
matplotlib.rcParams["xtick.major.width"] = 1.0
matplotlib.rcParams["ytick.major.width"] = 1.0
matplotlib.rcParams['xtick.major.pad']   = 6
matplotlib.rcParams['ytick.major.pad']   = 6
matplotlib.rcParams['xtick.major.size']  = 6
matplotlib.rcParams['ytick.major.size']  = 6

class Garc:
    
    # setting a class attribute for all the colors that can be used
    color_dict: Dict[str, str] = {
        "light red": "#ff8a80",
        "light pink": "#ff80ab",
        "soft magenta":"#ea80fc",
        "light violet": "#b388ff",
        "light blue": "#82b1ff",
        "light cyan": "#84ffff",
        "pale cyan": "#a7ffeb",
        "lime green": "#b9f6ca",
        "light green": "#ccff90",
        "light yellow": "#ffe57f",
        "light orange": "#ffd180",
        "grayish red": "#bcaaa4",
        "light gray": "#eeeeee",
        "grayish blue": "#b0bec5",
        "red": "#ff5252",
        "neon pink": "#ff4081",
        "bright magenta": "#e040fb",
        "purple": "#7c4dff",
        "soft blue": "#536dfe",
        "blue": "#448aff",
        "cyan": "#18ffff",
        "aquamarine": "#64ffda",
        "green yellow": "#b2ff59",
        "yellow": "#ffff00",
        "sunglow": "#ffcc33",
        "sandy brown": "#f4a460",
        "cinereous": "#98817b",
        "gainsboro": "#e0e0e0",
        "gray": "#90a4ae"
    }
    # setting a class attribute for an id
    _arcnum = 0

    @classmethod
    def display_colors(cls) -> List:
        """Class method that shows the possible color options that the user could choose from.
        Returns
        _______
        List[str]
            returns a list of strings that are the keys to the color_dict class attribute that 
            describe what colors can be chosen
        """

        return list(color_dict.keys())

    def __setitem__(self, key, item):
        self.__dict__[key] = item

    def __getitem__(self, key):
        return self.__dict__[key] 
    
    def change_face_color(self, color: str) -> None:
        """Function to change the face color for a specific Garc 
        Parameters
        __________
        color : str
            string that list the color that the user wants to change the facecolor to
        """

        error_message: str = f"color, {color}, not found. Please check the acceptable colors using the class method .display_colors()"

        self.facecolor = Garc.color_dict.get(color, sys.exit(error_message))

    def change_edge_color(self, color: str) -> None:
        """Function to change the edgecolor for a specific Garc 
        Parameters
        __________
        color : str
            string that list the color that the user wants to change the edgecolor to
        """

        error_message: str = f"color, {color}, not found. Please check the acceptable colors using the class method .display_colors()"

        self.edgecolor = Garc.color_dict.get(color, sys.exit(error_message))

    def show_attribute(self) -> None:
        """Function that will show all the attributes of a certain object"""
        print("The attributes of the class are:\n")

        # print out each attribute of the class by accessing the attribute dictionary
        for key in self.__dict__.keys():
            print(key)

    def random_facecolor(self,) -> None:    
        """Function for if the user wants to use random colors for the facecolor"""
        self.facecolor = Garc.color_dist[Garc._arcnum % len(list(Garc.color_dist.keys()))]

    def __init__(self, arc_id: str=None, record: str=None, size: int =1000, interspace: int =3): 

        # setting an id using either the id passed by the user or the class attribute for the id
        if arc_id == None:
            self.arc_id = str(Garc._arcnum) 
        else:
            self.arc_id = arc_id

        # establishing initial attributes of the class
        self._parental_gcircle = None
        self.interspace  = 2 * np.pi * (interspace / 360)
        self.raxis_range: Tuple[int, int] = (500, 550)
        # The face color defaults to grey
        self.facecolor: str = Garc.color_dict["gray"]
        # setting an initial edge color
        self.edgecolor: str = "#303030"
        self.linewidth: float = 0.75

        #The label will initial be set to the arc_id but the user can change this by accessing the attribute
        self.label: str = self.arc_id
        self.labelposition: int = 0
        self.labelsize: int = 10
        self.label_visible: bool = False

        if record is None:
            self.record = None
            self.size: int = 1000
        
        elif type(record) == Bio.SeqRecord.SeqRecord:
            self.record = record
            self.size   = len(str(self.record.seq))
        
        elif type(record) == str:
            match = re.fullmatch("[a-zA-Z]{1,2}_?[0-9]{5,6}", record)
            if os.path.exists(record) == True:
                self.record = SeqIO.read(value, format="genbank")  
            
            if match is None:
                raise ValueError("Incorrect value for NCBI accession number.") 
            else:
                url = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=gbwithparts&id={}&withparts=on".format(record) 
            outb = io.BytesIO()
            outs = io.StringIO()
            headers = {"User-Agent": "Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:47.0) Gecko/20100101 Firefox/47.0"} 
            request = urllib.request.Request(url, headers=headers) 
            
            with urllib.request.urlopen(request) as u:
                outb.write(u.read())
            outs.write(outb.getvalue().decode())
            
            with tempfile.TemporaryFile(mode="w+") as o:
                content = outs.getvalue()
                o.write(content)
                o.seek(0)  
                record = SeqIO.parse(o,"genbank")
                record = next(record)
            self.record = record 
            self.size = len(str(self.record.seq))
        else:
            self.record = None
            self.size = size
        
        Garc._arcnum += 1
    
    def calc_density(self, positions, window_size=1000):
        densities = [] 
        positions.sort()
        for i in range(0, self.size, window_size): 
            source = tuple(range(i, i+window_size))
            amount = 0 
            for pos in positions:
                if type(pos) == int:
                    if pos in source:
                        amount += 1 
                elif type(pos) == tuple:
                    if pos[0] <= source[-1] and pos[1] >= source[0]:
                        amount += 1
                else:
                    raise ValueError("List elements should be int type or tuple consiting of two int values")
            densities.append(amount) 

        source = tuple(range(i,self.size))
        amount = 0 
        for pos in positions:
            if type(pos) == int:
                if pos in source:
                    amount += 1 
            elif type(pos) == tuple:
                if pos[0] <= source[-1] and pos[1] >= source[0]:
                    amount += 1
            else:
                raise ValueError("List elements should be int type or tuple consiting of two int values")
        densities.append(amount*((self.size-i)/window_size))     
        return densities 

    def calc_nnratio(self, n1="G", n2="C", window_size=1000, step_size=None):
        if self.record is None:
            raise ValueError("self.record is None, please specify record value")
        
        if step_size is None:
            step_size = window_size
        
        seq = str(self.record.seq)
        gc_amounts = []
        for i in range(0, len(seq), step_size):
            if n2 is None:
                gc_amount = seq[i:i+window_size].upper().count(n1) * 1.0 / window_size
            else:
                gc_amount = (seq[i:i+window_size].upper().count(n1) + seq[i:i+window_size].upper().count(n2)) * 1.0 / window_size
            gc_amounts.append(gc_amount)
        if n2 is None:
            gc_amounts.append(seq[i:].upper().count(n1) * 1.0 / (len(seq)-i))
        else:
            gc_amounts.append((seq[i:].upper().count(n1) + seq[i:i+window_size].upper().count(n2)) * 1.0 / (len(seq)-i))
        
        self["{}{}_ratio".foramt(n1,n2)] = gc_amounts
        gc_amounts = np.array(gc_amounts)
        return gc_amounts

    def calc_nnskew(self, n1="G", n2="C", window_size=1000, step_size=None):
        #(G-C)/(G+C) 
        if self.record is None:
            raise ValueError("self.record is None, please specify record value")
        
        if step_size is None:
            step_size = window_size
        
        seq = str(self.record.seq) 
        gc_skews = []
        for i in range(0, len(seq), step_size):
            gc_skew = (seq[i:i+window_size].upper().count(n1) - seq[i:i+window_size].upper().count(n2)) * 1.0 / (seq[i:i+window_size].upper().count(n1) + seq[i:i+window_size].upper().count(n2)) * 1.0
            gc_skews.append(gc_skew)
        
        gc_skews.append((seq[i:].upper().count(n1) - seq[i:].upper().count(n2)) * 1.0 / (seq[i:].upper().count(n1) + seq[i:].upper().count(n2)) * 1.0)
        self["{}{}_skew".format(n1,n2)] = gc_skews
        gc_skews = np.array(gc_skews)
        return gc_skews 
  
