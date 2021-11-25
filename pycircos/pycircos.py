import os 
import re
import io 
import sys
import math
import urllib
import tempfile
import requests
import collections
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.path    as mpath
import matplotlib.patches as mpatches
from Bio import SeqIO
import Bio
from typing import List, Dict, Tuple

matplotlib.rcParams["figure.max_open_warning"] = 0
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
    #list100 = ["#ffcdd2","#f8bbd0","#e1bee7","#d1c4e9","#c5cae9","#bbdefb","#b3e5fc","#b2ebf2","#b2dfdb","#c8e6c9","#dcedc8","#f0f4c3","#fff9c4","#ffecb3","#ffe0b2","#ffccbc","#d7ccc8","#cfd8dc",
    colorlist = ["#ff8a80","#ff80ab","#ea80fc","#b388ff","#8c9eff","#82b1ff","#84ffff","#a7ffeb","#b9f6ca","#ccff90","#f4ff81","#ffff8d","#ffe57f","#ffd180","#ff9e80","#bcaaa4","#eeeeee","#b0bec5",
                 "#ff5252","#ff4081","#e040fb","#7c4dff","#536dfe","#448aff","#18ffff","#64ffda","#69f0ae","#b2ff59","#eeff41","#ffff00","#ffd740","#ffab40","#ff6e40","#a1887f","#e0e0e0","#90a4ae"]
    _arcnum = 0
    def __setitem__(self, key, item):
        self.__dict__[key] = item

    def __getitem__(self, key):
        return self.__dict__[key] 

    def __init__(self, arc_id=None, record=None, size=1000, interspace=3, raxis_range=(500, 550), facecolor=None, edgecolor="#303030", linewidth=0.75, label=None, labelposition=0, labelsize=10, label_visible=False): 
        self._parental_gcircle = None
        if arc_id == None:
            self.arc_id = str(Garc._arcnum) 
        else:
            self.arc_id = arc_id

        if record is None:
            self.record = None
            self.size = size
        
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
        
        if facecolor is None:
            facecolor = Garc.colorlist[Garc._arcnum % len(Garc.colorlist)] 
        self.interspace  = 2 * np.pi * (interspace / 360)
        self.raxis_range = raxis_range 
        self.facecolor   = facecolor 
        self.edgecolor   = edgecolor
        self.linewidth   = linewidth
        
        if label is None:
            self.label = arc_id
        else:
            self.label = label

        self.label_visible = label_visible
        self.labelposition = labelposition
        self.labelsize = labelsize
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
  
