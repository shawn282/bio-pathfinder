#!/usr/bin/python

import os
import sys
#from chemconvert import mol2svg
import svg
from chemistry import ChemGraph
from Tkinter import Message, Button, Tk
from tkFileDialog import askopenfilename
from chemconstants import *
from chemconvert import molfile2graph

def main():
    root = Tk()
    if (len(sys.argv) > 1):
        mol_filename = sys.argv[1]
    else:
        #print >> sys.stderr, "Syntax: %s <mol file> [<svg file>]" % sys.argv[0]
        mol_filename = askopenfilename(parent=root, title="Pick a MOL file to plot", filetypes=[("MOL files", "*.mol"),("All Files", "*")])
        if (mol_filename == ''):
            sys.exit(-1)

    if (len(sys.argv) > 2):
        svg_filename = sys.argv[2]
    else:
        svg_filename = None
        
    if (not os.path.isfile(mol_filename)):
        print >> sys.stderr, "File not found: " + mol_filename
        sys.exit(-2)
    
    #file = open(mol_filename)
    try:
        G = molfile2graph(mol_filename)
    except Exception, strerror:
        root.title("ERROR")
        Message(root, text=strerror, width=100).pack()
        Button(root, text="OK", command=root.destroy).pack()
        root.mainloop()
        sys.exit(-1)

    scene = svg.Scene(400, 400, 20)
    G.svg(scene)
    scene.set_attribute("height", 500)
    h = G.hash()
    scene.add(svg.Line((0, 400), (400, 400), stroke_width=3, stroke_color=black))
    scene.add(svg.Text((20, 450), h, font_size=20))
    print h
    #scene = mol2svg(file.read(), 500, 500, 20)
        
    if (svg_filename != None):
        scene.write_svg(svg_filename)
    else:
        scene.display()

if __name__ == '__main__': main()