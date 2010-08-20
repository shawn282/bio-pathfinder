#!/usr/bin/python

from chemconvert import smiles2graph, hash2graph
from chemistry import ChemException
import sys
import math
import os
import svg
import Tkinter
from tkFileDialog import asksaveasfile

_graph = None

def read_smiles():
    error_stringvar.set("")
    try:
        smiles = input_stringvar.get()
        graph = smiles2graph(smiles)
        error_stringvar.set("Ready.\nSMILES = %s" % smiles)
        graph_list.append(graph)
    except ChemException, strerror:
        error_stringvar.set("ERROR :" + str(strerror))
        graph = None
    root.update()
    return

def read_hash():
    try:
        hash = input_stringvar.get()
        graph = hash2graph(hash)
        graph.initialize_pos()
        error_stringvar.set("Ready.\nHASH = %s" % hash)
        graph_list.append(graph)
    except ChemException, strerror:
        error_stringvar.set("ERROR :" + str(strerror))
        root.update()
        graph = None
    root.update()
    return

def get_svg():
    if (graph_list == []):
        error_stringvar.set("Not ready yet, please read a new graph")
        root.update()
        return None
    scene = svg.Scene(400, 400, font_size=24)
    graph_list[-1].svg(scene)
    return scene

def plot():
    scene = get_svg()
    if (scene == None):
        return
    scene.display()
    root.destroy()

def save_as():
    scene = get_svg()
    if (scene == None):
        return
    
    svg_filename = asksaveasfile(parent=root, title="Save the graph to an SVG file...")
    if (svg_filename != None):
        svg_filename.write(str(scene))

def main():
    global root
    global error_stringvar
    global input_stringvar
    global graph_list
    
    root = Tkinter.Tk()   # create a root window
    root.title("Enter Hash string:")
    root.setvar("width", 400)
   
    input_stringvar = Tkinter.StringVar()
    input_stringvar.set("")
    
    error_stringvar = Tkinter.StringVar()
    error_stringvar.set("")
    
    graph_list = []
    
    Tkinter.Entry(root, width=50, textvariable=input_stringvar).pack(padx=10, pady=10)
    Tkinter.Button(root, text="Read as HASH", command=lambda: read_hash()).pack(pady=10)
    Tkinter.Button(root, text="Read as SMILES", command=lambda: read_smiles()).pack(pady=10)
    Tkinter.Button(root, text="Plot", command=lambda: plot()).pack(pady=10)
    Tkinter.Button(root, text="Save as ...", command=lambda: save_as()).pack(pady=10)
    Tkinter.Button(root, text="Quit", command=lambda: root.destroy()).pack(pady=10)
    Tkinter.Message(root, width=300, textvariable=error_stringvar).pack(padx=10, pady=10)
    root.mainloop() # create an event loop

if __name__ == '__main__': main()
