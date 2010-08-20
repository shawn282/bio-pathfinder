import util
import sys
import os
import re

inkscape = util.find_executable("inkscape", ["/Applications/Inkscape.app/Contents/Resources/bin/inkscape", "/usr/bin/inkscape"])
def svg2eps(svg_fname, eps_fname):
    os.popen(inkscape + " " + svg_fname + " --export-eps=" + eps_fname)

def svg2png(svg_fname, png_fname, width="500"):
    os.popen(inkscape + " " + svg_fname + " --export-png=" + png_fname + " --export-width=" + width)

def svg2pdf(svg_fname, pdf_fname):
    os.popen(inkscape + " " + svg_fname + " --export-pdf=" + pdf_fname)

## MAIN
if (len(sys.argv) < 2):
    raise Exception("Syntax: batch_svg2image <path> <output format>")

for fname in os.listdir(sys.argv[1]):
    if (fname[len(fname)-4:] == ".svg"):
        base_fname = sys.argv[1] + "/" + fname[0:len(fname)-4]
        print "converting", base_fname
        if (len(sys.argv) < 3 or sys.argv[2] == "png"):
            svg2png(base_fname + ".svg", base_fname + ".png")
        elif (sys.argv[2] == "eps"):
            svg2eps(base_fname + ".svg", base_fname + ".eps")
        elif (sys.argv[2] == "pdf"):
            svg2pdf(base_fname + ".svg", base_fname + ".pdf")
        else:
            raise Exception("unknown output figure format: " + sys.argv[2])