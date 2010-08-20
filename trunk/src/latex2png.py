#!/usr/bin/python

from urllib import urlopen
import sys

sys.stdout.write("Enter a LaTeX equation:\n")
latex_str = sys.stdin.readline()

proxies = {"http" : "http://wwwproxy.weizmann.ac.il:8080"}
url = "http://www.sitmo.com/gg/latex/latex2png.2.php?z=100&eq=" + latex_str
outfile = "/home/eladn/Desktop/latex.png"

instream = urlopen(url, None, proxies)
file = open(outfile, "w")
file.write(instream.read())
file.close()
