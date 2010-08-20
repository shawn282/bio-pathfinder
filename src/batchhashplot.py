#!/usr/bin/python

import sys
import os
import util
from chemconvert import hash2graph
from html_writer import HtmlWriter
from svg import Scene

html = HtmlWriter("../results/hash_list.html")
util._mkdir("../results/hash_list")

for line in util.parse_text_file(sys.argv[1]):
	print line
	graph = hash2graph(line)
	graph.initialize_pos()
	scene = graph.svg(Scene(200, 200, font_size=12))    
	html.write_svg(scene, "../results/hash_list/" + line)

html.display()
