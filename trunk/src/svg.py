#!/usr/bin/python

"""\
SVG.py - Construct/display SVG scenes.

The following code is a lightweight wrapper around SVG files. The metaphor
is to construct a scene, add objects to it, and then write it to a file
to display it.

This program uses ImageMagick to display the SVG files. ImageMagick also 
does a remarkable job of converting SVG files into other formats.
"""

import os
import math
import sys
from subprocess import Popen
import util
import platform

# Command to execute to display images.
display_prog = None
for file in ['c:\Program Files\Inkscape\inkview.exe', \
             '/Applications/svg Detective.app/Contents/MacOS/svg Detective', \
             '/Applications/Safari.app/Contents/MacOS/Safari', \
             '/usr/bin/inkview']:
    if (os.path.isfile(file)):
        display_prog = file
        print >> sys.stderr, "Using SVG viewer: %s" % display_prog
        break
if (display_prog == None):
    print >> sys.stderr, "Cannot find a SVG viewer, will not display SVG images"    

class SvgObject:
    def __init__(self):
        self.type = None
        self.attributes = {}
        self.xml_value = ""
        self.min_corner = (0,0)
        self.max_corner = (0,0)
        
    def set_attribute(self, attrib, value):
        if (value != None):
            self.attributes[attrib] = value
    
    def get_attribute(self, attrib):
        return self.attributes.get(attrib, None)
    
    def colorstr_hex(self, rgb): return "#%x%x%x" % (rgb[0]/16, rgb[1]/16, rgb[2]/16)
    
    def set_attribute_color(self, attrib, rgb):
        if (rgb != None):
            self.set_attribute(attrib, "#%x%x%x" % (rgb[0]/16, rgb[1]/16, rgb[2]/16))
    
    def pan(self, offset):
        self.min_corner = util.add_vector(self.min_corner, offset)
        self.max_corner = util.add_vector(self.max_corner, offset)
        pass
    
    def to_string(self):
        line = "\t<" + self.type
        for (attrib, value) in self.attributes.iteritems():
            if (value != None):
                line += " " + attrib + "=\"" + str(value) + "\""
        line += ">" + self.xml_value + "</" + self.type + ">" 
        return line
    
    def __repr__(self):
        return self.to_string()
    
    def clone(self):
        from copy import deepcopy
        return deepcopy(self)
    
    def get_min_corner(self):
        return self.min_corner

    def get_max_corner(self):
        return self.max_corner

class SvgCollection(SvgObject):
    def __init__(self):
        self.type = "collection"
        self.attributes = {}
        self.xml_value = ""
        self.items = []
        return
        
    def set_attribute(self, attrib, value):
        for item in self.items:
            item.set_attribute(attrib, value)
    
    def get_attribute(self, attrib):
        raise Exception("Cannot use get_attribute on an SvgCollection")
    
    def pan(self, offset):
        for item in self.items:
            item.pan(offset)
            
    def get_min_corner(self):
        return util.min_vector([item.get_min_corner() for item in self.items])

    def get_max_corner(self):
        return util.max_vector([item.get_max_corner() for item in self.items])

    def to_string(self):
        line = ""
        for item in self.items:
            line += item.to_string() + "\n"
        return line
    
    def add(self, item):
        if (isinstance(item, SvgObject)):
            self.items.append(item)
        else:
            raise Exception("Can only add objects of type SvgObject to a SvgCollection")
        
    def __iter__(self):
        return iter(self.items)


global default_font_size; default_font_size=12
class Scene(SvgObject):
    def __init__(self, width=400, height=400, font_size=default_font_size):
        SvgObject.__init__(self)
        self.type = "svg"
        self.svgname = None
        self.draw_border = False

        self.set_attribute("xmlns:xlink", "http://www.w3.org/1999/xlink")
        self.set_attribute("xmlns", "http://www.w3.org/2000/svg")
        self.set_attribute("color-rendering", "auto")
        self.set_attribute("color-interpolation", "auto")
        self.set_attribute("text-rendering", "auto")
        self.set_attribute("image-rendering", "auto")
        self.set_attribute("shape-rendering", "auto")
        
        self.set_attribute("width", width)
        self.set_attribute("height", height)

        self.set_attribute("fill-opacity", 1)
        self.set_attribute_color("fill", (0,0,0))
        self.set_attribute_color("stroke", (0,0,0))
        self.set_attribute("stroke-width", 1)
        self.set_attribute("stroke-opacity", 1)
        self.set_attribute("stroke-linecap", "square")
        self.set_attribute("stroke-miterlimit", 10)
        self.set_attribute("stroke-linejoin", "miter")
        self.set_attribute("stroke-dashoffset", 0)
        self.set_attribute("stroke-dasharray", "none")

        self.set_attribute("font-size", font_size)
        self.set_attribute("font-weight", "normal")
        self.set_attribute("font-family", "monospace")
        self.set_attribute("font-style", "normal")

        # add a black square frame around the scene
        self.items = SvgCollection()
        #self.items.append(Rectangle((0,0), width, height, (255,255,255), 1, (0,0,0), 1))
        return

    def width(self):
        return self.get_attribute("width")

    def height(self):
        return self.get_attribute("height")
    
    def font_size(self):
        return self.get_attribute("font-size")

    def add(self, other, offset=(0,0)):
        """Adds one item (or appends another scene) to the scene
           It is possible to add an offset to the new items
        """
        if (isinstance(other, Scene)):
            for item in other.items:
                newitem = item.clone()
                newitem.pan(offset)
                self.items.add(newitem)
        elif (isinstance(other, SvgObject)):
            newitem = other.clone()
            newitem.pan(offset)
            self.items.add(newitem)

    def pan(self, offset):
        for item in self.items:
            item.pan(offset)
            
    def justify(self, margin=20):
        (min_x, min_y) = self.items.get_min_corner()
        self.pan((margin - min_x, margin - min_y))

        (width, height) = self.items.get_max_corner()
        self.set_attribute("width", int(width + 2*margin))
        self.set_attribute("height", int(height + 2*margin))

    def border(self, flag=True):
        self.draw_border = flag

    def to_string(self):
        line = "<?xml version=\"1.0\"?>\n"
        line += "\n"
        line += "<!DOCTYPE svg PUBLIC '-//W3C//DTD SVG 1.0//EN' 'http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd'>\n"
        self.xml_value = "\n\t<g>"
        for item in self.items:
            self.xml_value += "\n" + str(item)
        
        if (self.draw_border):
            self.xml_value += "\n" + str(Rectangle((1,1), self.width()-2, self.height()-2, fill_opacity=0, stroke_opacity=1))
        
        self.xml_value += "\n\t</g>"
        line += SvgObject.to_string(self)
        
        return line

    def write_svg(self, filename="untitled.svg"):
        file = open(filename,'w')
        file.write(str(self))
        file.close()
        self.svgname = filename
        return

    def embed_in_html(self, html_file, html_file_path, relative_path):
        html_file.write("<object type=\"image/svg+xml\" name=\"omap\" width=\"%d\" height=\"%d\" " % (self.width(), self.height()));
        html_file.write("data=\"" + relative_path + ".svg\"></object>\n")
        html_file.flush()
        self.write_svg(html_file_path + "/" + relative_path + ".svg")

    def display(self):
        if (self.svgname == None):
            display_filename = util.generate_temp_filename(".svg")
        else:
            display_filename = self.svgname

        file = open(display_filename, 'w')
        file.write(str(self))
        file.close()
            
        p = Popen([display_prog, display_filename])
        os.waitpid(p.pid, 0)
        os.remove(display_filename)
    
class Line(SvgObject):
    def __init__(self, start, end, stroke_color=(0,0,0), stroke_width=1, stroke_opacity=1):
        SvgObject.__init__(self)
        self.type = "line"
        
        self.set_attribute("x1", start[0])
        self.set_attribute("y1", start[1])
        self.set_attribute("x2", end[0])
        self.set_attribute("y2", end[1])
        self.set_attribute_color("stroke", stroke_color)
        self.set_attribute("stroke-width", stroke_width)
        self.set_attribute("stroke-opacity", stroke_opacity)
        
        self.min_corner = util.min_vector([start, end])
        self.max_corner = util.max_vector([start, end])

    def pan(self, offset):
        SvgObject.pan(self, offset)
        self.set_attribute("x1", self.get_attribute("x1") + offset[0])
        self.set_attribute("y1", self.get_attribute("y1") + offset[1])
        self.set_attribute("x2", self.get_attribute("x2") + offset[0])
        self.set_attribute("y2", self.get_attribute("y2") + offset[1])

class Arrow(SvgCollection):
    def __init__(self, start, end, arrow_size=10, stroke_color=(0,0,0), stroke_width=1, stroke_opacity=1):
        SvgCollection.__init__(self)
        self.type = "arrow"
        self.add(Line(start, end, stroke_color, stroke_width, stroke_opacity))
        
        arrow_direction = util.direction(start, end)
        d1 = util.scale_vector(util.rotate_2D_vector(arrow_direction, math.pi *  5/6), arrow_size)
        d2 = util.scale_vector(util.rotate_2D_vector(arrow_direction, math.pi * -5/6), arrow_size)
        e1 = util.add_vector(end, d1)
        e2 = util.add_vector(end, d2)
        self.add(Line(end, e1, stroke_color, stroke_width, stroke_opacity))
        self.add(Line(end, e2, stroke_color, stroke_width, stroke_opacity))

class ChemicalArrow(SvgCollection):
    def __init__(self, start, end, forward=True, backward=True, width=10, stroke_color=(0,0,0), stroke_width=1, stroke_opacity=1):
        SvgCollection.__init__(self)
        self.type = "chemical_arrow"

        (offset, scalar, alpha) = util.get_affine_transform(start, end)
        points = [(0, -0.1),(1, -0.1),(0.75, -0.3),(1, 0.1),(0, 0.1),(0.25, 0.3)]
        points = [util.affine_2D_transform(v, offset, scalar, alpha) for v in points]
        if (forward):
            self.add(Line(points[0], points[1], stroke_color, stroke_width, stroke_opacity))
            self.add(Line(points[1], points[2], stroke_color, stroke_width, stroke_opacity))
        if (backward):
            self.add(Line(points[3], points[4], stroke_color, stroke_width, stroke_opacity))
            self.add(Line(points[4], points[5], stroke_color, stroke_width, stroke_opacity))

class Circle(SvgObject):
    def __init__(self, center, radius, fill_color=(255,255,255), fill_opacity=1, stroke_color=(0,0,0), stroke_width=1, stroke_opacity=1):
        SvgObject.__init__(self)
        self.type = "circle"
        
        self.set_attribute("cx", center[0])
        self.set_attribute("cy", center[1])
        self.set_attribute("r", radius)
        self.set_attribute_color("fill", fill_color)
        self.set_attribute("fill-opacity", fill_opacity)
        self.set_attribute_color("stroke", stroke_color)
        self.set_attribute("stroke-width", stroke_width)
        self.set_attribute("stroke-opacity", stroke_opacity)

        self.min_corner = (center[0] - radius, center[1] - radius)
        self.max_corner = (center[0] + radius, center[1] + radius)

    def pan(self, offset):
        SvgObject.pan(self, offset)
        self.set_attribute("cx", self.get_attribute("cx") + offset[0])
        self.set_attribute("cy", self.get_attribute("cy") + offset[1])

class Rectangle(SvgObject):
    def __init__(self, origin, width, height, fill_color=(255,255,255), fill_opacity=1, stroke_color=(0,0,0), stroke_width=1, stroke_opacity=0):
        SvgObject.__init__(self)
        self.type = "rect"
        
        self.set_attribute("x", origin[0])
        self.set_attribute("y", origin[1])
        self.set_attribute("height", height)
        self.set_attribute("width", width)
        self.set_attribute_color("fill", fill_color)
        self.set_attribute("fill-opacity", fill_opacity)
        self.set_attribute_color("stroke", stroke_color)
        self.set_attribute("stroke-width", stroke_width)
        self.set_attribute("stroke-opacity", stroke_opacity)
        
        self.min_corner = origin
        self.max_corner = (origin[0] + width, origin[1] + height)

    def pan(self, offset):
        SvgObject.pan(self, offset)
        self.set_attribute("x", self.get_attribute("x") + offset[0])
        self.set_attribute("y", self.get_attribute("y") + offset[1])

class Text(SvgObject):
    def __init__(self, origin, text, font_size=None, fill_color=None, stroke_color=None, stroke_width=0):
        SvgObject.__init__(self)
        self.type = "text"
        self.xml_value = text

        self.set_attribute("x", origin[0])
        self.set_attribute("y", origin[1])
        self.set_attribute_color("fill", fill_color)
        self.set_attribute("font-size", font_size)
        self.set_attribute_color("stroke", stroke_color)
        self.set_attribute("stroke-width", stroke_width)

        if (font_size == None):
            font_size = default_font_size
        self.min_corner = util.add_vector(origin, (0, -font_size))
        self.max_corner = util.add_vector(origin, (0, 0)) # I don't know how to asses the pixel length of the text

    def pan(self, offset):
        SvgObject.pan(self, offset)
        self.set_attribute("x", self.get_attribute("x") + offset[0])
        self.set_attribute("y", self.get_attribute("y") + offset[1])
        
def test():
    scene = Scene()
    scene.font_size = 24
    scene.add(Rectangle((100,100),200,200,(0,255,255),1,(255,0,0),30))
    scene.add(Line((200,200),(200,300),(0,0,0),5))
    scene.add(Line((200,200),(300,200),(0,0,0),5))
    scene.add(Line((200,200),(100,200),(0,0,0),5))
    scene.add(Arrow((200,200),(200,150),stroke_color=(0,0,0),stroke_width=2))
    scene.add(Circle((200,200),30,(0,0,255),0.1,(0,255,0),3))
    scene.add(Circle((200,300),30,(0,255,0),0.2,(255,0,0),5))
    scene.add(Circle((300,200),30,(255,0,0),0.3,(0,0,255),2))
    scene.add(Circle((100,200),30,(255,255,0),0.5,(0,0,0),1))
    scene.add(Circle((200,100),30,(255,0,255),1,(0,0,0),7))
    #scene.add(Text((50,50),"Testing SVG", 24, (0,0,0), (255,0,0), 2))
    scene.add(Text((50,50),"Testing SVG"))
    scene.add(ChemicalArrow((200,200),(210,140),stroke_color=(0,0,0),stroke_width=2))
    scene.justify()
    scene.border()
    scene.display()
    print >> sys.stderr, "Finished!"
    return

if __name__ == '__main__': test()
