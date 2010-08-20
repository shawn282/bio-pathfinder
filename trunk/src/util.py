#!/usr/bin/python

# util.py contains a set of functions for dealing with the file-system, and I/O of text files

import os
import re
import bag
import random
import math
import sys

#####################################################################################
## File System helper files
#####################################################################################

if (sys.argv[0] != '' and os.path.dirname(sys.argv[0]) != ''):
    os.chdir(os.path.dirname(sys.argv[0]))

def get_progdir():
    (progdir, filename) = os.path.split(sys.argv[0])
    if (progdir == ""):
        return "."
    else:
        return progdir

def generate_temp_filename(extension=".tmp"):
    return "%d%s" % (random.randint(0, 99999), extension)

def _mkdir(newdir):
    """works the way a good mkdir should :)
        - already exists, silently complete
        - regular file in the way, raise an exception
        - parent directory(ies) does not exist, make them as well
    """
    if os.path.isdir(newdir):
        pass
    elif os.path.isfile(newdir):
        raise OSError("a file with the same name as the desired " \
                      "dir, '%s', already exists." % newdir)
    else:
        head, tail = os.path.split(newdir)
        if head and not os.path.isdir(head):
            _mkdir(head)
        #print "_mkdir %s" % repr(newdir)
        if tail:
            os.mkdir(newdir)

def embed_link(html_file, html_file_path, relative_path, text):
    html_file.write("<a href=\"" + relative_path + ".html\">" + text + "</a>")
    html_file.flush()
    return open(html_file_path + "/" + relative_path + ".html", "w")

def text2hex(s):
    h = ""
    for i in range(len(s)):
        c = ord(s[i])
        h = h + ("%x" % c)
    return h

def parse_text_file(filename):
    f = open(filename, "r")
    parsed_lines = []
    for line in f.readlines():
        parsed_line = line.rstrip().lstrip()
        if (parsed_line != "" and parsed_line[0] != '#'):
            parsed_lines.append(parsed_line)
    f.close()
    return parsed_lines
 
def write_text_file(lines, filename):
    f = open(filename, "w")
    f.write("\n".join(lines))
    f.close()

def parse_dat(datfile, key_field_name, value_field_name):
    print "Reading DAT file '%s' ... " % datfile
    datfile = open(datfile, "r")
    dict = {}
    key = ""
    values = []
    
    for line in datfile:
        line = line.strip('\n')

        if (line == '//'):
            # a new section starts here, add the key-value pair to the table, and reset them
            if (key != ""):
                dict[key] = values
            key = ""
            values = []
        elif (re.search("^/", line)):
            pass # this line is an extension of the previous line, I am ignoring it for now, since it happens mostly in comments
        elif (re.search("^#", line)):
            pass # this is a remark line, do nothing
        else:
            [field, field_value] = re.split(" - ", line, 1)
            if (field == key_field_name):
                key = field_value
            elif (field == value_field_name):
                values.append(field_value)
    return dict

def readlines_until_mark(file, mark="//"):
    """returns a list containing the lines of the file (without the line-feeds)
       stops when reaching the end of the file, or the given mark.
       Note: doesn't include the mark in the returned list
    """
    lines = []
    while True:
        line = file.readline()
        if (line == "" and lines == []): # reached eof
            return None
        stripped_line = line.rstrip()
        if (line == "" or stripped_line == mark):
            return lines
        lines.append(stripped_line)

def read_next_dat_section(file):
    lines = readlines_until_mark(file, "//")
    if (lines == None):
        return None
    
    dict = {}
    for line in lines:
        if (line[0] == '/'):
            pass # this line is an extension of the previous line, I am ignoring it for now, since it happens mostly in comments
        elif (line[0] == '#'):
            pass # this is a remark line, do nothing
        else:
            try:
                [field, field_value] = line.split(" - ", 1)
            except ValueError, msg:
                raise Exception("syntax error in DAT file: " + line)
            dict[field] = field_value
    return dict

def find_executable(progname, fname_list, critical=False):
    """ Given a list of possible locations, finds the executable file and returns
        that location. If the file is not found, raises an exception (if critical=True)
        or just continues without it.
    """
    for fname in fname_list:
        if (os.path.isfile(fname)):
            sys.stderr.write("Mapping " + progname + " to " + fname + "\n")
            return fname
    if (critical):
        raise Exception("cannot locate the executable for: " + progname)
    else:
        sys.stderr.write("Cannot locate the program: " + progname + ", trying to manage without it.\n")
        return "echo"

#####################################################################################
## User Interface
#####################################################################################

def menu(items):
    sys.stderr.flush()
    for i in range(len(items)):
        print "%2d) %s" % (i+1, items[i])
    
    while (True):
        print "Choose a number between 1 and %d: " % len(items),
        try:
            x = sys.stdin.readline()
            index = int(x)-1
            return items[index]
        except ValueError:
            print "You've entered an invalid input, try again"
            
    
#####################################################################################
## General purpose geometry methods
#####################################################################################
def str2floatvector(s):
    """ Parses a string of comma separated numbers into a tuple of floats
    """
    return tuple([float(x) for x in s.split()])

def str2intvector(s):
    """ Parses a string of comma separated numbers into a tuple of integers
    """
    return tuple([int(x) for x in s.split()])

def abs(v):
    new_v = []
    for i in range(len(v)):
        if (v[i] < 0):
            new_v.append(-v[i])
        else:
            new_v.append(v[i])
    return tuple(new_v)

def pow(v, exponent):
    return tuple([x**exponent for x in v])

def norm(v, exponent=2):
    return sum(pow(abs(v), exponent)) ** (1.0/exponent)    

def subtract(v1, v2):
    """returns the vector (v1 - v2)
    """
    if (len(v1) != len(v2)):
        raise Exception("vector lengths do not match")
    new_v = []
    for i in range(len(v1)):
        new_v.append(v1[i] - v2[i])
    return tuple(new_v)

def add_vector(v1, v2):
    """returns the vector (v1 + v2)
    """
    if (len(v1) != len(v2)):
        raise Exception("vector lengths do not match")
    new_v = []
    for i in range(len(v1)):
        new_v.append(v1[i] + v2[i])
    return tuple(new_v)

def average_vector_list(v_list):
    if (len(v_list) == 0):
        raise Exception("cannot average a vector list of size 0")
    v = v_list[0]
    for i in range(1, len(v_list)):
        v = add_vector(v, v_list[i])
    return scale_vector(v, 1.0 / len(v_list))
    
def distance(pos1, pos2, exponent=2):
    return norm(subtract(pos1, pos2), exponent)

def direction(start, end):
    length = distance(start, end)
        
    new_v = []
    for i in range(len(start)):
        if (length > 0):
            new_v.append((end[i]-start[i])/length)
        else:
            new_v.append(0)
    return tuple(new_v)

def get_2D_angle(start, end):
    """Returns the angle (in Radians) of the vector between start and end, relative to the (1,0) vector
    """
    x = end[0] - start[0]
    y = end[1] - start[1]
    
    if (x == 0):
        if (y >= 0):
            angle = math.pi/2
        else:
            angle = math.pi*3/2
    else:
        angle = math.atan(float(y) / float(x))
        if (x < 0):
            angle += math.pi
    
    return angle % (2 * math.pi)

def radial_2D_direction(start, end):
    return (distance(start, end), get_2D_angle(start, end))

def lin2rad_2D(v):
    return radial_2D_direction((0,0), v)

def rad2lin_2D(radius, alpha):
    return (radius * math.cos(alpha), radius * math.sin(alpha))

def rotate_2D_vector(v, alpha):
    """Rotates a 2D vector by an angle alpha (in Radians)
    """
    if (len(v) != 2):
        raise Exception("This method only accepts 2D vectors")
    
    new_x = v[0] * math.cos(alpha) - v[1] * math.sin(alpha)
    new_y = v[0] * math.sin(alpha) + v[1] * math.cos(alpha)
    return (new_x, new_y)

def scale_vector(v, scalar):
    new_v = []
    for i in range(len(v)):
        new_v.append(v[i] * scalar)
    return tuple(new_v)

def affine_2D_transform(v, offset=(0.0,0.0), scalar=1.0, alpha=0.0):
    """Find the affine transform that maps:
       (0,0) -> offset
       (1,0) -> offset + (scalar*cos(alpha), scalar*sin(alpha))
    """
    if (len(v) != 2):
        raise Exception("This method only accepts 2D vectors")
    new_x = scalar * (v[0] * math.cos(alpha) - v[1] * math.sin(alpha)) + offset[0]
    new_y = scalar * (v[0] * math.sin(alpha) + v[1] * math.cos(alpha)) + offset[1]
    return (new_x, new_y)

def get_affine_transform(new_zero, new_one):
    offset = new_zero
    scalar = distance(new_zero, new_one)
    alpha  = get_2D_angle(new_zero, new_one)
    return (offset, scalar, alpha)

def min_vector(v_list):
    if (v_list == None):
        return None
    dim = len(v_list[0])
    ans = []
    for d in range(dim):
        ans.append(min([v[d] for v in v_list]))
    return tuple(ans)

def max_vector(v_list):
    if (v_list == None):
        return None
    dim = len(v_list[0])
    ans = []
    for d in range(dim):
        ans.append(max([v[d] for v in v_list]))
    return tuple(ans)

def test():
    print str2intvector('1 2 3 4   5   6   ')
    
    for i in range(-20, 20):
        alpha = 2*math.pi*i/10
        print alpha%(2*math.pi) - get_2D_angle((0,0), (math.cos(alpha), math.sin(alpha)))
    
    return
    
if __name__ == '__main__': test()
