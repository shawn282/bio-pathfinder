#!/usr/bin/python

from chemistry import ChemGraph, deduce_permutation
from chemconvert import compound2graph
import chemconstants
import bag
from svg import *
import util
import html_writer
from pathfinder import *

def main():
    util._mkdir("../results/svg")
    #substrate = '_GAP + _GAP'
    #substrate = '_GLC-6P_OCF'
    #(substrate, target) = ('L-Serine', '_3-PG')
    (substrate, target) = ('3-Phospho-Hydroxypyruvate', 'L-Serine')
    G_subs = compounds2graph(substrate)
    h_target = compounds2graph(target).hash()
    h_subs = G_subs.hash()

    product_html_writer = html_writer.HtmlWriter("../results/valid_products.html")
    product_html_writer.write("<center>Substrate : " + strip_hash_chirality(h_subs) + "</center><br>")
    product_html_writer.write("<center>Target : " + strip_hash_chirality(h_target) + "</center><br>")
    product_html_writer.write_svg(G_subs.svg(Scene(400, 200, 10)), "svg/" + h_subs)
    product_html_writer.write("<p>")
    
    for (h_substrate, G_product, reaction) in generate_new_compounds({h_subs : G_subs}, False):
        h_product = G_product.hash()
        product_html_writer.write("<br><center>" + reaction + " : " + strip_hash_chirality(h_product) + "</center><br>")

        if (strip_hash_chirality(h_product) == strip_hash_chirality(h_target)):
            product_html_writer.write("<center>Target Found !!!</center><br>")
        product_html_writer.write_svg(G_product.svg(Scene(400, 200, 10)), "svg/" + h_product)

    product_html_writer.write("</p>")

    product_html_writer.close()
    return
        
if __name__ == '__main__': main()
