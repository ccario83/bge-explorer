#! /usr/bin/python

######################################################################
#                                                                    #
#           PLOT ARROWS FOR GENE CLUSTER GIVEN A GenBank FILE        #
#                           Peter Cimermancic                        #
#                               April 2010                           #
#                                                                    #
#                         modified by Clint Cario                    #
#                               May 2014                             #
######################################################################

import sys, copy
import numpy as np
from math import floor
import colorsys
import json

from pprint import pprint

# --- Draw arrow for gene
def arrow(X,Y,L,H,strand,h,l,color='#eeeeee',fill_opacity=1.0,stroke_width=1.5,title='undefined',description='undefined'):
    '''
    SVG code for arrow:
        - (X,Y) ... upper left (+) or right (-) corner of the arrow
        - L ... arrow length
        - H ... arrow height
        - strand
        - h ... arrow head edge width
        - l ... arrow head length
        - color
        - strand
    the edges are ABCDEFG starting from (X,Y) like:

      A           B|\
(X,Y) _____________| \ 
     |                \C
    G|_____________   /
      F           E| / 
                   |/
                     D
     the arrow is reversed for '-' strands. If the orientation is unknown:

      
(X,Y)A_________________ B
     |                 |
     |_________________|
    D                   C


    '''
    
    if strand == '+' or strand == 1:
        A = [X,Y]
        B = [X+L-l,Y]
        C = [X+L-l,Y-h]
        D = [X+L,Y+H/2]
        E = [X+L-l,Y+H+h+1]
        F = [X+L-l,Y+H]
        G = [X,Y+H]

        if L < l:
            # squeeze arrow if length shorter than head length
            B = [X,Y]
            C = [X,Y-h]
            D = [X+L,Y+H/2]
            E = [X,Y+H+h+1]
            F = [X,Y+H]

        line = """
  <polygon points="%i,%i %i,%i %i,%i %i,%i %i,%i %i,%i %i,%i" style="fill:%s; fill-opacity:%f; stroke:#000000;stroke-width:%.2f">
    <title>%s</title>
    <desc>%s</desc>
  </polygon> """ % (A[0],A[1],B[0],B[1],C[0],C[1],D[0],D[1],E[0],E[1],F[0],F[1],G[0],G[1],color,fill_opacity,stroke_width,title,description)
    
    elif strand == '-' or strand == -1:
        A = [X+L,Y]
        B = [X+l,Y]
        C = [X+l,Y-h]
        D = [X,Y+H/2]
        E = [X+l,Y+H+h+1]
        F = [X+l,Y+H]
        G = [X+L,Y+H]

        if L < l:
            # squeeze arrow if length shorter than head length
            B = [X+L,Y]
            C = [X+L,Y-h]
            D = [X,Y+H/2]
            E = [X+L,Y+H+h+1]
            F = [X+L,Y+H]

        line = """
  <polygon points="%i,%i %i,%i %i,%i %i,%i %i,%i %i,%i %i,%i" style="fill:%s; fill-opacity:%f; stroke:#000000;stroke-width:%.2f">
    <title>%s</title>
    <desc>%s</desc>
  </polygon>""" % (A[0],A[1],B[0],B[1],C[0],C[1],D[0],D[1],E[0],E[1],F[0],F[1],G[0],G[1],color,fill_opacity,stroke_width,title,description)
    
    else:
        A = [X,Y]
        B = [X+L,Y]
        C = [X+L,Y+H]
        D = [X,Y+H]
        
        line = """
  <polygon points="%i,%i %i,%i %i,%i %i,%i" style="fill:%s; fill-opacity:%f; stroke:#000000;stroke-width:%.2f">
    <title>%s</title>
    <desc>%s</desc>
  </polygon>""" % (A[0],A[1],B[0],B[1],C[0],C[1],D[0],D[1],color,fill_opacity,stroke_width,title,description)

    return line


def rectangle(X,Y,l,h,pfam='text should be here',pfamde='',color='#dddddd'):
    '''
    Draw the rectangle for annotation
    '''
    line = """\n  <rect x="%i" y="%i" width="%i" height="%i" style="fill:%s">
        <title>%s</title><desc>%s</desc>
      </rect>""" % (X,Y,l,h,color,pfam,pfamde)
    return line


def line(X,Y,L):
    '''
    Draw a line below genes
    '''
    line =     """\n  <line x1="%i" y1="%i" x2="%i" y2="%i" style="stroke:rgb(130,130,130);stroke-width:2"/>""" % (X,Y,X+L,Y)
    return line


def SVG(attribs,show_domains,ArrowHeight=20,HeadEdge=8,HeadLength=10,marginX=100,marginY=30,scaling=100.0,gene_colors=None,color_key=None):
    '''
    Create the main SVG document:
        - read in GenBank documnet
        - find genes, start and stop positions, and strands
        - write the SVG files
    '''
    genes = attribs['genes']
    #domains = attribs['domains'] if 'domains' in attribs.keys() else None
    domains = attribs['domains'] if show_domains else None

    # Use empty strings if the following don't exist prevents key errors
    for g in genes:
        not_present = [i for i in ['strand','description','notes','locus'] if i not in g.keys()]
        for i in not_present:
            g[i] = ''

    # --- rotate if too many minuses
    maxs    = max([g['end'] for g in genes])
    strands =     [g['strand'] for g in genes]
    if (strands.count('-')+strands.count(-1)) > (strands.count('+')+strands.count(1)):
        for gene in genes:
            start,end = abs(gene['start']-maxs),abs(gene['end']-maxs)
            gene['start'] = end
            gene['end']   = start
            if (gene['strand']=='-' or gene['strand']==-1): gene['strand']=1
            elif (gene['strand']=='+' or gene['strand']==1): gene['strand']=-1
        # Do the same for the domains if they exist
        if domains != None:
            for domain in domains:
                start,end = abs(domain['gene_start']-maxs),abs(domain['gene_end']-maxs)
                domain['gene_start'] = end
                domain['gene_end']   = start

    # --- create SVG header
    all_text = ''
    # draw a line that corespond to cluster size
    mins,maxs = min([g['start'] for g in genes]),max([g['end'] for g in genes])
    cluster_size = maxs-mins
    all_text += line(marginX,marginY+ArrowHeight/2,cluster_size / scaling)

    if show_domains:
        domains = attribs['domains']
        mins = min([d['gene_start'] for d in domains])
        for domain in domains:
            #print "%10d, %10d, %10d, %10d, %10d, %10d, %10d" % (domain['domain_id'], domain['gene_start'], domain['gene_end'], domain['domain_start'], domain['domain_end'], domain['start'], domain['end'])
            start = (domain['gene_start']+domain['domain_start']-mins)/scaling
            end = (domain['gene_start']+domain['domain_end']-mins)/scaling
            # write arrow to SVG file
            all_text += arrow(start+marginX,marginY,end-start,ArrowHeight,None,HeadEdge,HeadLength,"#000000",1.0,1.0,domain['pfam_id'],"")

    mins = min([g['start'] for g in genes])
    for gene in genes:
        #print gene['gene_id'], gene['start'], gene['end'], gene['description'] 
        strand = gene['strand']
        start = (gene['start']-mins)/scaling
        end = (gene['end']-mins)/scaling
        try:
            color = gene_colors[gene[color_key]]
        except:
            color = "#eeeeee"
        pfams = "-".join([d['pfam_id'] for d in attribs['domains'] if d['gene']==gene['gene_id']])
        popup = str(gene['gene_id']) + ": " + gene['locus'] + ": " + gene['description'] + ": " + pfams
        # write arrow to SVG file
        all_text += arrow(start+marginX,marginY,end-start,ArrowHeight,strand,HeadEdge,HeadLength,color,0.8,1.5,popup,gene['notes'])


    
    return all_text


def generate_palette(N=15):
    HSV_tuples = [(x*1.0/N, 0.85, 0.85) for x in range(N)]
    RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
    RGB_tuples = map(lambda x: tuple(map(lambda y: int(y * 255),x)),RGB_tuples)
    HEX_tuples = map(lambda x: tuple(map(lambda y: chr(y).encode('hex'),x)), RGB_tuples)
    HEX_tuples = map(lambda x: "".join(x), HEX_tuples)
    HEX_tuples = map(lambda x: "#"+x.upper(), HEX_tuples)
    return HEX_tuples


def experimental_colors(bgcs):
    # Get all pfams in each gene
    pairs = []
    for id_,b in bgcs.iteritems():
        if 'domains' in b.keys():
            for d in b['domains']:
                pairs.append((d['gene'],d['pfam_id']))
    # Collapse
    gene_loci = {}
    for pair in pairs:
        id_, pf = pair
        if id_ in gene_loci.keys():
            gene_loci[id_].append(pf)
        else: 
            gene_loci[id_] = [pf]
    # Sort and Join as string (now unique key)
    for k,d in gene_loci.iteritems():
        d.sort()
        gene_loci[k] = "-".join(d)
    # uniquify 
    unq = list(set([d for k,d in gene_loci.iteritems()]))
    colors = generate_palette() if (len(unq) >= 15) else generate_palette(len(unq))
    # Assign a color map for pfam_id => hex color code
    palette = {}
    for idx,id_ in enumerate(unq):
        value = colors[idx%len(colors)]
        palette[id_] = value
    # Make mapping gene_id => color instead of pfam list
    new_palette = {}   
    for k,d in gene_loci.iteritems():
        new_palette[k] = palette[d]
    return new_palette


def assign_colors(bgcs, color_key, default_term, generic_terms):
    # Get a list of unique locus tags from all bgcs/genes
    gene_loci = []
    # id_,b = [(a,b) for a,b in bgcs.iteritems() if a==10545][0]
    for id_,b in bgcs.iteritems():
        # g = b['genes'][0]
        for g in b['genes']:
            # While looping over the genes, may as well load json notes
            try:
                g['notes'] = json.loads(g['notes'])
            except:
                g['notes'] = 'no additional info found in local db'
            # This is whatever you want to color by, make be the same in the SVG function when calling arrow(color=XXX)!
            if color_key not in g.keys():
                g[color_key] = default_term
            gene_loci.append(g[color_key])
    gene_loci = list(set(gene_loci))
    # Define colors
    colors = generate_palette() if (len(gene_loci) >= 15) else generate_palette(len(gene_loci))

    # Assign a color map for pfam_id => hex color code
    palette = {}
    for idx,id_ in enumerate(gene_loci):
        value = colors[idx%len(colors)]
        palette[id_] = value
    for term in generic_terms:
        palette[term] = "#eeeeee"
    return palette


def arrower(bgcs, show_domains=False, color_key='pfams'):
    # what to color by
    palette = None
    if color_key == 'description':
        default_term = 'hypothetical protein'
        generic_terms = ['hypothetical protein','conserved protein domain']
        palette = assign_colors(bgcs, color_key, default_term, generic_terms)
    elif color_key == 'pfams':
        color_key = 'gene_id'
        palette = experimental_colors(bgcs)
    
    # --- SVG string  (took out width="21.59cm" height="27.94cm")
    header = """<?xml version="1.0" standalone="no"?>
    <!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
    <svg id="arrow" height="%d" xmlns='http://www.w3.org/2000/svg' xmlns:xlink='http://www.w3.org/1999/xlink' onload='Init(evt)'>"""

    all_text = "\n  <script xlink:href='ToolTip.js' type='text/ecmascript'><![CDATA[]]></script>"
    H = 20/2.          # arrow width
    E = 3/2.           # arrow head edge size
    l = 10/2.          # arrow head edge length
    X = 20             # left-site margins
    Y = 50             # top-site margins
    S = 75.0           # scaling factor
    F = 14             # font size

    ## bgcs = {'BGC_ID' => {'genome_name', start_locus_tag', 'end_locus_tag', 'domains' => [{'locus_description', 'pfam_id', start', 'end'}]}}
    ## For now, locus_description is unknown and  = "" (empty string)
    # Debug
    # bgc_id, attribs = [(a,b) for a,b in bgcs.iteritems() if a==10545][0]
    for bgc_id, attribs in bgcs.iteritems():
        all_text += "\n\n<g id='arrow_%i'><text x='20' y='%i' font-family='Arial' font-size='12'>%s [%s-%s]</text>" %\
                    (bgc_id, Y/2.-8, attribs['genome_name'], attribs['start_locus_tag'], attribs['end_locus_tag'])
        all_text += SVG(attribs,show_domains,ArrowHeight=H,HeadEdge=E,HeadLength=l,marginX=X,marginY=Y/2.,scaling=S,gene_colors=palette,color_key=color_key)
        all_text += "</g>"
        Y += 90

    all_text += """
      <g id='ToolTip' opacity='0.8' visibility='hidden' pointer-events='none'>
        <rect id='tipbox' x='0' y='5' width='88' height='20' rx='2' ry='2' fill='white' stroke='black' />
        <text id='tipText' x='5' y='20' font-family='Arial' font-size='12'>
          <tspan id='tipTitle' x='5' font-weight='bold'><![CDATA[]]></tspan>
          <tspan id='tipDesc' x='5' dy='15' fill='blue'><![CDATA[]]></tspan>
        </text>
      </g>"""

    all_text = header%int(Y/2) + all_text + "\n</svg>"

    return all_text
    #output = open('APE_%i.xml' % cnt,'w')
    #output.write(all_text)
    #output.close()


'''
Experimental palette generator
def generate_palette(N):
    HSV_tuples = [((x%25)*1.0/25, 1.0-(floor(x/25)*.2), 1.0-(floor(x/25)*.2)) for x in range(N)]
    RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
    RGB_tuples = map(lambda x: tuple(map(lambda y: int(y * 255),x)),RGB_tuples)
    HEX_tuples = map(lambda x: tuple(map(lambda y: chr(y).encode('hex'),x)), RGB_tuples)
    HEX_tuples = map(lambda x: "".join(x), HEX_tuples)
    HEX_tuples = map(lambda x: "#"+x.upper(), HEX_tuples)
    return HEX_tuples


Correct usage: python <script_name>.py <GenBank File>\n
\tAdditional Flags (all in pixels units): 
\t  -H <ArrowHeight> -E <HeadEdge> -l <HeadLength> -X <marginX> -Y <marginY> -S <scaling>\n
\t ArrowHeight ... the width of the arrow central part
\t HeadEdge    ... additional width of the head
\t HeadLength  ... head length
\t marginX     ... left-site margins
\t marginY     ... top-site margins
\t scaling     ... scaling of px per bp (100 means 100bp/px)\n
After running the script, the SVG script will be printed on the screen. To save
it directly into the file, one can use command: 
\tpython arrower_michael.py <file> <flags> > output_file.svg
File can be then further edited in Adobe Illustrator.
''' 

