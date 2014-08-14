# SHN = FR0000000000
# SHN = 112233445566

import sys
import os
import glob
import fiona
import subprocess
from shapely.geometry import Polygon
from shapely.geometry import Point
from shapely.geometry import asShape
from shapely.geometry import mapping
# from shapely.ops import cascaded_union
from shapely.ops import unary_union

#-- Global variables
FOLDER   = "/Users/hugo/temp/France/"
PPREPAIR = "/Users/hugo/projects/pprepair/pprepair"
    

def main():
    os.chdir(FOLDER)
    lsFiles = []
    for fIn in glob.glob('*.shp'): 
        lsFiles.append(fIn)

    # 1. validate each PP first
    # for f in lsFiles:
    #     validatepp(f)

    # 2. hierarchical validation
    for i in range(3,1,-1):
        j = i - 1

        c = fiona.open(lsFiles[i-1], 'r')
        dSHN = {}
        for f in c:
            code = f['properties']['SHN'][:(j*2)]
            if code not in dSHN:
                dSHN[code] = []
            dSHN[code].append(asShape(f['geometry']))
        print len(dSHN.keys())
        print dSHN.keys()






        sys.exit()
        # polys.append(asShape(f['geometry']))
        union = unary_union(polys)
        if union.geom_type != 'Polygon' or union.is_valid == False:
            print "ERROR: spatial extent not a single valid polygon. Abort."
            sys.exit()

        with fiona.open('out.shp', 
                        'w',
                        driver=c.driver,
                        crs=c.crs,
                        schema=c.schema) as output:
            f['geometry'] = mapping(union) 
            output.write(f)





def validatepp(fIn):
    print "=====", fIn, "====="
    cmd = PPREPAIR + " -v -i " + fIn
    op = subprocess.Popen(cmd.split(' '), 
                          stdout=subprocess.PIPE, 
                          stderr=subprocess.PIPE)
    R = op.poll()
    if R:
       res = op.communicate()
       raise ValueError(res[1])
    o =  op.communicate()[0]
    if o.find('NOT valid.') != -1:
        print "ERROR: pp is invalid"
        print o
    else:
        print "\tvalid."


def do(infile):
    c = fiona.open(infile, 'r')
    #-- Find and print the number of geometries
    print "Number of features:", len(c)
    #-- store the geometries in a list
    lsLines = []
    for each in c:
        lsLines.append(asShape(each['geometry']))

    querypt = Point(79152, 445909)
    dist = 1e9
    closest = 0
    for i,l in enumerate(lsLines):
        d = l.distance(querypt)
        if d < dist:
            dist = d
            closest = i
    print "distance:", dist
    print "ID", closest
    print lsLines[closest]

    print "done."

    
      
if __name__ == "__main__":
    main()
