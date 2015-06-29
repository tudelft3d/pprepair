# Hierarchical Edge-Matching

import sys
import os
import glob
import fiona
import shutil
import subprocess
from shapely.geometry import Polygon
from shapely.geometry import Point
from shapely.geometry import asShape
from shapely.geometry import mapping
from shapely.ops import unary_union

tmpfolder = "/Users/hugo/temp/000/"

ds1 = "/Users/hugo/Dropbox/data/edge-matching/hierarchical/cbs2009/gem_2009_gn3.shp"
ds2 = "/Users/hugo/Dropbox/data/edge-matching/hierarchical/cbs2009/wijk_2009_gn3.shp"

PPREPAIR = "/Users/hugo/projects/pprepair/pprepair"


geoms1 = {}
geoms2 = {}
c1 = 0
c2 = 0
onefeature1 = 0
onefeature2 = 0

def main():
    readinputs()
    invalids = []
    for gmcode in geoms1:
        create_subset_shp(gmcode)
        if (validatepp(gmcode) == False):
            invalids.append(gmcode)
    print len(invalids), "/", len(geoms1)
    print "done."


def create_subset_shp(gmcode):
    os.chdir(tmpfolder)
    if not os.path.exists(tmpfolder):
        os.mkdir(tmpfolder)
    else:
        shutil.rmtree(tmpfolder)
        os.mkdir(tmpfolder)
    with fiona.open(tmpfolder+'e.shp', 
                    'w',
                    driver=c1.driver,
                    crs=c1.crs,
                    schema=c1.schema) as output:
        onefeature1['geometry'] = mapping(geoms1[gmcode]) 
        output.write(onefeature1)
    output.close()

    with fiona.open(tmpfolder+'p.shp', 
                    'w',
                    driver=c2.driver,
                    crs=c2.crs,
                    schema=c2.schema) as output:
        for g in geoms2[gmcode]:
            onefeature2['geometry'] = mapping(g) 
            output.write(onefeature2)
    output.close()


def readinputs():
    global c1, c2
    global onefeature1, onefeature2
    c1 = fiona.open(ds1, 'r')
    for f in c1:
        geoms1[f['properties']['GM_CODE']] = asShape(f['geometry'])
    onefeature1 = f
    c2 = fiona.open(ds2, 'r')
    for f in c2:
        code = f['properties']['GM_CODE']
        if code in geoms2:
            geoms2[code].append(asShape(f['geometry']))
        else:
            geoms2[code] = [asShape(f['geometry'])]
    onefeature2 = f


def validatepp(gmcode):
    print "=====", gmcode, "====="
    cmd = []
    cmd.append(PPREPAIR)
    cmd.append("-i")
    cmd.append("/Users/hugo/temp/000/p.shp")
    cmd.append("-e")
    cmd.append("/Users/hugo/temp/000/e.shp")
    cmd.append("-v")

    op = subprocess.Popen(cmd, 
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
        return False
    else:
        print "\tvalid."
        return True


if __name__ == "__main__":
    main()
