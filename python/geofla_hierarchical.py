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

ds1 = "/Users/hugo/Dropbox/data/edge-matching/hierarchical/FRA_adm/FRA_adm1.shp"
ds2 = "/Users/hugo/Dropbox/data/edge-matching/hierarchical/FRA_adm/FRA_adm2.shp"

PPREPAIR = "/Users/hugo/projects/pprepair/pprepair"


def main():
    c = fiona.open(ds1, 'r')
    
    #-- Find and print the number of geometries
    print "Number of features:", len(c)
    
    #-- Find and print the attributes
    print "Attributes:"
    for f in c.schema['properties']:
        print "\t", f

    ids = []
    for f in c:
        ids.append(f['properties']['ID_0'])
    uniqueids = set(ids)
    print len(uniqueids)
    print len(ids)

    
    # #-- store the geometries in a list
    # lsPolys = []
    # for f in c:
    #     lsPolys.append(asShape(f['geometry']))



# def main():
#     os.chdir(FOLDER)
#     lsFiles = []
#     for fIn in glob.glob('*.shp'): 
#         lsFiles.append(fIn)

#     # 1. validate each PP first
#     valid = True
#     for f in lsFiles:
#         if (validatepp(f) == False):
#             valid = False
#     if valid == False:
#         sys.exit()

#     # 2. hierarchical validation
#     for i in range(6, 1, -1):
#         j = i - 1
#         ci = fiona.open(lsFiles[i-1], 'r')
#         dSHN = {}
#         for f in ci:
#             code = f['properties']['SHN'][:(j*2)]
#             if code not in dSHN:
#                 dSHN[code] = []
#             dSHN[code].append(f)
#         # print len(dSHN)
#         cj = fiona.open(lsFiles[j-1], 'r')
#         for key in dSHN:
#             print key
#             # create folder
#             fname = "%s_%s" % (lsFiles[i-1][:-4], key)
#             if os.path.exists(fname):
#                 shutil.rmtree(fname)
#             os.mkdir(fname)
#             # write the subset of the polygons
#             name = "%s/%s.shp" % (fname, key)
#             with fiona.open(name, 
#                             'w',
#                             driver=ci.driver,
#                             crs=ci.crs,
#                             schema=ci.schema) as output:
#                 for f in dSHN[key]:
#                     output.write(f)
#             output.close()

#             polygons = []
#             for f in dSHN[key]:
#                 polygons.append(asShape(f['geometry']))
#             if len(polygons) == 0:
#                 print "ERROR: no matching parent. Abort"
#                 sys.exit()
#             union = unary_union(polygons)
#             if union.is_valid == False:
#                 print "ERROR: spatial extent not a valid polygon. Abort."
#                 sys.exit()
#             # print union.geom_type, "\n"
#             name = "%s/%s_extent.shp" % (fname, key)
#             with fiona.open(name, 
#                             'w',
#                             driver=cj.driver,
#                             crs=cj.crs,
#                             schema=cj.schema) as output:
#                 f['geometry'] = mapping(union) 
#                 output.write(f)



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
        return False
    else:
        print "\tvalid."
        return True


if __name__ == "__main__":
    main()
