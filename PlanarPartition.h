/*
 Copyright (c) 2009-2012, 
 Gustavo Adolfo Ken Arroyo Ohori    g.a.k.arroyoohori@tudelft.nl
 Hugo Ledoux                        h.ledoux@tudelft.nl
 Martijn Meijers                    b.m.meijers@tudelft.nl
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met: 
 
 1. Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer. 
 2. Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution. 
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef PLANARPARTITION_H
#define PLANARPARTITION_H

#include "IOWorker.h"

class PlanarPartition {
public:
    // Constructors and destructors, initialisation
	PlanarPartition();
	~PlanarPartition();
    
    // Operations
    bool addToTriangulation(const char *file, unsigned int schemaIndex = 0);
    
    bool tagTriangulation();
    bool makeAllHolesValid();
    bool addAllowedHole(Point p);
    bool addAllowedHoles(const char *file);
    bool splitRegions(double ratio);
    
    bool checkValidity();
    bool repairTrianglesByNumberOfNeighbours(bool alsoUniverse);
	bool repairTrianglesByAbsoluteMajority(bool alsoUniverse);
	bool repairTrianglesByLongestBoundary(bool alsoUniverse);
	bool repairRegionsByLongestBoundary(bool alsoUniverse);
	bool repairRegionsByRandomNeighbour(bool alsoUniverse);
	bool repairByPriorityList(const char *file);
    bool repairEdgeMatching(const char *file);
    
    bool matchSchemata();
    
    bool reconstructPolygons(bool removeVertices = false);
    
    bool exportPolygons(const char *file, bool withProvenance);
    bool exportTriangulation(const char *file, bool withNumberOfTags, bool withFields, bool withProvenance);
    
    void printInfo();
    
private:  // Comment to have access to the triangulation and other data structures from outside
	// Internal states
	enum State {
		CREATED,
		TRIANGULATED,
		TAGGED,
		REPAIRED,
		RECONSTRUCTED
	};
    State state;
    
    // I/O handler
    IOWorker io;
    
    // Generated stuff
	Triangulation triangulation;
    TaggingVector edgesToTag;
    std::vector<std::pair<PolygonHandle *, Polygon> > outputPolygons;
};

#endif