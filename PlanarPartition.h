/*
 Copyright (c) 2009-2022,
 Ken Arroyo Ohori    k.ohori@tudelft.nl
 Hugo Ledoux         h.ledoux@tudelft.nl
 Martijn Meijers     b.m.meijers@tudelft.nl
 All rights reserved.
 
 This file is part of pprepair: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 Licensees holding a valid commercial license may use this file in
 accordance with the commercial license agreement provided with
 the software.
 
 This file is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
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
