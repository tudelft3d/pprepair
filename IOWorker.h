/*
 Copyright (c) 2009-2013,
 Ken Arroyo Ohori    g.a.k.arroyoohori@tudelft.nl
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

#ifndef IOWORKER_H
#define IOWORKER_H

#include "FaceInfo.h"
#include "definitions/CGALDefinitions.h"

class IOWorker {
public:
    // Construction
    IOWorker();
    
    // Main operations
    bool addToTriangulation(Triangulation &triangulation, TaggingVector &edgesToTag, const char *file, unsigned int schemaIndex);
    bool tagTriangulation(Triangulation &triangulation, TaggingVector &edgesToTag, bool spatialExtent = false);
    bool makeAllHolesValid(Triangulation &triangulation);
    void removeAllExtentTags(Triangulation &triangulation);
    bool splitRegions(Triangulation &triangulation, double ratio);
    bool repairTrianglesByNumberOfNeighbours(Triangulation &triangulation, bool alsoUniverse);
	bool repairTrianglesByAbsoluteMajority(Triangulation &triangulation, bool alsoUniverse);
	bool repairTrianglesByLongestBoundary(Triangulation &triangulation, bool alsoUniverse);
	bool repairRegionsByLongestBoundary(Triangulation &triangulation, bool alsoUniverse);
	bool repairRegionsByRandomNeighbour(Triangulation &triangulation, bool alsoUniverse);
	bool repairByPriorityList(Triangulation &triangulation, const char *file);
    bool repairEdgeMatching(Triangulation &triangulation, const char *file);
    void repairSpatialExtent(Triangulation &triangulation);
    bool matchSchemata(Triangulation &triangulation);
    void removeConstraints(Triangulation &triangulation);
    void removeVertices(Triangulation &triangulation);
    bool reconstructPolygons(Triangulation &triangulation, std::vector<std::pair<PolygonHandle *, Polygon> > &outputPolygons);
    bool exportPolygons(std::vector<std::pair<PolygonHandle *, Polygon> > &outputPolygons, const char *file, bool withProvenance);
    bool exportTriangulation(Triangulation &t, const char *file, bool withNumberOfTags, bool withFields, bool withProvenance);
    
    // Printing functions
    void insertToStream(std::ostream &ostr, OGRFeatureDefn *layerDefinition, unsigned int indentation = 0, int schemaIndex = -1);
    void insertToStream(std::ostream &ostr, const OGRFieldType &ft);
    void insertToStream(std::ostream &ostr, const OGRwkbGeometryType &gt);
    void insertTriangulationInfo(std::ostream &ostr, const Triangulation &t);
    
private:
    // Data structures
    struct FieldDescriptor {
		char *file;
		int layer;
		int field;
		
		FieldDescriptor(char *fl, int l, int fld) {
			file = fl;
			layer = l;
			field = fld;
		}
		
		bool operator<(const FieldDescriptor &fe) const {
			//std::cout << "COMP {" << std::endl << "\t" << file << "\t" << layer << "\t" << field << std::endl << "\t" << fe.file << "\t" << fe.layer << "\t" << fe.field << std::endl;
			if (file < fe.file) return true;
			if (file > fe.file) return false;
			if (layer < fe.layer) return true;
			if (layer > fe.layer) return false;
			if (field < fe.field) return true;
			return false;
		}
	};
    
    struct FieldDefinition {
		char *name;
		OGRFieldType type;
		OGRJustification justification;
		int width;
		int precision;
		
		FieldDefinition(const char *n, OGRFieldType t, OGRJustification j, int w, int p) {
			name = new char[strlen(n)+1];
			strcpy(name, n);
			type = t;
			justification = j;
			width = w;
			precision = p;
		}
		
		~FieldDefinition() {
			delete name;
		}
		
		bool matches(FieldDefinition *f) {
			// For the moment, only exact matches.
			//  if we had reference data, we could match them together (change values to the ones of the reference and return true)
			if (strcmp(name, f->name) != 0) return false;
			if (type != f->type) return false;
			if (justification != f->justification) return false;
			if (width != f->width) return false;
			if (precision != f->precision) return false;
			return true;
		}
	};
    
    struct FieldComparator {
		bool operator() (Field * const &f1, Field * const &f2) const {
			return (*f1) < (*f2);
		}
	};
    
    // What is kept from input
	std::vector<char *> fileNames;
    std::vector<PolygonHandle *> polygons;
    OGRFieldType schemaFieldType;
    std::vector<FieldDefinition *> fields;
    std::map<FieldDescriptor, unsigned int> fieldEquivalencies;
    OGRSpatialReference spatialReference;
    
    // Internal special tags
	PolygonHandle universe;
    PolygonHandle extent;
    
    // Cached values
    Triangulation::Face_handle startingSearchFace, startingSearchFaceInRing;  // faces that are expected to be close to the next point to be added
    
    // Helper functions
    void expandTriangleIntoRegion(Triangulation::Finite_faces_iterator &currentFace, Triangulation &triangulation, std::set<Triangulation::Face_handle> &facesInRegion, std::set<Triangulation::Face_handle> &processedFaces);
    unsigned int removeDuplicateVertices(std::list<Point> &ring);
    std::vector<Ring *> splitRing(Ring &ring);
    void testRings(std::vector<Ring *> &outerRings, std::vector<Ring *> &innerRings, std::vector<std::vector<Ring> > &classification, long fid);
    void copyFields(OGRFeature *ogrfeature, PolygonHandle *handle);
    void tagStack(std::stack<Triangulation::Face_handle> &positiveStack, std::stack<Triangulation::Face_handle> &negativeStack, PolygonHandle *positiveHandle, PolygonHandle *negativeHandle);
    void tagStack(std::stack<Triangulation::Face_handle> &stack, PolygonHandle *handle);
    std::list<Triangulation::Vertex_handle> *getBoundary(Triangulation::Face_handle face, int edge, PolygonHandle *polygon);
    void addtoCount(std::map<PolygonHandle *, unsigned int> &count, PolygonHandle *ph);
    void addToLength(std::map<PolygonHandle *, double> &lengths, PolygonHandle *ph, double length);
};

#endif