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
    bool tagTriangulation(Triangulation &triangulation, TaggingVector &edgesToTag);
    bool makeAllHolesValid(Triangulation &triangulation);
    bool splitRegions(Triangulation &triangulation, double ratio);
    bool repairTrianglesByNumberOfNeighbours(Triangulation &triangulation, bool alsoUniverse);
	bool repairTrianglesByAbsoluteMajority(Triangulation &triangulation, bool alsoUniverse);
	bool repairTrianglesByLongestBoundary(Triangulation &triangulation, bool alsoUniverse);
	bool repairRegionsByLongestBoundary(Triangulation &triangulation, bool alsoUniverse);
	bool repairRegionsByRandomNeighbour(Triangulation &triangulation, bool alsoUniverse);
	bool repairByPriorityList(Triangulation &triangulation, const char *file);
    bool repairEdgeMatching(Triangulation &triangulation, const char *file);
    bool matchSchemata(Triangulation &triangulation);
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
    
    // Internal special tags
	PolygonHandle universe;
    
    // Cached values
    Triangulation::Face_handle startingSearchFace, startingSearchFaceInRing;  // faces that are expected to be close to the next point to be added
    
    // Helper functions
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