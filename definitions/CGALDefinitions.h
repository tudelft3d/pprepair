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

#ifndef CGALDEFINITIONS_H
#define CGALDEFINITIONS_H

// Compile-time options
#define EXACT_CONSTRUCTIONS       // Exact arithmetic: memory and processing time increase
//#define TRIANGULATION_HIERARCHY   // Faster point location algorithm: more memory

// CGAL kernel
#ifdef EXACT_CONSTRUCTIONS
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#else
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#endif

// CGAL values
#include <CGAL/enum.h>

// CGAL classes
#include <CGAL/Polygon_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_hierarchy_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

// VertexInfo for SAFE
#ifdef USE_VERTEX_INFO
#include <VertexInfo.h>
typedef std::vector<VertexInfo> RingInfo;
#endif

// Kernel
#ifdef EXACT_CONSTRUCTIONS
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
#else
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
#endif


#include <CGAL/squared_distance_2.h>

// Low level stuff
#ifdef TRIANGULATION_HIERARCHY
#ifdef USE_VERTEX_INFO
typedef CGAL::Triangulation_vertex_base_with_info_2<VertexInfo,K> TVB;
#else
typedef CGAL::Triangulation_vertex_base_2<K> TVB;
#endif
typedef CGAL::Triangulation_hierarchy_vertex_base_2<TVB> VB;
#else
typedef CGAL::Triangulation_vertex_base_2<K> VB;
#endif
typedef CGAL::Constrained_triangulation_face_base_2<K> FB;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo, K, FB> FBWI;
typedef CGAL::Triangulation_data_structure_2<VB, FBWI> TDS;
typedef CGAL::Exact_predicates_tag PT;
typedef CGAL::Exact_intersections_tag IT;
#ifdef EXACT_CONSTRUCTIONS
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, IT> CDT;
#else
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, PT> CDT;
#endif
#ifdef TRIANGULATION_HIERARCHY
typedef CGAL::Triangulation_hierarchy_2<CDT> CDTH;
typedef CGAL::Constrained_triangulation_plus_2<CDTH> Triangulation;
#else
typedef CGAL::Constrained_triangulation_plus_2<CDT> Triangulation;
#endif

// Other types, for easy reading
typedef Triangulation::Point Point;
typedef Triangulation::Segment Segment;
typedef CGAL::Polygon_2<K> Ring;

// Non CGAL types
typedef std::vector<std::pair<std::vector<Triangulation::Vertex_handle>, std::vector<std::vector<Triangulation::Vertex_handle> > > > TaggingVector;

// Polygon type to avoid CGAL's Polygon_with_holes_2
class Polygon {
public:
  typedef std::vector<Ring>::const_iterator Hole_const_iterator;
  
  Polygon(Ring &outer, std::vector<Ring>::iterator innerBegin, std::vector<Ring>::iterator innerEnd) {
    if (outer.is_counterclockwise_oriented() == true) {
      outer.reverse_orientation();
//      std::cout << "orientation reversed." << std::endl;
    }
    outerRing = outer;
    for (std::vector<Ring>::iterator f = innerBegin; f != innerEnd; f++) {
      if (f->is_counterclockwise_oriented() == false) {
        f->reverse_orientation();
//        std::cout << "orientation inner reversed." << std::endl;
      }
    }
    innerRings = std::vector<Ring>(innerBegin, innerEnd);
  }
  
  const Ring &outer_boundary() const {
    return outerRing;
  }
  
  Hole_const_iterator holes_begin() const {
    return innerRings.begin();
  }
  
  Hole_const_iterator holes_end() const {
    return innerRings.end();
  }
  
private:
  Ring outerRing;
  std::vector<Ring> innerRings;
};

#endif