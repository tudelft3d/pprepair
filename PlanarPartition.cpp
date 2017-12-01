/*
 Copyright (c) 2009-2014,
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

#include "PlanarPartition.h"


PlanarPartition::PlanarPartition() {
	// Registers drivers for all supported formats in OGR
	OGRRegisterAll();
	// Set internal states
	state = CREATED;
  hasExtent = false;
//  _bbox
	// std::cout precision (for debugging)
	std::cout.setf(std::ios::fixed,std::ios::floatfield);
	std::cout.precision(6);
}


PlanarPartition::~PlanarPartition() {
	triangulation.clear();
}

bool PlanarPartition::hasSpatialExtent() {
	return hasExtent;
}

int PlanarPartition::noPolygons() {
  return polygons.size();
}

bool PlanarPartition::addOGRdatasetExtent(std::string &file) {
  if (state > TRIANGULATED) {
    std::cerr << "Error: The triangulation has already been tagged. It cannot be modified!" << std::endl;
		return false;
	}
  std::cout << "Adding spatial extent dataset: " << file << std::endl;
  #if GDAL_VERSION_MAJOR < 2
    OGRDataSource *dataSource = OGRSFDriverRegistrar::Open(file.c_str(), false);
  #else
    GDALDataset *dataSource = (GDALDataset*) GDALOpenEx(file.c_str(), GDAL_OF_READONLY, NULL, NULL, NULL);
  #endif
	if (dataSource == NULL) {
		std::cerr << "Error: Could not open file." << std::endl;
		return false;
	}
  OGRLayer *dataLayer = dataSource->GetLayer(0);
  std::cout << "\t(" << dataLayer->GetFeatureCount(true) << " polygons)" << std::endl;
  dataLayer->ResetReading();
  OGREnvelope bbox;
  dataLayer->GetExtent(&bbox);
  _bbox.Merge(bbox);

  //-- create a bigger bbox with the polygon as a hole
  double shift = 1000;
  OGRLinearRing oring;
  oring.addPoint(_bbox.MinX - shift, _bbox.MinY - shift);
  oring.addPoint(_bbox.MinX - shift, _bbox.MaxY + shift);
  oring.addPoint(_bbox.MaxX + shift, _bbox.MaxY + shift);
  oring.addPoint(_bbox.MaxX + shift, _bbox.MinY - shift);
  oring.addPoint(_bbox.MinX - shift, _bbox.MinY - shift);
  OGRPolygon bboxpolygon;
  bboxpolygon.addRing(&oring);
  
  OGRGeometry* se_geom = static_cast<OGRGeometry*>(&bboxpolygon);
  OGRFeature *feature;
  while ((feature = dataLayer->GetNextFeature()) != NULL) {
    switch(feature->GetGeometryRef()->getGeometryType()) {
      case wkbPolygon: 
      case wkbMultiPolygon: {
        if (feature->GetGeometryRef()->IsValid() == false) {
          std::cerr << "Error: spatial extent polygon invalid." << std::endl;
          return false;
        }
        se_geom = se_geom->Difference(feature->GetGeometryRef());
        break;
      }
      default: {
        std::cout << feature->GetGeometryRef()->getGeometryType() << std::endl;
        std::cerr << "Error: spatial extent must be a (Multi)Polygon." << std::endl;
        return false;
      }
    }
  }
  //-- create the feature 
  dataLayer->ResetReading();
  feature = dataLayer->GetNextFeature();
  feature->SetGeometry(se_geom);
  std::vector<OGRFeature*> ls;
  ls.push_back(feature->Clone());
  allFeatureDefns.push_back(feature->GetDefnRef());
  addFeatures(ls, true);
  #if GDAL_VERSION_MAJOR < 2
  	OGRDataSource::DestroyDataSource(dataSource);
  #else
  	GDALClose(dataSource);
  #endif
  hasExtent = true;
  return true;
}


bool PlanarPartition::addOGRdataset(std::string &file, bool skipvalideach) {
  // Check if we have already made changes to the triangulation
  if (state > TRIANGULATED) {
    std::cerr << "Error: The triangulation has already been tagged. It cannot be modified!" << std::endl;
		return false;
	}
  std::cout << "Reading input dataset: " << file << std::endl;
  std::vector<OGRFeature*> lsInputFeatures;
  if (getOGRFeatures(file, lsInputFeatures) == false)
    return false;
  //-- keep track of all the OGRFeatureDefn that come in, for PL/EM repair
  if (lsInputFeatures.size() > 0) {
    allFeatureDefns.push_back(lsInputFeatures[0]->GetDefnRef());
  }
  
  if (skipvalideach == false) {
    if (validateSinglePolygons(lsInputFeatures) == false)
      return false;
    std::cout << "\tDone, all polygons are valid." << std::endl;
  }
  else {
    std::cout << "\t->Individual validation of polygons skipped (AT YOUR OWN RISK!)" << std::endl;
  }
  if (addFeatures(lsInputFeatures) == false)
    return false;
  return true;
}


bool PlanarPartition::duplicateVerticesInPolygon(OGRPolygon* geometry) {
  bool valid = true;
  //-- 1. check for duplicate vertices
  //-- oring
  for (int p = 0; p < (geometry->getExteriorRing()->getNumPoints() - 1); p++) {
    Point a(geometry->getExteriorRing()->getX(p), geometry->getExteriorRing()->getY(p));
    Point b(geometry->getExteriorRing()->getX(p+1), geometry->getExteriorRing()->getY(p+1));
    if (a == b) {
      valid = false;
      std::cout << "Invalid polygon: duplicate vertices." << std::endl;
      std::cout << "Problem at (" << geometry->getExteriorRing()->getX(p) << ", " << geometry->getExteriorRing()->getY(p) << ")" << std::endl;
      break;
    }
  }
  //-- irings
  for (int r = 0; r < geometry->getNumInteriorRings(); r++) {
    for (int p = 0; p < (geometry->getInteriorRing(r)->getNumPoints() - 1); p++) {
      Point a(geometry->getInteriorRing(r)->getX(p), geometry->getInteriorRing(r)->getY(p));
      Point b(geometry->getInteriorRing(r)->getX(p+1), geometry->getInteriorRing(r)->getY(p+1));
      if (a == b) {
        valid = false;
        std::cout << "Invalid polygon: duplicate vertices." << std::endl;
        std::cout << "Problem at (" << geometry->getInteriorRing(r)->getX(p) << ", " << geometry->getInteriorRing(r)->getY(p) << ")" << std::endl;
        break;
      }
    }
  }
  return valid;
}


bool PlanarPartition::validateSinglePolygons(std::vector<OGRFeature*> &lsOGRFeatures) {
  std::cout << "\tValidating individually every polygon..." << std::endl;
  bool allvalid = true;
  int idno = 1;
  for (std::vector<OGRFeature*>::iterator it = lsOGRFeatures.begin() ; it != lsOGRFeatures.end(); ++it) {
    switch((*it)->GetGeometryRef()->getGeometryType()) {
      case wkbPolygon:
      case wkbPolygon25D: {
        OGRPolygon *geometry = (OGRPolygon *)(*it)->GetGeometryRef();
        if (geometry->IsValid() == false) {
          allvalid = false;
          std::cout << "--> Polygon #" << idno << std::endl;
        }
        else {
          if (duplicateVerticesInPolygon(geometry) == false) {
            allvalid = false;
            std::cout << "--> Polygon #" << idno << std::endl;
          }
        }
        break;
      }
      case wkbMultiPolygon:
      case wkbMultiPolygon25D: {
        OGRMultiPolygon *geometry = static_cast<OGRMultiPolygon *>((*it)->GetGeometryRef());
        for (int cur = 0; cur < geometry->getNumGeometries(); cur++) {
          OGRPolygon *p = (OGRPolygon *)geometry->getGeometryRef(cur);
          if (p->IsValid() == false) {
            allvalid = false;
            std::cout << "--> MultiPolygon #" << idno << std::endl;
          }
          else {
            if (duplicateVerticesInPolygon(p) == false) {
              allvalid = false;
              std::cout << "--> MultiPolygon #" << idno << std::endl;
            }
          }
        }
        break;
      }
      default: {
        std::cout << "UNKNOWN GEOMETRY TYPE, skipping feature." << std::endl;
      }
    }
    idno++;
  }
  return allvalid;
}

Polygon PlanarPartition::OGRPolygon2CGAL(OGRPolygon* geometry) {
  std::vector<std::list<Point> > outerRingsList;
  std::vector<std::list<Point> > innerRingsList;
  
  outerRingsList.push_back(std::list<Point>());
  // oring
  for (int currentPoint = 0; currentPoint < (geometry->getExteriorRing()->getNumPoints() - 1); currentPoint++)
    outerRingsList.back().push_back(Point(geometry->getExteriorRing()->getX(currentPoint),
                                          geometry->getExteriorRing()->getY(currentPoint)));
  // irings
  innerRingsList.reserve(geometry->getNumInteriorRings());
  for (int currentRing = 0; currentRing < geometry->getNumInteriorRings(); currentRing++) {
    innerRingsList.push_back(std::list<Point>());
    for (int currentPoint = 0; currentPoint < (geometry->getInteriorRing(currentRing)->getNumPoints() - 1); currentPoint++) {
      innerRingsList.back().push_back(Point(geometry->getInteriorRing(currentRing)->getX(currentPoint),
                                            geometry->getInteriorRing(currentRing)->getY(currentPoint)));
    }
  }
  Ring oring(outerRingsList[0].begin(), outerRingsList[0].end());
  outerRingsList.clear();
  std::vector<Ring> irings;
  for (unsigned int currentRing = 0; currentRing < innerRingsList.size(); currentRing++) {
    irings.push_back(Ring(innerRingsList[currentRing].begin(), innerRingsList[currentRing].end()));
    innerRingsList[currentRing].clear();
  }
  return Polygon(oring, irings.begin(), irings.end());
}


bool PlanarPartition::addFeatures(std::vector<OGRFeature*> &lsOGRFeatures, bool spatialextent) {
  std::cout << "\tAdding the polygons to the PP..." << std::endl;
  int i = 1;
  for (std::vector<OGRFeature*>::iterator f = lsOGRFeatures.begin() ; f != lsOGRFeatures.end(); f++) {
    if (i % 100 == 0)
      std::cout << "\tpolygon #" << i << std::endl;
    i++;
    std::vector<Polygon> polygonsVector;
    switch((*f)->GetGeometryRef()->getGeometryType()) {
      case wkbPolygon:
      case wkbPolygon25D: {
        OGRPolygon *geometry = static_cast<OGRPolygon *>((*f)->GetGeometryRef());
        polygonsVector.push_back(OGRPolygon2CGAL(geometry));
        break;
      }
      case wkbMultiPolygon:
      case wkbMultiPolygon25D: {
        OGRMultiPolygon *geometry = static_cast<OGRMultiPolygon *>((*f)->GetGeometryRef());
        // Check each polygon of the MultiPolygon
        for (int currentPolygon = 0; currentPolygon < geometry->getNumGeometries(); currentPolygon++) {
          OGRPolygon *g = static_cast<OGRPolygon *>(geometry->getGeometryRef(currentPolygon));
          polygonsVector.push_back(OGRPolygon2CGAL(g));
        }
        break;
      }
      default:
        std::cerr << "\tFeature #" << (*f)->GetFID() << ": unsupported type (";
        std::cerr << "). Skipped." << std::endl;
        continue;
        break;
    }
    
    for (std::vector<Polygon>::iterator currentPolygon = polygonsVector.begin(); currentPolygon != polygonsVector.end(); currentPolygon++) {
      PolygonHandle *handle;
      if (spatialextent == true)
        handle = &_extenttag;
      else
        handle = new PolygonHandle(*f);
      polygons.push_back(handle);

      // Create edges vector for this handle
      edgesToTag.push_back(std::pair<std::vector<Triangulation::Vertex_handle>, std::vector<std::vector<Triangulation::Vertex_handle> > >());

      // Insert edges into the triangulation and edges vector
      for (Ring::Edge_const_iterator currentEdge = currentPolygon->outer_boundary().edges_begin();
           currentEdge != currentPolygon->outer_boundary().edges_end();
           ++currentEdge) {
        Triangulation::Vertex_handle sourceVertex = triangulation.insert(currentEdge->source(), startingSearchFace);
        startingSearchFace = triangulation.incident_faces(sourceVertex);
        Triangulation::Vertex_handle targetVertex = triangulation.insert(currentEdge->target(), startingSearchFace);
        triangulation.insert_constraint(sourceVertex, targetVertex);
        startingSearchFace = triangulation.incident_faces(targetVertex);
        edgesToTag.back().first.push_back(sourceVertex);
      }
      for (Polygon::Hole_const_iterator currentRing = currentPolygon->holes_begin();
           currentRing != currentPolygon->holes_end();
           ++currentRing) {
        edgesToTag.back().second.push_back(std::vector<Triangulation::Vertex_handle>());
        for (Ring::Edge_const_iterator currentEdge = currentRing->edges_begin(); currentEdge != currentRing->edges_end(); ++currentEdge) {
          Triangulation::Vertex_handle sourceVertex = triangulation.insert(currentEdge->source(), startingSearchFace);
          startingSearchFace = triangulation.incident_faces(sourceVertex);
          Triangulation::Vertex_handle targetVertex = triangulation.insert(currentEdge->target(), startingSearchFace);
          triangulation.insert_constraint(sourceVertex, targetVertex);
          startingSearchFace = triangulation.incident_faces(targetVertex);
          edgesToTag.back().second.back().push_back(sourceVertex);
        }
      }
    }
    // Free (some) memory
    polygonsVector.clear();
    (*f)->SetGeometry(NULL); //-- we keep that to replace the PolygonHandle
  }
  if (triangulation.number_of_faces() > 0) {
    state = TRIANGULATED;
  }
  return true;
}


bool PlanarPartition::getOGRFeatures(std::string file, std::vector<OGRFeature*> &lsOGRFeatures) {
  #if GDAL_VERSION_MAJOR < 2
    OGRDataSource *dataSource = OGRSFDriverRegistrar::Open(file.c_str(), false);
  #else
    GDALDataset *dataSource = (GDALDataset*) GDALOpenEx(file.c_str(), GDAL_OF_READONLY, NULL, NULL, NULL);
  #endif
	if (dataSource == NULL) {
		std::cerr << "Error: Could not open file." << std::endl;
		return false;
	}
	int numberOfLayers = dataSource->GetLayerCount();
  for (int currentLayer = 0; currentLayer < numberOfLayers; currentLayer++) {
    OGRLayer *dataLayer = dataSource->GetLayer(currentLayer);
    //-- get bbox of dataset and update global one
    OGREnvelope bbox;
    dataLayer->GetExtent(&bbox);
    _bbox.Merge(bbox);
    dataLayer->ResetReading();
    unsigned int numberOfPolygons = dataLayer->GetFeatureCount(true);
    std::cout << "\tReading layer #" << currentLayer+1 << " (" << numberOfPolygons << " polygons)" << std::endl;
    
    OGRFeature *feature;
    while ((feature = dataLayer->GetNextFeature()) != NULL) {
      switch(feature->GetGeometryRef()->getGeometryType()) {
        case wkbPolygon:
        case wkbPolygon25D:
        case wkbMultiPolygon:
        case wkbMultiPolygon25D:{
          lsOGRFeatures.push_back(feature->Clone());
          break;
        }
        default: {
          std::cout << "UNKNOWN GEOMETRY TYPE, skipping feature." << std::endl;
        }
      }
    }
  }
  // Free OGR data source
  #if GDAL_VERSION_MAJOR < 2
    OGRDataSource::DestroyDataSource(dataSource);
  #else
    GDALClose(dataSource);
  #endif
//  std::cout << "\tdone." << std::endl;
  return true;
}



bool PlanarPartition::buildPP() {
  if (state < TRIANGULATED) {
    std::cout << "No triangulation to tag!" << std::endl;
    return false;
  } if (state > TRIANGULATED) {
    std::cout << "Triangulation already tagged!" << std::endl;
    return false;
  }
  
  std::cout << "Building the PP (tagging the triangles)...";
  std::stack<Triangulation::Face_handle> stack;
  Triangulation::Vertices_in_constraint_iterator previousVertex, currentVertex;
  Triangulation::Face_handle currentFace;
  int incident;
  bool sameOrder;
  // Add all edges of a polygon
  for (unsigned int currentPolygon = 0; currentPolygon < edgesToTag.size(); ++currentPolygon) {
    
    // Outer boundary
    for (unsigned int currentEdge = 0; currentEdge < edgesToTag[currentPolygon].first.size(); ++currentEdge) {
      previousVertex = triangulation.vertices_in_constraint_begin(edgesToTag[currentPolygon].first[currentEdge],
                                                                  edgesToTag[currentPolygon].first[(currentEdge+1)%edgesToTag[currentPolygon].first.size()]);
      // Check if the returned order is the same
      if ((*previousVertex)->point() == edgesToTag[currentPolygon].first[currentEdge]->point())
        sameOrder = true;
      else
        sameOrder = false;
      currentVertex = previousVertex;
      ++currentVertex;
      while (currentVertex != triangulation.vertices_in_constraint_end(edgesToTag[currentPolygon].first[currentEdge],
                                                                       edgesToTag[currentPolygon].first[(currentEdge+1)%edgesToTag[currentPolygon].first.size()])) {
        if (sameOrder) {
          if (!triangulation.is_edge(*previousVertex, *currentVertex, currentFace, incident)) {
            std::cout << "\tError: Cannot find adjoining face to an edge from the edge list!" << std::endl;
            return false;
          }
        } else {
          if (!triangulation.is_edge(*currentVertex, *previousVertex, currentFace, incident)) {
            std::cout << "\tError: Cannot find adjoining face to an edge from the edge list!" << std::endl;
            return false;
          }
        } previousVertex = currentVertex;
        ++currentVertex;
        stack.push(currentFace);
      }
    }
    // Free memory for boundaries
    edgesToTag[currentPolygon].first.clear();
    
    // Inner boundaries
    for (unsigned int currentRing = 0; currentRing < edgesToTag[currentPolygon].second.size(); ++currentRing) {
      for (unsigned int currentEdge = 0; currentEdge < edgesToTag[currentPolygon].second[currentRing].size(); ++currentEdge) {
        previousVertex = triangulation.vertices_in_constraint_begin(edgesToTag[currentPolygon].second[currentRing].at(currentEdge),
                                                                    edgesToTag[currentPolygon].second[currentRing].at((currentEdge+1)%edgesToTag[currentPolygon].second[currentRing].size()));
        // Check if the returned order is the same
        if ((*previousVertex)->point() == edgesToTag[currentPolygon].second[currentRing].at(currentEdge)->point()) sameOrder = true;
        else sameOrder = false;
        currentVertex = previousVertex;
        ++currentVertex;
        while (currentVertex != triangulation.vertices_in_constraint_end(edgesToTag[currentPolygon].second[currentRing].at(currentEdge),
                                                                         edgesToTag[currentPolygon].second[currentRing].at((currentEdge+1)%edgesToTag[currentPolygon].second[currentRing].size()))) {
          if (sameOrder) {
            if (!triangulation.is_edge(*previousVertex, *currentVertex, currentFace, incident)) {
              std::cout << "Error: No edge found!" << std::endl;
              return false;
            }
          } else {
            if (!triangulation.is_edge(*currentVertex, *previousVertex, currentFace, incident)) {
              std::cout << "Error: No edge found!" << std::endl;
              return false;
            }
          } previousVertex = currentVertex;
          currentVertex++;
          stack.push(currentFace);
        }
      }
    }
    // Free memory for inner boundary
    edgesToTag[currentPolygon].second.clear();
    // Expand the tags: special handling of the spatialExtent tag if needed
    tagStack(stack, polygons[currentPolygon]);
//    if ( (hasExtent == true) && (currentPolygon+1 == edgesToTag.size()) )
//      tagStack(stack, &_extenttag);
//    else
//      tagStack(stack, polygons[currentPolygon]);
  }
  
  // Free remaining memory
  edgesToTag.clear();
  // Tag the universe
  currentFace = triangulation.infinite_face();
  stack.push(currentFace);
  tagStack(stack, &universetag);
  
  state = TAGGED;
  isValid();
  std::cout << " done." << std::endl;
  return true;
}


void PlanarPartition::removeAllExtentTags() {
	for (Triangulation::Finite_faces_iterator currentFace = triangulation.finite_faces_begin(); currentFace != triangulation.finite_faces_end(); ++currentFace) {
		if (currentFace->info().hasTag(&_extenttag)) {
      currentFace->info().removeAllTags();
			currentFace->info().addTag(&universetag);
		}
	}
}

void PlanarPartition::tagStack(std::stack<Triangulation::Face_handle> &stack, PolygonHandle *handle) {
	while (!stack.empty()) {
		Triangulation::Face_handle currentFace = stack.top();
		stack.pop();
    if (!currentFace->info().hasTag(handle))
      currentFace->info().addTag(handle);
		if (!currentFace->neighbor(0)->info().hasTag(handle) && !currentFace->is_constrained(0)) {
			currentFace->neighbor(0)->info().addTag(handle);
			stack.push(currentFace->neighbor(0));
		}
    if (!currentFace->neighbor(1)->info().hasTag(handle) && !currentFace->is_constrained(1)) {
			currentFace->neighbor(1)->info().addTag(handle);
			stack.push(currentFace->neighbor(1));
		}
    if (!currentFace->neighbor(2)->info().hasTag(handle) && !currentFace->is_constrained(2)) {
			currentFace->neighbor(2)->info().addTag(handle);
			stack.push(currentFace->neighbor(2));
		}
	}
}


void PlanarPartition::repairSpatialExtent() {
	// Use a temporary vector to make it deterministic and order independent
	std::vector<std::pair<Triangulation::Face_handle, PolygonHandle *> > facesToRepair;
	std::set<Triangulation::Face_handle> processedFaces;
	for (Triangulation::Finite_faces_iterator currentFace = triangulation.finite_faces_begin(); currentFace != triangulation.finite_faces_end(); ++currentFace) {
		if (!currentFace->info().hasOneTag() && !processedFaces.count(currentFace)) {
			//-- Expand this triangle into a complete region
			std::set<Triangulation::Face_handle> facesInRegion;
      expandTriangleIntoRegion(currentFace, facesInRegion, processedFaces);
			PolygonHandle *tagToAssign = NULL;
      bool errorrelatedtoextent = false;
      if (currentFace->info().hasNoTags()) {
        //-- a gap
        bool extentgap = false;
        //-- check if the gap is neighbouring to an extent
        for (std::set<Triangulation::Face_handle>::const_iterator cur = facesInRegion.begin(); cur != facesInRegion.end(); cur++) {
          for (int i = 0; i < 3; i++) {
            if ( (*cur)->neighbor(i)->info().hasTag(&_extenttag) == true) {
              extentgap = true;
              break;
            }
          }
        }
        if (extentgap == true) {
          while (true) {
       			//-- Find a random tag among the direct neighbours of the region
            std::set<Triangulation::Face_handle>::iterator randomFace = facesInRegion.begin();
            std::advance(randomFace, rand()%facesInRegion.size());
            int neighbourIndex = rand()%3;
            unsigned int numberOfTags = (*randomFace)->neighbor(neighbourIndex)->info().numberOfTags();
            if (numberOfTags == 0) continue;
            if (numberOfTags == 1) {
              tagToAssign = (*randomFace)->neighbor(neighbourIndex)->info().getTags();
              if ( (tagToAssign != &universetag) && (tagToAssign != &_extenttag) )
                break;
            }
            else {
              std::list<PolygonHandle *>::const_iterator randomTag = static_cast<MultiPolygonHandle *>((*randomFace)->neighbor(neighbourIndex)->info().getTags())->getHandles()->begin();
              std::advance(randomTag, rand()%numberOfTags);
              tagToAssign = *randomTag;
              if ( (tagToAssign != &universetag) && (tagToAssign != &_extenttag) )
                break;
            }
          }
        }
        errorrelatedtoextent = extentgap;
      }
      else {
        //-- an overlap
        std::set<Triangulation::Face_handle>::iterator cur = facesInRegion.begin();
        if ( (*cur)->info().hasTag(&_extenttag) ) {
          tagToAssign = &_extenttag;
          errorrelatedtoextent = true;
        }
      }
			// Assign the region to the random tag
      if (errorrelatedtoextent == true) {
        for (std::set<Triangulation::Face_handle>::iterator currentFaceInRegion = facesInRegion.begin(); currentFaceInRegion != facesInRegion.end(); ++currentFaceInRegion) {
          facesToRepair.push_back(std::pair<Triangulation::Face_handle, PolygonHandle *>(*currentFaceInRegion, tagToAssign));
        }
      }
		}
	}
	// Re-tag faces in the vector
	for (std::vector<std::pair<Triangulation::Face_handle, PolygonHandle *> >::iterator currentFace = facesToRepair.begin(); currentFace != facesToRepair.end(); ++currentFace) {
		currentFace->first->info().removeAllTags();
		currentFace->first->info().addTag(currentFace->second);
	}
}


void PlanarPartition::expandTriangleIntoRegion(Triangulation::Finite_faces_iterator &currentFace,
                                               std::set<Triangulation::Face_handle> &facesInRegion,
                                               std::set<Triangulation::Face_handle> &processedFaces) {
  facesInRegion.insert(currentFace);
  std::stack<Triangulation::Face_handle> facesToProcess;
  facesToProcess.push(currentFace);
  while (facesToProcess.size() > 0) {
    Triangulation::Face_handle currentFaceInStack = facesToProcess.top();
    facesToProcess.pop();
    processedFaces.insert(currentFaceInStack);
    if (!currentFaceInStack->neighbor(0)->info().hasOneTag() && !facesInRegion.count(currentFaceInStack->neighbor(0)) &&
        !triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(currentFaceInStack, 0))) {
      facesInRegion.insert(currentFaceInStack->neighbor(0));
      facesToProcess.push(currentFaceInStack->neighbor(0));
    } if (!currentFaceInStack->neighbor(1)->info().hasOneTag() && !facesInRegion.count(currentFaceInStack->neighbor(1)) &&
          !triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(currentFaceInStack, 1))) {
      facesInRegion.insert(currentFaceInStack->neighbor(1));
      facesToProcess.push(currentFaceInStack->neighbor(1));
    } if (!currentFaceInStack->neighbor(2)->info().hasOneTag() && !facesInRegion.count(currentFaceInStack->neighbor(2)) &&
          !triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(currentFaceInStack, 2))) {
      facesInRegion.insert(currentFaceInStack->neighbor(2));
      facesToProcess.push(currentFaceInStack->neighbor(2));
    }
  }
  
}


bool PlanarPartition::repairPL(const std::string &file, bool alsoUniverse) {
  
//-- 1. Fetch the priority list in the file
  std::ifstream priofile(file.c_str(), std::ifstream::in);
  if (!priofile)
  {
    std::cout << "Priority file could not be opened." << std::endl;
		return false;
  }
  //-- each polygon must have the attribute used for repair, otherwise abort repair
  std::string att;
  std::getline(priofile, att);
  for (std::vector<OGRFeatureDefn*>::const_iterator it = allFeatureDefns.begin();
       it != allFeatureDefns.end();
       it++) {
    if ( (*it)->GetFieldIndex(att.c_str()) == -1) {
      std::cout << "File " << (*it)->GetName() << " doesn't have the attribute " << att << std::endl;
      return false;
    }
  }
  std::map<std::string, unsigned int> priorityMap;
  unsigned int c = 0;
	while (!priofile.eof()) {
    std::string value;
    std::getline(priofile, value);
    if (value != "") {
      priorityMap[value] = c;
      c++;
    }
  }
  priofile.close();

  // Use a temporary vector to make it deterministic and order independent
	std::vector<std::pair<Triangulation::Face_handle, PolygonHandle *> > facesToRepair;
	std::set<Triangulation::Face_handle> processedFaces;
	for (Triangulation::Finite_faces_iterator currentFace = triangulation.finite_faces_begin(); currentFace != triangulation.finite_faces_end(); ++currentFace) {
		if (!currentFace->info().hasOneTag() && !processedFaces.count(currentFace)) {
			// Expand this triangle into a complete region
			std::set<Triangulation::Face_handle> facesInRegion;
			facesInRegion.insert(currentFace);
			std::stack<Triangulation::Face_handle> facesToProcess;
			facesToProcess.push(currentFace);
			while (facesToProcess.size() > 0) {
				Triangulation::Face_handle currentFaceInStack = facesToProcess.top();
				facesToProcess.pop();
				processedFaces.insert(currentFaceInStack);
				if (!currentFaceInStack->neighbor(0)->info().hasOneTag() && !facesInRegion.count(currentFaceInStack->neighbor(0)) &&
            !triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(currentFaceInStack, 0))) {
					facesInRegion.insert(currentFaceInStack->neighbor(0));
					facesToProcess.push(currentFaceInStack->neighbor(0));
				} if (!currentFaceInStack->neighbor(1)->info().hasOneTag() && !facesInRegion.count(currentFaceInStack->neighbor(1)) &&
              !triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(currentFaceInStack, 1))) {
					facesInRegion.insert(currentFaceInStack->neighbor(1));
					facesToProcess.push(currentFaceInStack->neighbor(1));
				} if (!currentFaceInStack->neighbor(2)->info().hasOneTag() && !facesInRegion.count(currentFaceInStack->neighbor(2)) &&
              !triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(currentFaceInStack, 2))) {
					facesInRegion.insert(currentFaceInStack->neighbor(2));
					facesToProcess.push(currentFaceInStack->neighbor(2));
				}
			}
			
			// Find the tag with the highest priority
			PolygonHandle *tagToAssign = NULL;
			unsigned int priorityOfTag = UINT_MAX;
      std::map<std::string, unsigned int>::const_iterator itatt;
			for (std::set<Triangulation::Face_handle>::iterator currentFaceInRegion = facesInRegion.begin(); currentFaceInRegion != facesInRegion.end(); ++currentFaceInRegion) {
      //-- Gaps --
        if ((*currentFaceInRegion)->info().hasNoTags()) {
          for (int j = 0; j <= 2; j++) {
            if (!(*currentFaceInRegion)->neighbor(j)->info().hasNoTags()) {
              if ((*currentFaceInRegion)->neighbor(j)->info().hasOneTag() && (*currentFaceInRegion)->neighbor(j)->info().getTags() != &universetag) {
                std::string v = (*currentFaceInRegion)->neighbor(j)->info().getTags()->getValueAttributeAsString(att);
                itatt = priorityMap.find(v);
                if ( (itatt != priorityMap.end()) && (itatt->second < priorityOfTag) ) {
                  priorityOfTag = itatt->second;
                  tagToAssign = (*currentFaceInRegion)->neighbor(j)->info().getTags();
                }
              }
              else {
                MultiPolygonHandle *handle = static_cast<MultiPolygonHandle *>((*currentFaceInRegion)->neighbor(j)->info().getTags());
                for (std::list<PolygonHandle *>::const_iterator currentTag = handle->getHandles()->begin(); currentTag != handle->getHandles()->end(); ++currentTag) {
                  if (*currentTag == &universetag)
                    continue;
                  std::string v = (*currentTag)->getValueAttributeAsString(att);
                  itatt = priorityMap.find(v);
                  if ( (itatt != priorityMap.end()) && (itatt->second < priorityOfTag) ) {
                    priorityOfTag = itatt->second;
                    tagToAssign = *currentTag;
                  }
                }
              }
            }
          }
        }
			//-- Overlap
				else {
					if ((*currentFaceInRegion)->info().hasOneTag() && (*currentFaceInRegion)->info().getTags() != &universetag) {
            std::string v = (*currentFaceInRegion)->info().getTags()->getValueAttributeAsString(att);
            itatt = priorityMap.find(v);
            if ( (itatt != priorityMap.end()) && (itatt->second < priorityOfTag) ) {
              priorityOfTag = itatt->second;
              tagToAssign = (*currentFaceInRegion)->info().getTags();
            }
					}
          else {
						MultiPolygonHandle *handle = static_cast<MultiPolygonHandle *>((*currentFaceInRegion)->info().getTags());
						for (std::list<PolygonHandle *>::const_iterator currentTag = handle->getHandles()->begin(); currentTag != handle->getHandles()->end(); ++currentTag) {
							if (*currentTag == &universetag)
                continue;
              std::string v = (*currentTag)->getValueAttributeAsString(att);
              itatt = priorityMap.find(v);
              if ( (itatt != priorityMap.end()) && (itatt->second < priorityOfTag) ) {
                priorityOfTag = itatt->second;
                tagToAssign = *currentTag;
              }
						}
					}
				}
			}
			
			// Assign the tag to the triangles in the region
			for (std::set<Triangulation::Face_handle>::iterator currentFaceInRegion = facesInRegion.begin(); currentFaceInRegion != facesInRegion.end(); ++currentFaceInRegion) {
				facesToRepair.push_back(std::pair<Triangulation::Face_handle, PolygonHandle *>(*currentFaceInRegion, tagToAssign));
			}
		}
	}
	
	// Re-tag faces in the vector
	for (std::vector<std::pair<Triangulation::Face_handle, PolygonHandle *> >::iterator currentFace = facesToRepair.begin();
       currentFace != facesToRepair.end();
       ++currentFace) {
		currentFace->first->info().removeAllTags();
		currentFace->first->info().addTag(currentFace->second);
	}
  return true;
}



bool PlanarPartition::repairEM_attribute(std::map<std::string, unsigned int> &priorityMap,
                                         std::string &att,
                                         bool alsoUniverse) {
  
  // Use a temporary vector to make it deterministic and order independent
  std::vector<std::pair<Triangulation::Face_handle, PolygonHandle *> > facesToRepair;
  std::set<Triangulation::Face_handle> processedFaces;
  for (Triangulation::Finite_faces_iterator currentFace = triangulation.finite_faces_begin(); currentFace != triangulation.finite_faces_end(); ++currentFace) {
    if (!currentFace->info().hasOneTag() && !processedFaces.count(currentFace)) {
      // Expand this triangle into a complete region
      std::set<Triangulation::Face_handle> facesInRegion;
      facesInRegion.insert(currentFace);
      std::stack<Triangulation::Face_handle> facesToProcess;
      facesToProcess.push(currentFace);
      while (facesToProcess.size() > 0) {
        Triangulation::Face_handle currentFaceInStack = facesToProcess.top();
        facesToProcess.pop();
        processedFaces.insert(currentFaceInStack);
        if (!currentFaceInStack->neighbor(0)->info().hasOneTag() && !facesInRegion.count(currentFaceInStack->neighbor(0)) &&
            !triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(currentFaceInStack, 0))) {
          facesInRegion.insert(currentFaceInStack->neighbor(0));
          facesToProcess.push(currentFaceInStack->neighbor(0));
        } if (!currentFaceInStack->neighbor(1)->info().hasOneTag() && !facesInRegion.count(currentFaceInStack->neighbor(1)) &&
              !triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(currentFaceInStack, 1))) {
          facesInRegion.insert(currentFaceInStack->neighbor(1));
          facesToProcess.push(currentFaceInStack->neighbor(1));
        } if (!currentFaceInStack->neighbor(2)->info().hasOneTag() && !facesInRegion.count(currentFaceInStack->neighbor(2)) &&
              !triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(currentFaceInStack, 2))) {
          facesInRegion.insert(currentFaceInStack->neighbor(2));
          facesToProcess.push(currentFaceInStack->neighbor(2));
        }
      }
      
      // Find the tag with the highest priority
      PolygonHandle *tagToAssign = NULL;
      unsigned int priorityOfTagg = 0;
      unsigned int priorityOfTago = UINT_MAX;
      std::map<std::string, unsigned int>::const_iterator itatt;
      for (std::set<Triangulation::Face_handle>::iterator currentFaceInRegion = facesInRegion.begin(); currentFaceInRegion != facesInRegion.end(); ++currentFaceInRegion) {
        //-- Gaps --
        if ((*currentFaceInRegion)->info().hasNoTags()) {
          for (int j = 0; j <= 2; j++) {
            if (!(*currentFaceInRegion)->neighbor(j)->info().hasNoTags()) {
              if ((*currentFaceInRegion)->neighbor(j)->info().hasOneTag() && (*currentFaceInRegion)->neighbor(j)->info().getTags() != &universetag) {
                std::string v = (*currentFaceInRegion)->neighbor(j)->info().getTags()->getValueAttributeAsString(att);
                itatt = priorityMap.find(v);
                if ( (itatt != priorityMap.end()) && (itatt->second >= priorityOfTagg) ) {
                  priorityOfTagg = itatt->second;
                  tagToAssign = (*currentFaceInRegion)->neighbor(j)->info().getTags();
                }
              }
              else {
                MultiPolygonHandle *handle = static_cast<MultiPolygonHandle *>((*currentFaceInRegion)->neighbor(j)->info().getTags());
                for (std::list<PolygonHandle *>::const_iterator currentTag = handle->getHandles()->begin(); currentTag != handle->getHandles()->end(); ++currentTag) {
                  if (*currentTag == &universetag)
                    continue;
                  std::string v = (*currentTag)->getValueAttributeAsString(att);
                  itatt = priorityMap.find(v);
                  if ( (itatt != priorityMap.end()) && (itatt->second >= priorityOfTagg) ) {
                    priorityOfTagg = itatt->second;
                    tagToAssign = *currentTag;
                  }
                }
              }
            }
          }
        }
        //-- Overlap
        else {
          if ((*currentFaceInRegion)->info().hasOneTag() && (*currentFaceInRegion)->info().getTags() != &universetag) {
            std::string v = (*currentFaceInRegion)->info().getTags()->getValueAttributeAsString(att);
            itatt = priorityMap.find(v);
            if ( (itatt != priorityMap.end()) && (itatt->second < priorityOfTago) ) {
              priorityOfTago = itatt->second;
              tagToAssign = (*currentFaceInRegion)->info().getTags();
            }
          }
          else {
            MultiPolygonHandle *handle = static_cast<MultiPolygonHandle *>((*currentFaceInRegion)->info().getTags());
            for (std::list<PolygonHandle *>::const_iterator currentTag = handle->getHandles()->begin(); currentTag != handle->getHandles()->end(); ++currentTag) {
              if (*currentTag == &universetag)
                continue;
              std::string v = (*currentTag)->getValueAttributeAsString(att);
              itatt = priorityMap.find(v);
              if ( (itatt != priorityMap.end()) && (itatt->second < priorityOfTago) ) {
                priorityOfTago = itatt->second;
                tagToAssign = *currentTag;
              }
            }
          }
        }
      }
      
      // Assign the tag to the triangles in the region
      for (std::set<Triangulation::Face_handle>::iterator currentFaceInRegion = facesInRegion.begin(); currentFaceInRegion != facesInRegion.end(); ++currentFaceInRegion) {
        facesToRepair.push_back(std::pair<Triangulation::Face_handle, PolygonHandle *>(*currentFaceInRegion, tagToAssign));
      }
    }
  }
  
  // Re-tag faces in the vector
  for (std::vector<std::pair<Triangulation::Face_handle, PolygonHandle *> >::iterator currentFace = facesToRepair.begin();
       currentFace != facesToRepair.end();
       ++currentFace) {
    currentFace->first->info().removeAllTags();
    currentFace->first->info().addTag(currentFace->second);
  }
  return true;
}

bool PlanarPartition::repairEM_dataset(std::map<std::string, unsigned int> &priorityMap, double splitregions, bool alsoUniverse) {
  if (splitregions > 0)
    return repairEM_dataset_add_constraints(priorityMap, splitregions, alsoUniverse);
  else
    return repairEM_dataset_without_constraints(priorityMap, alsoUniverse);
}

bool PlanarPartition::repairEM_dataset_without_constraints(std::map<std::string, unsigned int> &priorityMap, bool alsoUniverse) {
  // Use a temporary vector to make it deterministic and order independent
  std::vector<std::pair<Triangulation::Face_handle, PolygonHandle *> > facesToRepair;
  std::set<Triangulation::Face_handle> processedFaces;
  for (Triangulation::Finite_faces_iterator currentFace = triangulation.finite_faces_begin(); currentFace != triangulation.finite_faces_end(); ++currentFace) {
    if (!currentFace->info().hasOneTag() && !processedFaces.count(currentFace)) {
      // Expand this triangle into a complete region
      std::set<Triangulation::Face_handle> facesInRegion;
      facesInRegion.insert(currentFace);
      std::stack<Triangulation::Face_handle> facesToProcess;
      facesToProcess.push(currentFace);
      while (facesToProcess.size() > 0) {
        Triangulation::Face_handle currentFaceInStack = facesToProcess.top();
        facesToProcess.pop();
        processedFaces.insert(currentFaceInStack);
        if (!currentFaceInStack->neighbor(0)->info().hasOneTag() && !facesInRegion.count(currentFaceInStack->neighbor(0)) &&
              !triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(currentFaceInStack, 0))) {
          facesInRegion.insert(currentFaceInStack->neighbor(0));
          facesToProcess.push(currentFaceInStack->neighbor(0));
        }
        if (!currentFaceInStack->neighbor(1)->info().hasOneTag() && !facesInRegion.count(currentFaceInStack->neighbor(1)) &&
              !triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(currentFaceInStack, 1))) {
          facesInRegion.insert(currentFaceInStack->neighbor(1));
          facesToProcess.push(currentFaceInStack->neighbor(1));
        }
        if (!currentFaceInStack->neighbor(2)->info().hasOneTag() && !facesInRegion.count(currentFaceInStack->neighbor(2)) &&
              !triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(currentFaceInStack, 2))) {
          facesInRegion.insert(currentFaceInStack->neighbor(2));
          facesToProcess.push(currentFaceInStack->neighbor(2));
        }
      }
      
      // Find the tag with the highest priority
      PolygonHandle *tagToAssign = NULL;
      unsigned int priorityOfTagg = 0;
      unsigned int priorityOfTago = UINT_MAX;
      std::map<std::string, unsigned int>::const_iterator itatt;
      for (std::set<Triangulation::Face_handle>::iterator currentFaceInRegion = facesInRegion.begin(); currentFaceInRegion != facesInRegion.end(); ++currentFaceInRegion) {
        //-- Gaps --
        if ((*currentFaceInRegion)->info().hasNoTags()) {
          for (int j = 0; j <= 2; j++) {
            if (!(*currentFaceInRegion)->neighbor(j)->info().hasNoTags()) {
              if ((*currentFaceInRegion)->neighbor(j)->info().hasOneTag() && (*currentFaceInRegion)->neighbor(j)->info().getTags() != &universetag) {
                std::string v = (*currentFaceInRegion)->neighbor(j)->info().getTags()->getDSName();
//                std::string v = (*currentFaceInRegion)->neighbor(j)->info().getTags()->getValueAttributeAsString(att);
                itatt = priorityMap.find(v);
                if ( (itatt != priorityMap.end()) && (itatt->second >= priorityOfTagg) ) {
                  priorityOfTagg = itatt->second;
                  tagToAssign = (*currentFaceInRegion)->neighbor(j)->info().getTags();
                }
              }
              else {
                MultiPolygonHandle *handle = static_cast<MultiPolygonHandle *>((*currentFaceInRegion)->neighbor(j)->info().getTags());
                for (std::list<PolygonHandle *>::const_iterator currentTag = handle->getHandles()->begin(); currentTag != handle->getHandles()->end(); ++currentTag) {
                  if (*currentTag == &universetag)
                    continue;
                  std::string v = (*currentTag)->getDSName();
                  itatt = priorityMap.find(v);
                  if ( (itatt != priorityMap.end()) && (itatt->second >= priorityOfTagg) ) {
                    priorityOfTagg = itatt->second;
                    tagToAssign = *currentTag;
                  }
                }
              }
            }
          }
        }
        //-- Overlap
        else {
          if ((*currentFaceInRegion)->info().hasOneTag() && (*currentFaceInRegion)->info().getTags() != &universetag) {
            std::string v = (*currentFaceInRegion)->info().getTags()->getDSName();
            itatt = priorityMap.find(v);
            if ( (itatt != priorityMap.end()) && (itatt->second < priorityOfTago) ) {
              priorityOfTago = itatt->second;
              tagToAssign = (*currentFaceInRegion)->info().getTags();
            }
          }
          else {
            MultiPolygonHandle *handle = static_cast<MultiPolygonHandle *>((*currentFaceInRegion)->info().getTags());
            for (std::list<PolygonHandle *>::const_iterator currentTag = handle->getHandles()->begin(); currentTag != handle->getHandles()->end(); ++currentTag) {
              if (*currentTag == &universetag)
                continue;
              std::string v = (*currentTag)->getDSName();
              itatt = priorityMap.find(v);
              if ( (itatt != priorityMap.end()) && (itatt->second < priorityOfTago) ) {
                priorityOfTago = itatt->second;
                tagToAssign = *currentTag;
              }
            }
          }
        }
      }
      
      // Assign the tag to the triangles in the region
      for (std::set<Triangulation::Face_handle>::iterator currentFaceInRegion = facesInRegion.begin(); currentFaceInRegion != facesInRegion.end(); ++currentFaceInRegion) {
        facesToRepair.push_back(std::pair<Triangulation::Face_handle, PolygonHandle *>(*currentFaceInRegion, tagToAssign));
      }
    }
  }
  
  // Re-tag faces in the vector
  for (std::vector<std::pair<Triangulation::Face_handle, PolygonHandle *> >::iterator currentFace = facesToRepair.begin();
       currentFace != facesToRepair.end();
       ++currentFace) {
    currentFace->first->info().removeAllTags();
    currentFace->first->info().addTag(currentFace->second);
  }
  return true;
}


bool PlanarPartition::repairEM_dataset_add_constraints(std::map<std::string, unsigned int> &priorityMap, double splitregions, bool alsoUniverse) {
  
  add_extra_constraints_in_gaps(splitregions);
  
  // Use a temporary vector to make it deterministic and order independent
  std::vector<std::pair<Triangulation::Face_handle, PolygonHandle *> > facesToRepair;
  std::set<Triangulation::Face_handle> processedFaces;
  bool triedthemall = false;
  while (isValid() == false) {
   std::cout << "=====start round=====" << std::endl;
    facesToRepair.clear();
    processedFaces.clear();
    for (Triangulation::Finite_faces_iterator currentFace = triangulation.finite_faces_begin(); currentFace != triangulation.finite_faces_end(); ++currentFace) {
      if (!currentFace->info().hasOneTag() && !processedFaces.count(currentFace)) {
        // Expand this triangle into a complete region
        std::set<Triangulation::Face_handle> facesInRegion;
        facesInRegion.insert(currentFace);
        std::stack<Triangulation::Face_handle> facesToProcess;
        facesToProcess.push(currentFace);
        while (facesToProcess.size() > 0) {
          Triangulation::Face_handle currentFaceInStack = facesToProcess.top();
          facesToProcess.pop();
          processedFaces.insert(currentFaceInStack);
          if (!currentFaceInStack->neighbor(0)->info().hasOneTag() && !facesInRegion.count(currentFaceInStack->neighbor(0)) &&
                !triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(currentFaceInStack, 0))) {
            facesInRegion.insert(currentFaceInStack->neighbor(0));
            facesToProcess.push(currentFaceInStack->neighbor(0));
          }
          if (!currentFaceInStack->neighbor(1)->info().hasOneTag() && !facesInRegion.count(currentFaceInStack->neighbor(1)) &&
                !triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(currentFaceInStack, 1))) {
            facesInRegion.insert(currentFaceInStack->neighbor(1));
            facesToProcess.push(currentFaceInStack->neighbor(1));
          }
          if (!currentFaceInStack->neighbor(2)->info().hasOneTag() && !facesInRegion.count(currentFaceInStack->neighbor(2)) &&
                !triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(currentFaceInStack, 2))) {
            facesInRegion.insert(currentFaceInStack->neighbor(2));
            facesToProcess.push(currentFaceInStack->neighbor(2));
          }
        }
        
        // Find the tag with the highest priority
//        std::cout << "one region" << std::endl;
        PolygonHandle *tagToAssign = NULL;
        unsigned int priorityOfTagg = 0;
        unsigned int priorityOfTago = UINT_MAX;
        std::map<std::string, unsigned int>::const_iterator itatt;
        std::set<std::string> adjds;
        bool gap = false;
        for (std::set<Triangulation::Face_handle>::iterator currentFaceInRegion = facesInRegion.begin(); currentFaceInRegion != facesInRegion.end(); ++currentFaceInRegion) {
          //-- Gaps --
          if ((*currentFaceInRegion)->info().hasNoTags()) {
            gap = true;
            for (int j = 0; j <= 2; j++) {
              if (!(*currentFaceInRegion)->neighbor(j)->info().hasNoTags()) {
                if ( (*currentFaceInRegion)->neighbor(j)->info().hasOneTag() && 
                     (*currentFaceInRegion)->neighbor(j)->info().getTags() != &universetag ) {
                  std::string v = (*currentFaceInRegion)->neighbor(j)->info().getTags()->getDSName();
                  adjds.insert(v);
                  itatt = priorityMap.find(v);
                  if ( (itatt != priorityMap.end()) && (itatt->second >= priorityOfTagg) ) {
                    priorityOfTagg = itatt->second;
                    tagToAssign = (*currentFaceInRegion)->neighbor(j)->info().getTags();
                  }
                }
                else {
                  MultiPolygonHandle *handle = static_cast<MultiPolygonHandle *>((*currentFaceInRegion)->neighbor(j)->info().getTags());
                  for (std::list<PolygonHandle *>::const_iterator currentTag = handle->getHandles()->begin(); currentTag != handle->getHandles()->end(); ++currentTag) {
                    if (*currentTag == &universetag)
                      continue;
                    std::string v = (*currentTag)->getDSName();
                    adjds.insert(v);
                    itatt = priorityMap.find(v);
                    if ( (itatt != priorityMap.end()) && (itatt->second >= priorityOfTagg) ) {
                      priorityOfTagg = itatt->second;
                      tagToAssign = *currentTag;
                    }
                  }
                }
              }
            }
          }
          //-- Overlap
          else {
            if ((*currentFaceInRegion)->info().hasOneTag() && (*currentFaceInRegion)->info().getTags() != &universetag) {
              std::string v = (*currentFaceInRegion)->info().getTags()->getDSName();
              itatt = priorityMap.find(v);
              if ( (itatt != priorityMap.end()) && (itatt->second < priorityOfTago) ) {
                priorityOfTago = itatt->second;
                tagToAssign = (*currentFaceInRegion)->info().getTags();
              }
            }
            else {
              MultiPolygonHandle *handle = static_cast<MultiPolygonHandle *>((*currentFaceInRegion)->info().getTags());
              for (std::list<PolygonHandle *>::const_iterator currentTag = handle->getHandles()->begin(); currentTag != handle->getHandles()->end(); ++currentTag) {
                if (*currentTag == &universetag)
                  continue;
                std::string v = (*currentTag)->getDSName();
                itatt = priorityMap.find(v);
                if ( (itatt != priorityMap.end()) && (itatt->second < priorityOfTago) ) {
                  priorityOfTago = itatt->second;
                  tagToAssign = *currentTag;
                }
              }
            }
          }
        }
        if (gap == false) {
          for (std::set<Triangulation::Face_handle>::iterator currentFaceInRegion = facesInRegion.begin(); currentFaceInRegion != facesInRegion.end(); ++currentFaceInRegion) {
            facesToRepair.push_back(std::pair<Triangulation::Face_handle, PolygonHandle *>(*currentFaceInRegion, tagToAssign));
          }
        }
        if (gap == true) {
          if (adjds.size() > 1) {
            for (std::set<Triangulation::Face_handle>::iterator currentFaceInRegion = facesInRegion.begin(); currentFaceInRegion != facesInRegion.end(); ++currentFaceInRegion) {
              facesToRepair.push_back(std::pair<Triangulation::Face_handle, PolygonHandle *>(*currentFaceInRegion, tagToAssign));
            }
          }
          else { //-- only 1 dataset
            itatt = priorityMap.find(tagToAssign->getDSName());
            if (itatt->second > 0) {
              for (std::set<Triangulation::Face_handle>::iterator currentFaceInRegion = facesInRegion.begin(); currentFaceInRegion != facesInRegion.end(); ++currentFaceInRegion) {
                if (triedthemall)
                  std::cout << tagToAssign->getDSName() << std::endl;
                facesToRepair.push_back(std::pair<Triangulation::Face_handle, PolygonHandle *>(*currentFaceInRegion, tagToAssign));
              }
            }
          }
        }
      }
    }
    
    //-- if no changes in the round of all faces, then add labels and try again.
    std::cout << "#:" << facesToRepair.size() << std::endl;
    if (facesToRepair.size() == 0)
      triedthemall = true;
    else
      triedthemall = false;
    // Re-tag faces in the vector
    for (std::vector<std::pair<Triangulation::Face_handle, PolygonHandle *> >::iterator currentFace = facesToRepair.begin();
         currentFace != facesToRepair.end();
         ++currentFace) {
      currentFace->first->info().removeAllTags();
      currentFace->first->info().addTag(currentFace->second);
    }
  }
  return true;
}


bool PlanarPartition::repairRN(bool alsoUniverse) {
	bool repaired = true;
  // Use a temporary vector to make it deterministic and order independent
	std::vector<std::pair<Triangulation::Face_handle, PolygonHandle *> > facesToRepair;
	std::set<Triangulation::Face_handle> processedFaces;
	for (Triangulation::Finite_faces_iterator currentFace = triangulation.finite_faces_begin(); currentFace != triangulation.finite_faces_end(); ++currentFace) {
		if (!currentFace->info().hasOneTag() && !processedFaces.count(currentFace)) {
			// Expand this triangle into a complete region
			std::set<Triangulation::Face_handle> facesInRegion;
			facesInRegion.insert(currentFace);
			std::stack<Triangulation::Face_handle> facesToProcess;
			facesToProcess.push(currentFace);
			while (facesToProcess.size() > 0) {
				Triangulation::Face_handle currentFaceInStack = facesToProcess.top();
				facesToProcess.pop();
				processedFaces.insert(currentFaceInStack);
				if (!currentFaceInStack->neighbor(0)->info().hasOneTag() && !facesInRegion.count(currentFaceInStack->neighbor(0)) &&
            !triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(currentFaceInStack, 0))) {
					facesInRegion.insert(currentFaceInStack->neighbor(0));
					facesToProcess.push(currentFaceInStack->neighbor(0));
				} if (!currentFaceInStack->neighbor(1)->info().hasOneTag() && !facesInRegion.count(currentFaceInStack->neighbor(1)) &&
              !triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(currentFaceInStack, 1))) {
					facesInRegion.insert(currentFaceInStack->neighbor(1));
					facesToProcess.push(currentFaceInStack->neighbor(1));
				} if (!currentFaceInStack->neighbor(2)->info().hasOneTag() && !facesInRegion.count(currentFaceInStack->neighbor(2)) &&
              !triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(currentFaceInStack, 2))) {
					facesInRegion.insert(currentFaceInStack->neighbor(2));
					facesToProcess.push(currentFaceInStack->neighbor(2));
				}
			}
			// Find a random tag
			PolygonHandle *tagToAssign;
			while (true) {
				std::set<Triangulation::Face_handle>::iterator randomFace = facesInRegion.begin();
				std::advance(randomFace, rand() % facesInRegion.size());
				int neighbourIndex = (rand() % 3);
				unsigned int numberOfTags = (*randomFace)->neighbor(neighbourIndex)->info().numberOfTags();
				if (numberOfTags == 0)
          continue;
				if (numberOfTags == 1) {
					tagToAssign = (*randomFace)->neighbor(neighbourIndex)->info().getTags();
					if (alsoUniverse || tagToAssign != &universetag)
            break;
				}
        else {
					std::list<PolygonHandle *>::const_iterator randomTag = static_cast<MultiPolygonHandle *>((*randomFace)->neighbor(neighbourIndex)->info().getTags())->getHandles()->begin();
					std::advance(randomTag, rand()%numberOfTags);
					tagToAssign = *randomTag;
					if (alsoUniverse || tagToAssign != &universetag)
            break;
				}
			}
			// Assign the region to the random tag
			for (std::set<Triangulation::Face_handle>::iterator currentFaceInRegion = facesInRegion.begin(); currentFaceInRegion != facesInRegion.end(); ++currentFaceInRegion) {
				facesToRepair.push_back(std::pair<Triangulation::Face_handle, PolygonHandle *>(*currentFaceInRegion, tagToAssign));
			}
		}
	}
	// Re-tag faces in the vector
	for (std::vector<std::pair<Triangulation::Face_handle, PolygonHandle *> >::iterator currentFace = facesToRepair.begin(); currentFace != facesToRepair.end(); ++currentFace) {
		currentFace->first->info().removeAllTags();
		currentFace->first->info().addTag(currentFace->second);
	}
	return repaired;
}



bool PlanarPartition::repair(const std::string &method, bool alsoUniverse, const std::string &priofile, double splitregions) {
	if (state < TAGGED) {
		std::cout << "Triangulation not yet tagged. Cannot repair!" << std::endl;
		return false;
	}
  if (state > TAGGED) {
		std::cout << "PP already valid!" << std::endl;
		return false;
	}
  
  //-- repair spatial extent if present
  if (hasExtent == true)
    repairSpatialExtent();
	
  //-- repair with the specific chosen method
  time_t thisTime = time(NULL);
  bool repaired;
  if ( (method == "RN") || (method == "fix") ) {
    std::cout << "Repairing by random neighbour..." << std::endl;
    repaired = repairRN(alsoUniverse);
  }
  else if (method == "LB") {
    std::cout << "Repairing by longest boundary..." << std::endl;
    repaired = repairLB(alsoUniverse);
  }
  else if (method == "PL") {
    std::cout << "Repairing by priority list..." << std::endl;
    repaired = repairPL(priofile, alsoUniverse);
  }
  else if (method == "EM") {
    std::cout << "Repairing with edge-matching method...";
    std::map<std::string, unsigned int> priorityMap;
    std::string attr;
    if (getPriorityList(priofile, priorityMap, attr) == false) {
      return false;
    }
    if (attr == "datasets") {
      std::cout << "datasets." << std::endl;

      repaired = repairEM_dataset(priorityMap, splitregions, alsoUniverse);
    }
    else {
      std::cout << "with attributes." << std::endl;
      repaired = repairEM_attribute(priorityMap, attr, alsoUniverse);
    }
  }
  
	if (repaired) {
		std::cout << "Repair successful (" << time(NULL)-thisTime << " s). All polygons are now valid." << std::endl;
	}
  else {
		std::cout << "Repair of all polygons not possible (" << time(NULL)-thisTime << " s)." << std::endl;
	}
	//-- handling of tags for spatial extent
  if (hasExtent == true)
    removeAllExtentTags();
  if (repaired)
    state = REPAIRED;
	return repaired;
}


bool PlanarPartition::getPriorityList(const std::string &file, std::map<std::string, unsigned int> &priorityMap, std::string &att) {
  //-- 1. Fetch the priority list in the file
  std::ifstream priofile(file.c_str(), std::ifstream::in);
  if (!priofile)
  {
    std::cout << "Priority file could not be opened." << std::endl;
    return false;
  }
  //-- 2. type of edge-matching: 1st line either "datasets" OR "an_attribute"
  std::getline(priofile, att);
  if (att != "datasets") {
    //-- each polygon must have the attribute used for repair, otherwise abort repair
    for (std::vector<OGRFeatureDefn*>::const_iterator it = allFeatureDefns.begin();
         it != allFeatureDefns.end();
         it++) {
      if ( (*it)->GetFieldIndex(att.c_str()) == -1) {
        std::cout << "File " << (*it)->GetName() << " doesn't have the attribute " << att << std::endl;
        return false;
      }
    }
  }
  unsigned int c = 0;
  while (!priofile.eof()) {
    std::string value;
    std::getline(priofile, value);
    if (value != "") {
      priorityMap[value] = c;
      c++;
    }
  }
  priofile.close();
  return true;
}

bool PlanarPartition::repairLB(bool alsoUniverse) {
	bool repaired = true;
	// Use a temporary vector to make it deterministic and order independent
	std::vector<std::pair<Triangulation::Face_handle, PolygonHandle *> > facesToRepair;
	std::set<Triangulation::Face_handle> processedFaces;
	for (Triangulation::Finite_faces_iterator currentFace = triangulation.finite_faces_begin(); currentFace != triangulation.finite_faces_end(); ++currentFace) {
		std::map<PolygonHandle *, double> tagBoundaryLength;
		if (!currentFace->info().hasOneTag() && !processedFaces.count(currentFace)) {
			// Expand this triangle into a complete region
			std::set<Triangulation::Face_handle> facesInRegion;
			facesInRegion.insert(currentFace);
			std::stack<Triangulation::Face_handle> facesToProcess;
			facesToProcess.push(currentFace);
			while (facesToProcess.size() > 0) {
				Triangulation::Face_handle currentFaceInStack = facesToProcess.top();
				facesToProcess.pop();
				processedFaces.insert(currentFaceInStack);
				if (!currentFaceInStack->neighbor(0)->info().hasOneTag() && !facesInRegion.count(currentFaceInStack->neighbor(0)) &&
            !triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(currentFaceInStack, 0))) {
					facesInRegion.insert(currentFaceInStack->neighbor(0));
					facesToProcess.push(currentFaceInStack->neighbor(0));
				} if (!currentFaceInStack->neighbor(1)->info().hasOneTag() && !facesInRegion.count(currentFaceInStack->neighbor(1)) &&
              !triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(currentFaceInStack, 1))) {
					facesInRegion.insert(currentFaceInStack->neighbor(1));
					facesToProcess.push(currentFaceInStack->neighbor(1));
				} if (!currentFaceInStack->neighbor(2)->info().hasOneTag() && !facesInRegion.count(currentFaceInStack->neighbor(2)) &&
              !triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(currentFaceInStack, 2))) {
					facesInRegion.insert(currentFaceInStack->neighbor(2));
					facesToProcess.push(currentFaceInStack->neighbor(2));
				}
			}
			// Add up the boundary for each triangle and tag
			for (std::set<Triangulation::Face_handle>::iterator currentFaceInRegion = facesInRegion.begin(); currentFaceInRegion != facesInRegion.end(); ++currentFaceInRegion) {
				if (!facesInRegion.count((*currentFaceInRegion)->neighbor(0))) {
					addToLength(tagBoundaryLength, (*currentFaceInRegion)->neighbor(0)->info().getTags(), sqrt(CGAL::to_double(triangulation.segment(*currentFaceInRegion, 0).squared_length())));
				} if (!facesInRegion.count((*currentFaceInRegion)->neighbor(1))) {
					addToLength(tagBoundaryLength, (*currentFaceInRegion)->neighbor(1)->info().getTags(), sqrt(CGAL::to_double(triangulation.segment(*currentFaceInRegion, 1).squared_length())));
				} if (!facesInRegion.count((*currentFaceInRegion)->neighbor(2))) {
					addToLength(tagBoundaryLength, (*currentFaceInRegion)->neighbor(2)->info().getTags(), sqrt(CGAL::to_double(triangulation.segment(*currentFaceInRegion, 2).squared_length())));
				}
			}
			// Find the tag with longest boundary
			double maxLength = 0.0;
			std::map<PolygonHandle *, double>::iterator longest = tagBoundaryLength.end();
			for (std::map<PolygonHandle *, double>::iterator currentLength = tagBoundaryLength.begin(); currentLength != tagBoundaryLength.end(); ++currentLength) {
				if (currentLength->first != NULL && (alsoUniverse || currentLength->first != &universetag)) {
					if (currentLength->second > maxLength && (currentFace->info().hasTag(currentLength->first) || currentFace->info().hasNoTags())) {
						maxLength = currentLength->second;
						longest = currentLength;
					} else if (currentLength->second == maxLength) {
						longest = tagBoundaryLength.end();
					}
				}
			}
			// Assign the region to the tag with the longest boundary (if there is one)
			if (longest == tagBoundaryLength.end()) repaired = false;
			else {
				for (std::set<Triangulation::Face_handle>::iterator currentFaceInRegion = facesInRegion.begin(); currentFaceInRegion != facesInRegion.end(); ++currentFaceInRegion) {
					facesToRepair.push_back(std::pair<Triangulation::Face_handle, PolygonHandle *>(*currentFaceInRegion, longest->first));
				}
			}
		}
	}
	
	// Re-tag faces in the vector
	for (std::vector<std::pair<Triangulation::Face_handle, PolygonHandle *> >::iterator currentFace = facesToRepair.begin(); currentFace != facesToRepair.end(); ++currentFace) {
		currentFace->first->info().removeAllTags();
		currentFace->first->info().addTag(currentFace->second);
	}
	return repaired;
}


void PlanarPartition::addToLength(std::map<PolygonHandle *, double> &lengths, PolygonHandle *ph, double length) {
	if (ph == NULL)
    return;
	if (ph->isMultiPolygonHandle()) {
		MultiPolygonHandle *mph = static_cast<MultiPolygonHandle *>(ph);
		for (std::list<PolygonHandle *>::const_iterator currentPolygonHandle = mph->getHandles()->begin(); currentPolygonHandle != mph->getHandles()->end(); ++currentPolygonHandle) {
			if (lengths.count(*currentPolygonHandle))
        lengths[*currentPolygonHandle] += length;
			else
        lengths[*currentPolygonHandle] = length;
		}
	}
  else {
		if (lengths.count(ph))
      lengths[ph] += length;
		else
      lengths[ph] = length;
	}
}


bool PlanarPartition::makeAllHolesValid() {
  for (Triangulation::Finite_faces_iterator currentFace = triangulation.finite_faces_begin(); currentFace != triangulation.finite_faces_end(); ++currentFace) {
    if (currentFace->info().hasNoTags()) {
      currentFace->info().addTag(&universetag);
    }
  }
  return true;
}

bool PlanarPartition::isValid() {
  if (state < TAGGED) {
    std::cout << "Triangulation not yet tagged. Cannot check!" << std::endl;
    return false;
  } 
  if (state >= REPAIRED) {
    return true;
  }
  for (Triangulation::Finite_faces_iterator f = triangulation.finite_faces_begin(); f != triangulation.finite_faces_end(); ++f) {
    if (!(*f).info().hasOneTag()) {
      return false;  
    }
  }
  state = REPAIRED;
  return true;
}


bool PlanarPartition::add_extra_constraints_in_gaps(double splitregions) {
  std::cout << "Adding extra constraints (max length = " << splitregions << ") in gaps to improve edge-matching results...";
//-- 1. find and store split vertices with d>2 (incident to hole; incident to 1 dataset)
  std::set<Triangulation::Vertex_handle> vs_hole;
  Triangulation::Finite_vertices_iterator v = triangulation.finite_vertices_begin();
  while (v != triangulation.finite_vertices_end()) {
    Triangulation::Face_circulator ff = triangulation.incident_faces(v), curF = ff;
    int noconstraints = 0;
    bool holeinstar = false;
    bool infinite = false;
    do {
      int i = curF->index(v);
      if (curF->is_constrained(curF->ccw(i)) == true)
        noconstraints++;
      if (curF->info().isHole() == true)
        holeinstar = true;
      if (triangulation.is_infinite(curF) == true)
        infinite = true;
      curF++;
    } while (curF != ff);
    if ( (infinite == false) && (noconstraints > 2) && (holeinstar == true) ) {
        vs_hole.insert(v);
    }
    v++;
  }
//  int j = 1;
//  std::cout << "---SUMMARY hole-vertices---" << std::endl;
//  for (std::set<Triangulation::Vertex_handle>::iterator curv = vs_hole.begin(); curv != vs_hole.end(); curv++) {
//    std::cout << j << " -- " << (*curv)->point().x() << ":" << (*curv)->point().y() << std::endl;
//    j++;
//  }
  
//-- 2. add the constraint in hole being the shortest
//  std::cout << "CONSTRAINTS ADDED" << std::endl;
  Triangulation::Face_handle thef;
  bool found;
  for (std::set<Triangulation::Vertex_handle>::iterator curv = vs_hole.begin(); curv != vs_hole.end(); curv++) {
    Triangulation::Face_circulator ff = triangulation.incident_faces(*curv), curF = ff;
    //-- collect all edges
    K::FT edgelength = 1e6;
    found = false;
    //-- find the shortest edge with splitvertex
    do {
      if (curF->info().isHole() == true) {
        int i = curF->index(*curv);
        Triangulation::Vertex_handle vend = curF->vertex(curF->cw(i));
        if (curF->neighbor(curF->ccw(i))->info().isHole() == true) {
          if (vs_hole.count(vend) == 1) {
            K::FT d = CGAL::squared_distance((*curv)->point(), vend->point());
            if ( (d < edgelength) && (d < (splitregions*splitregions)) ) {
                edgelength = d;
                thef = curF;
                found = true;
            }
          }
        }
      }
      curF++;
    } while (curF != ff);
    if (found == true) {
      int i = thef->index(*curv);
      triangulation.insert_constraint(*curv, thef->vertex(thef->cw(i)));
//      std::cout << "CONSTRAINED EDGE (SPLITVERTEX) ADDED:" << std::endl;
//      std::cout << " A: " << (*curv)->point().x() << ":" << (*curv)->point().y() << std::endl;
//      std::cout << " B: " << (thef->vertex(thef->cw(i)))->point().x() << ":" << (thef->vertex(thef->cw(i)))->point().y() << std::endl;
    }
    else {
      edgelength = 1e6;
      found = false;
      //-- find the shortest edge WITHOUT splitvertex
      do {
        if (curF->info().isHole() == true) {
          int i = curF->index(*curv);
          Triangulation::Vertex_handle vend = curF->vertex(curF->cw(i));
          if (curF->neighbor(curF->ccw(i))->info().isHole() == true) {
            K::FT d = CGAL::squared_distance((*curv)->point(), vend->point());
            if ( (d < edgelength) && (d < (splitregions*splitregions)) ) {
              edgelength = d;
              thef = curF;
              found = true;
            }
          }
        }
        curF++;
      } while (curF != ff);
      if (found == true) {
        int i = thef->index(*curv);
        triangulation.insert_constraint(*curv, thef->vertex(thef->cw(i)));
//        std::cout << "CONSTRAINED EDGE ADDED:" << std::endl;
//        std::cout << " A: " << (*curv)->point().x() << ":" << (*curv)->point().y() << std::endl;
//        std::cout << " B: " << (thef->vertex(thef->cw(i)))->point().x() << ":" << (thef->vertex(thef->cw(i)))->point().y() << std::endl;
      }
    }
  }
  std::cout << "done" << std::endl;
  return true;
}


bool PlanarPartition::splitRegions(double ratio) {
	if (state < TAGGED) {
		std::cout << "Triangulation not yet tagged. Cannot split!" << std::endl;
		return false;
	} if (state > TAGGED) {
		std::cout << "Triangulation already repaired!" << std::endl;
		return false;
	}
	std::cout << "Splitting regions..." << std::endl;
	time_t thisTime = time(NULL);
	
	if (ratio <= 1.0)
    return false;
	
	double shortSide, longSide, thisSide;
	unsigned int whichSide, splits = 0;
	for (Triangulation::Finite_faces_iterator currentFace = triangulation.finite_faces_begin(); currentFace != triangulation.finite_faces_end(); ++currentFace) {
		// Check for the longest and shortest sides
		shortSide = longSide = sqrt(CGAL::to_double(triangulation.segment(currentFace, 0).squared_length()));
		whichSide = 0;
		thisSide = sqrt(CGAL::to_double(triangulation.segment(currentFace, 1).squared_length()));
		if (thisSide > longSide)
      longSide = thisSide;
		else if (thisSide < shortSide) {
			shortSide = thisSide;
			whichSide = 1;
		}
    thisSide = sqrt(CGAL::to_double(triangulation.segment(currentFace, 2).squared_length()));
		if (thisSide > longSide)
      longSide = thisSide;
		else if (thisSide < shortSide) {
			shortSide = thisSide;
			whichSide = 2;
		}
		// Add constrained edge if they exceed the long/short ratio
		if (longSide/shortSide >= ratio) {
			currentFace->set_constraint(whichSide, true);
			++splits;
		}
	}
	std::cout << "\t" << splits << " constrained edges added." << std::endl;
	std::cout << "\tRegions split (" << time(NULL)-thisTime << " s)." << std::endl;
	return true;
}



bool PlanarPartition::reconstructPolygons(bool removeVertices) {
  if (state < TAGGED) {
		std::cout << "Triangulation not tagged. Cannot reconstruct!" << std::endl;
		return false;
	} if (state < REPAIRED) std::cout << "Warning: Triangulation not yet repaired. There could be errors..." << std::endl;
	if (state > REPAIRED) {
		std::cout << "Polygons already reconstructed!" << std::endl;
		return false;
	}
	std::cout << "Reconstructing polygons (geometry)..." << std::endl;
	time_t thisTime = time(NULL);
  
  this->removeConstraints();
  if (removeVertices)
    this->removeVertices();
  
  for (Triangulation::Finite_faces_iterator seedingFace = triangulation.finite_faces_begin(); seedingFace != triangulation.finite_faces_end(); ++seedingFace) {
    PolygonHandle *currentTag = seedingFace->info().getOneTag();
    if (currentTag == NULL) continue;
    
    // STEP 1: Find a suitable seeding triangle (connected to the outer boundary)
    if (currentTag == &universetag) {
      seedingFace->info().removeAllTags();
      continue;
    } if (seedingFace->neighbor(0)->info().getOneTag() == currentTag &&
          seedingFace->neighbor(1)->info().getOneTag() == currentTag &&
          seedingFace->neighbor(2)->info().getOneTag() == currentTag) continue;
    
    // STEP 2: Get boundary
    seedingFace->info().removeAllTags();
    std::list<Triangulation::Vertex_handle> vertices;
    if (seedingFace->neighbor(2)->info().hasTag(currentTag)) {
      seedingFace->neighbor(2)->info().removeAllTags();
      std::list<Triangulation::Vertex_handle> *l2 = getBoundary(seedingFace->neighbor(2), seedingFace->neighbor(2)->index(seedingFace), currentTag);
      vertices.splice(vertices.end(), *l2);
      delete l2;
    } vertices.push_back(seedingFace->vertex(0));
    if (seedingFace->neighbor(1)->info().hasTag(currentTag)) {
      seedingFace->neighbor(1)->info().removeAllTags();
      std::list<Triangulation::Vertex_handle> *l1 = getBoundary(seedingFace->neighbor(1), seedingFace->neighbor(1)->index(seedingFace), currentTag);
      vertices.splice(vertices.end(), *l1);
      delete l1;
    } vertices.push_back(seedingFace->vertex(2));
    if (seedingFace->neighbor(0)->info().hasTag(currentTag)) {
      seedingFace->neighbor(0)->info().removeAllTags();
      std::list<Triangulation::Vertex_handle> *l0 = getBoundary(seedingFace->neighbor(0), seedingFace->neighbor(0)->index(seedingFace), currentTag);
      vertices.splice(vertices.end(), *l0);
      delete l0;
    } vertices.push_back(seedingFace->vertex(1));
    
    // STEP 3: Find cutting vertices
    std::set<Triangulation::Vertex_handle> visitedVertices;
    std::set<Triangulation::Vertex_handle> repeatedVertices;
    for (std::list<Triangulation::Vertex_handle>::iterator currentVertex = vertices.begin(); currentVertex != vertices.end(); ++currentVertex) {
      if (!visitedVertices.insert(*currentVertex).second) repeatedVertices.insert(*currentVertex);
    } visitedVertices.clear();
    
    // STEP 4: Cut and join rings in the correct order
    std::list<std::list<Triangulation::Vertex_handle> *> rings;
    std::stack<std::list<Triangulation::Vertex_handle> *> chainsStack;
    std::map<Triangulation::Vertex_handle, std::list<Triangulation::Vertex_handle> *> vertexChainMap;
    std::list<Triangulation::Vertex_handle> *newChain = new std::list<Triangulation::Vertex_handle>();
    
    // New vertex
    for (std::list<Triangulation::Vertex_handle>::iterator currentVertex = vertices.begin(); currentVertex != vertices.end(); ++currentVertex) {
      // New chain
      if (repeatedVertices.count(*currentVertex) > 0) {
        // Closed by itself
        if (newChain->front() == *currentVertex) {
          // Degenerate (insufficient vertices to be valid)
          if (newChain->size() < 3) delete newChain;
          else {
            std::list<Triangulation::Vertex_handle>::iterator secondElement = newChain->begin();
            ++secondElement;
            // Degenerate (zero area)
            if (newChain->back() == *secondElement) delete newChain;
            // Valid
            else rings.push_back(newChain);
          }
        }
        // Open by itself
        else {
          // Closed with others in stack
          if (vertexChainMap.count(*currentVertex)) {
            while (chainsStack.top() != vertexChainMap[*currentVertex]) {
              newChain->splice(newChain->begin(), *chainsStack.top());
              chainsStack.pop();
            } newChain->splice(newChain->begin(), *chainsStack.top());
            chainsStack.pop();
            vertexChainMap.erase(*currentVertex);
            // Degenerate (insufficient vertices to be valid)
            if (newChain->size() < 3) delete newChain;
            else {
              std::list<Triangulation::Vertex_handle>::iterator secondElement = newChain->begin();
              ++secondElement;
              // Degenerate (zero area)
              if (newChain->back() == *secondElement) delete newChain;
              // Valid
              else rings.push_back(newChain);
            }
          }
          // Open
          else {
            // Not first chain
            if (repeatedVertices.count(newChain->front()) > 0)
              vertexChainMap[newChain->front()] = newChain;
            chainsStack.push(newChain);
          }
        }
        newChain = new std::list<Triangulation::Vertex_handle>();
      }
      newChain->push_back(*currentVertex);
    }
    
    // Final ring
    while (chainsStack.size() > 0) {
      newChain->splice(newChain->begin(), *chainsStack.top());
      chainsStack.pop();
    }
    
    // Degenerate (insufficient vertices to be valid)
    if (newChain->size() < 3) {
      delete newChain;
    }
    else {
      std::list<Triangulation::Vertex_handle>::iterator secondElement = newChain->begin();
      ++secondElement;
      // Degenerate (zero area)
      if (newChain->back() == *secondElement)
        delete newChain;
      // Valid
      else
        rings.push_back(newChain);
    }
    
    if (chainsStack.size() > 0)
      std::cout << "Error: Stack has " << chainsStack.size() << " elements. Should be empty." << std::endl;
    
    // STEP 5: Make a polygon from this list and save it
    std::vector<Ring> innerRings;
    Ring outerRing;
    for (std::list<std::list<Triangulation::Vertex_handle> *>::iterator currentRing = rings.begin(); currentRing != rings.end(); ++currentRing) {
      Ring newRing;
      for (std::list<Triangulation::Vertex_handle>::iterator currentPoint = (*currentRing)->begin(); currentPoint != (*currentRing)->end(); ++currentPoint) {
        newRing.push_back((*currentPoint)->point());
      }
      if (newRing.is_clockwise_oriented())
        outerRing = newRing;
      else
        innerRings.push_back(newRing);
    }
    outputPolygons.push_back(std::pair<PolygonHandle *, Polygon>(currentTag, Polygon(outerRing, innerRings.begin(), innerRings.end())));
    // Free memory from the chains
    for (std::list<std::list<Triangulation::Vertex_handle> *>::iterator currentRing = rings.begin(); currentRing != rings.end(); ++currentRing) {
      delete *currentRing;
    }
  }
  state = RECONSTRUCTED;
	std::cout << "Polygons reconstructed (" << time(NULL)-thisTime << " s)." << std::endl;
  return true;
}


bool PlanarPartition::exportPolygonsSHP(std::string &folder) {
	if (state < RECONSTRUCTED) {
		std::cout << "Polygons have not been reconstructed yet. Nothing to export!" << std::endl;
		return false;
	}
	
	std::cout << "Exporting polygons..." << std::endl;
	time_t thisTime = time(NULL);
  
//-- 1. get all the diff OGRFeatureDefn
  std::set<OGRFeatureDefn*> allFDefs;
  for (std::vector<std::pair<PolygonHandle *, Polygon> >::iterator it = outputPolygons.begin();
       it != outputPolygons.end();
       ++it) {
    allFDefs.insert(it->first->feature->GetDefnRef());
  }

//-- 2. create new a SHP for each one
  #if GDAL_VERSION_MAJOR < 2
    std::map<OGRFeatureDefn*,OGRDataSource*> allshps;
    std::map<OGRFeatureDefn*,OGRDataSource*>::iterator allshpsit;
  #else
    std::map<OGRFeatureDefn*,GDALDataset*> allshps;
    std::map<OGRFeatureDefn*,GDALDataset*>::iterator allshpsit;
  #endif
  for (std::set<OGRFeatureDefn*>::iterator it = allFDefs.begin(); it != allFDefs.end(); ++it) {
    const char *driverName = "ESRI Shapefile";
    #if GDAL_VERSION_MAJOR < 2
  	  OGRSFDriver *driver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(driverName);
    #else
 	    GDALDriver *driver = GetGDALDriverManager()->GetDriverByName(driverName);
    #endif
    if (driver == NULL) {
      std::cout << "\tError: OGR Shapefile driver not found." << std::endl;
      return false;
    }
    std::string tmp = (*it)->GetName();
    std::string outname = folder + "/" + tmp + ".r.shp";
    #if GDAL_VERSION_MAJOR < 2
  	  OGRDataSource *dataSource = driver->Open(outname.c_str(), false);
    #else
      GDALDataset *dataSource = (GDALDataset*) GDALOpenEx(outname.c_str(), GDAL_OF_READONLY, NULL, NULL, NULL);
    #endif
    if (dataSource != NULL) {
      std::cout << "\tOverwriting file..." << std::endl;
      #if GDAL_VERSION_MAJOR < 2
        if (driver->DeleteDataSource(dataSource->GetName())!= OGRERR_NONE) {
      #else
        if (driver->Delete(outname.c_str())!= CE_None) {
      #endif
          std::cout << "\tError: Couldn't erase file with same name." << std::endl;
          return false;
      }
      #if GDAL_VERSION_MAJOR < 2
        OGRDataSource::DestroyDataSource(dataSource);
      #else
      	GDALClose(dataSource);
      #endif
    }
    std::cout << "\tWriting file " << outname << std::endl;
    #if GDAL_VERSION_MAJOR < 2
  	  dataSource = driver->CreateDataSource(outname.c_str(), NULL);
    #else
    	dataSource = driver->Create(outname.c_str(),0,0,0,GDT_Unknown,NULL);
    #endif
    if (dataSource == NULL) {
      std::cout << "\tError: Could not create file." << std::endl;
      return false;
    }
    allshps[*it] = dataSource;
    // TODO : deal with CRS
    //	OGRLayer *layer = dataSource->CreateLayer("polygons", spatialReference, wkbPolygon, NULL);
    OGRLayer *layer = dataSource->CreateLayer("polygons", NULL, wkbPolygon, NULL);
    if (layer == NULL) {
      std::cout << "\tError: Could not create layer." << std::endl;
      return false;
    }
    //-- create all the fields
    for (int i = 0; i < (*it)->GetFieldCount(); i++) {
      layer->CreateField((*it)->GetFieldDefn(i));
    }
  }

  for (std::vector<std::pair<PolygonHandle *, Polygon> >::iterator currentPolygon = outputPolygons.begin();
       currentPolygon != outputPolygons.end();
       ++currentPolygon) {
		OGRPolygon polygon;
		OGRLinearRing outerRing;
    OGRFeature* f = currentPolygon->first->feature;
    allshpsit = allshps.find(f->GetDefnRef());
    OGRLayer *layer = allshpsit->second->GetLayer(0);
  
    if (currentPolygon->second.outer_boundary().size() < 1)
      continue;
		for (Ring::Vertex_iterator currentVertex = currentPolygon->second.outer_boundary().vertices_begin();
         currentVertex != currentPolygon->second.outer_boundary().vertices_end();
         ++currentVertex) {
			outerRing.addPoint(CGAL::to_double(currentVertex->x()), CGAL::to_double(currentVertex->y()));
		}
    outerRing.addPoint(CGAL::to_double(currentPolygon->second.outer_boundary().vertex(0).x()),
                       CGAL::to_double(currentPolygon->second.outer_boundary().vertex(0).y()));
		polygon.addRing(&outerRing);
		for (Polygon::Hole_const_iterator currentRing = currentPolygon->second.holes_begin(); currentRing != currentPolygon->second.holes_end(); ++currentRing) {
			OGRLinearRing innerRing;
			for (Ring::Vertex_iterator currentVertex = currentRing->vertices_begin(); currentVertex != currentRing->vertices_end(); ++currentVertex) {
				innerRing.addPoint(CGAL::to_double(currentVertex->x()), CGAL::to_double(currentVertex->y()));
			}
      innerRing.addPoint(CGAL::to_double(currentRing->vertex(0).x()), CGAL::to_double(currentRing->vertex(0).y()));
			polygon.addRing(&innerRing);
		}
		
		f->SetGeometry(&polygon);
		// Create OGR feature
		if (layer->CreateFeature(f) != OGRERR_NONE)
      std::cout << "\tError: Could not create feature." << std::endl;
//		OGRFeature::DestroyFeature(f); //-- TODO: free the features? hmmmm...
	}
  //-- clear memory for all the created SHP
  for (allshpsit = allshps.begin(); allshpsit != allshps.end(); ++allshpsit) {
    #if GDAL_VERSION_MAJOR < 2
  	  OGRDataSource::DestroyDataSource(allshpsit->second);
    #else
 	    GDALClose(allshpsit->second);
    #endif
  }
	std::cout << "Polygons exported (" << time(NULL)-thisTime << " s)." << std::endl;
  return true;
}


bool PlanarPartition::exportTriangulation(std::string &file) {
	if (state < TRIANGULATED || state > REPAIRED) {
		std::cout << "No triangulation to export!" << std::endl;
		return false;
	}
  
  std::cout << "Exporting triangulation as a SHP..." << std::endl;
  
  const char *driverName = "ESRI Shapefile";
  #if GDAL_VERSION_MAJOR < 2
  	OGRSFDriver *driver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(driverName);
  #else
  	GDALDriver *driver = GetGDALDriverManager()->GetDriverByName(driverName);
  #endif
	if (driver == NULL) {
		std::cout << "Driver not found." << std::endl;
		return false;
	}
	#if GDAL_VERSION_MAJOR < 2
  	OGRDataSource *dataSource = driver->Open(file.c_str(), false);
 #else
 	  GDALDataset *dataSource = (GDALDataset*) GDALOpenEx(file.c_str(), GDAL_OF_READONLY, NULL, NULL, NULL);
 #endif
	if (dataSource != NULL) {
		std::cout << "Erasing current file..." << std::endl;
		#if GDAL_VERSION_MAJOR < 2
  		if (driver->DeleteDataSource(dataSource->GetName())!= OGRERR_NONE) {
    #else
    	if (driver->Delete(file.c_str())!= CE_None) {
    #endif
        std::cout << "Couldn't erase current file." << std::endl;
        return false;
		}
    #if GDAL_VERSION_MAJOR < 2
      OGRDataSource::DestroyDataSource(dataSource);
    #else
    	GDALClose(dataSource);
    #endif
	}
  std::cout << "Writing file " << file << "..." << std::endl;
  #if GDAL_VERSION_MAJOR < 2
    dataSource = driver->CreateDataSource(file.c_str(), NULL);
  #else
    dataSource = driver->Create(file.c_str(),0,0,0,GDT_Unknown,NULL);
  #endif
	if (dataSource == NULL) {
		std::cout << "Could not create file." << std::endl;
		return false;
	}
	OGRLayer *layer = dataSource->CreateLayer("triangles", NULL, wkbPolygon, NULL);
	if (layer == NULL) {
		std::cout << "Could not create layer." << std::endl;
		return false;
	}
  OGRFieldDefn numberOfTagsField("Tags", OFTInteger);
  if (layer->CreateField(&numberOfTagsField) != OGRERR_NONE) {
    std::cout << "Could not create field Tags." << std::endl;
    return false;
	}
  
  for (CDT::Finite_faces_iterator currentFace = triangulation.finite_faces_begin();
       currentFace != triangulation.finite_faces_end();
       ++currentFace) {
		OGRLinearRing ring;
		ring.addPoint(CGAL::to_double((*(*currentFace).vertex(0)).point().x()), CGAL::to_double((*(*currentFace).vertex(0)).point().y()), 0.0);
		ring.addPoint(CGAL::to_double((*(*currentFace).vertex(1)).point().x()), CGAL::to_double((*(*currentFace).vertex(1)).point().y()), 0.0);
		ring.addPoint(CGAL::to_double((*(*currentFace).vertex(2)).point().x()), CGAL::to_double((*(*currentFace).vertex(2)).point().y()), 0.0);
		ring.addPoint(CGAL::to_double((*(*currentFace).vertex(0)).point().x()), CGAL::to_double((*(*currentFace).vertex(0)).point().y()), 0.0);
		OGRPolygon polygon;
		polygon.addRing(&ring);
		OGRFeature *feature = OGRFeature::CreateFeature(layer->GetLayerDefn());
    if ((*currentFace).info().getTags() == NULL) {
      feature->SetField("Tags", 0);
    }
    else if ((*currentFace).info().getTags() != &universetag) {
      feature->SetField("Tags", (int)(*currentFace).info().numberOfTags());
    }
    else {
      feature->SetField("Tags", 0);
    }
  	feature->SetGeometry(&polygon);
		if (layer->CreateFeature(feature) != OGRERR_NONE)
      std::cout << "Could not create feature." << std::endl;
		OGRFeature::DestroyFeature(feature);
	}
  #if GDAL_VERSION_MAJOR < 2
    OGRDataSource::DestroyDataSource(dataSource);
  #else
    GDALClose(dataSource);
  #endif
	return true;
}


void PlanarPartition::getProblemRegionsAsOGR(std::vector<OGRGeometry*> &holes, std::vector<OGRGeometry*> &overlaps)
{
  std::vector<std::pair<Triangulation::Face_handle, PolygonHandle *> > facesToRepair;
  std::set<Triangulation::Face_handle> processedFaces;
  for (Triangulation::Finite_faces_iterator currentFace = triangulation.finite_faces_begin(); currentFace != triangulation.finite_faces_end(); ++currentFace) {
    if (!currentFace->info().hasOneTag() && !processedFaces.count(currentFace)) {
      // Expand this triangle into a complete region
      std::set<Triangulation::Face_handle> facesInRegion;
      facesInRegion.insert(currentFace);
      std::stack<Triangulation::Face_handle> facesToProcess;
      facesToProcess.push(currentFace);
      while (facesToProcess.size() > 0) {
        Triangulation::Face_handle currentFaceInStack = facesToProcess.top();
        facesToProcess.pop();
        processedFaces.insert(currentFaceInStack);
        if (!currentFaceInStack->neighbor(0)->info().hasOneTag() && !facesInRegion.count(currentFaceInStack->neighbor(0)) &&
            !triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(currentFaceInStack, 0))) {
          facesInRegion.insert(currentFaceInStack->neighbor(0));
          facesToProcess.push(currentFaceInStack->neighbor(0));
        } if (!currentFaceInStack->neighbor(1)->info().hasOneTag() && !facesInRegion.count(currentFaceInStack->neighbor(1)) &&
              !triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(currentFaceInStack, 1))) {
          facesInRegion.insert(currentFaceInStack->neighbor(1));
          facesToProcess.push(currentFaceInStack->neighbor(1));
        } if (!currentFaceInStack->neighbor(2)->info().hasOneTag() && !facesInRegion.count(currentFaceInStack->neighbor(2)) &&
              !triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(currentFaceInStack, 2))) {
          facesInRegion.insert(currentFaceInStack->neighbor(2));
          facesToProcess.push(currentFaceInStack->neighbor(2));
        }
      }
      //-- construct OGR polygon
      OGRMultiPolygon settriangles;
      OGRPolygon polygon;
      OGRLinearRing oring;
      for (std::set<Triangulation::Face_handle>::iterator curF = facesInRegion.begin(); curF != facesInRegion.end(); ++curF) {
        oring.addPoint(CGAL::to_double((*curF)->vertex(0)->point().x()), CGAL::to_double((*curF)->vertex(0)-> point().y()));
        oring.addPoint(CGAL::to_double((*curF)->vertex(1)->point().x()), CGAL::to_double((*curF)->vertex(1)-> point().y()));
        oring.addPoint(CGAL::to_double((*curF)->vertex(2)->point().x()), CGAL::to_double((*curF)->vertex(2)-> point().y()));
        oring.addPoint(CGAL::to_double((*curF)->vertex(0)->point().x()), CGAL::to_double((*curF)->vertex(0)-> point().y()));
        polygon.addRing(&oring);
        settriangles.addGeometryDirectly(polygon.clone());
        polygon.empty();
        oring.empty();
      }
      OGRGeometry* u;
      u = settriangles.UnionCascaded();
      if (currentFace->info().numberOfTags() == 0)
        holes.push_back(u->clone());
      else
        overlaps.push_back(u->clone());
    }
  }
}

void PlanarPartition::getListOverlappingPolygons(std::ostream &ostr) {
  std::map<int, std::set<int> > overlaps;
  for (Triangulation::Finite_faces_iterator curF = triangulation.finite_faces_begin(); curF != triangulation.finite_faces_end(); ++curF) {
    if ((*curF).info().isOverlap()) {
      std::list<int> tmp;
      MultiPolygonHandle *handle = static_cast<MultiPolygonHandle *>((*curF).info().getTags());
      for (std::list<PolygonHandle *>::const_iterator curtag = handle->getHandles()->begin(); curtag != handle->getHandles()->end(); ++curtag) {
        if (*curtag == &universetag)
          continue;
        int fid = (*curtag)->getFID();
        tmp.push_back(fid);
      }
      for (std::list<int>::iterator i = tmp.begin(); i != tmp.end(); i++) {
        std::map<int, std::set<int> >::iterator search = overlaps.find(*i);
        if (search != overlaps.end()) {
          (search->second).insert(tmp.begin(), tmp.end());
        }
        else {
          std::set<int> t(tmp.begin(), tmp.end());
          overlaps[*i] = t;
        }
      }
    }
  }
  for (std::map<int, std::set<int> >::iterator i = overlaps.begin(); i != overlaps.end(); i++) {
    ostr << i->first << ",";
    for (std::set<int>::iterator j = (i->second).begin(); j != (i->second).end(); j++) {
      if (*j != i->first)
        ostr << *j << ",";
    }
    ostr << std::endl;
  }
}

void PlanarPartition::printProblemRegions(std::ostream &ostr) {
  std::vector<std::pair<Triangulation::Face_handle, PolygonHandle *> > facesToRepair;
  std::set<Triangulation::Face_handle> processedFaces;
  unsigned int holes = 0;
  unsigned int overlaps = 0;
  for (Triangulation::Finite_faces_iterator currentFace = triangulation.finite_faces_begin(); currentFace != triangulation.finite_faces_end(); ++currentFace) {
    if (!currentFace->info().hasOneTag() && !processedFaces.count(currentFace)) {
      // Expand this triangle into a complete region
      std::set<Triangulation::Face_handle> facesInRegion;
      facesInRegion.insert(currentFace);
      std::stack<Triangulation::Face_handle> facesToProcess;
      facesToProcess.push(currentFace);
      while (facesToProcess.size() > 0) {
        Triangulation::Face_handle currentFaceInStack = facesToProcess.top();
        facesToProcess.pop();
        processedFaces.insert(currentFaceInStack);
        if (!currentFaceInStack->neighbor(0)->info().hasOneTag() && !facesInRegion.count(currentFaceInStack->neighbor(0)) &&
            !triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(currentFaceInStack, 0))) {
          facesInRegion.insert(currentFaceInStack->neighbor(0));
          facesToProcess.push(currentFaceInStack->neighbor(0));
        } if (!currentFaceInStack->neighbor(1)->info().hasOneTag() && !facesInRegion.count(currentFaceInStack->neighbor(1)) &&
              !triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(currentFaceInStack, 1))) {
          facesInRegion.insert(currentFaceInStack->neighbor(1));
          facesToProcess.push(currentFaceInStack->neighbor(1));
        } if (!currentFaceInStack->neighbor(2)->info().hasOneTag() && !facesInRegion.count(currentFaceInStack->neighbor(2)) &&
              !triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(currentFaceInStack, 2))) {
          facesInRegion.insert(currentFaceInStack->neighbor(2));
          facesToProcess.push(currentFaceInStack->neighbor(2));
        }
      }
      if (currentFace->info().numberOfTags() == 0)
        holes++;
      else
        overlaps++;
    }
  }
  ostr << "*** Problematic Regions ***" << std::endl <<
  "\tOverlaps: "         << overlaps << " regions" << std::endl <<
  "\tHoles:    "         << holes << " regions" << std::endl;
}


void PlanarPartition::printTriangulationInfo(std::ostream &ostr) {
  unsigned int tag0 = 0;
  unsigned int tag1 = 0;
  unsigned int tag2 = 0;
  unsigned int total;
  double tag0_a = 0; //-- areas
  double tag1_a = 0;
  double tag2_a = 0;
  for (Triangulation::Finite_faces_iterator curF = triangulation.finite_faces_begin(); curF != triangulation.finite_faces_end(); ++curF) {
    if ((*curF).info().hasNoTags()) {
      tag0++;
      tag0_a += CGAL::to_double(triangulation.triangle(curF).area());
    }
    else if ((*curF).info().hasOneTag()) {
      tag1++;
      tag1_a += CGAL::to_double(triangulation.triangle(curF).area());
    }
    else {
      tag2++;
      tag2_a += CGAL::to_double(triangulation.triangle(curF).area());
    }
  }
  
  ostr << "*** Triangulation ***" << std::endl <<
  "\tVertices: " << triangulation.number_of_vertices() << std::endl <<
  "\tEdges: " << triangulation.tds().number_of_edges() << std::endl <<
  "\tTriangles: " << triangulation.number_of_faces() << std::endl;
  total = tag0 + tag1 + tag2;
  ostr << "\tOk:       " << tag1 << " triangles" << std::endl <<
  "\tOverlaps: "         << tag2 << " triangles" << std::endl <<
  "\tHoles:    "         << tag0 << " triangles" << std::endl;
}


void PlanarPartition::reportProblemRegions(std::ostream &ostr, double thinness, double minSliverArea) {
  std::vector<OGRGeometry*> holes;
  std::vector<OGRGeometry*> overlaps;
  getProblemRegionsAsOGR(holes, overlaps);
  ostr << "*** Regions ***" << std::endl;
  //-- overlaps
  ostr << "\tOverlaps: " << overlaps.size()  << " regions(s)." << std::endl;
  for (std::vector<OGRGeometry*>::iterator g = overlaps.begin() ; g != overlaps.end(); g++) {
    OGRPolygon *tmp = static_cast<OGRPolygon*>(*g);
    std::cout << "\t\t- " << tmp->get_Area() << " unit^2" << std::endl;
  }
  
  //-- holes
  if (thinness < 0.0) {
    ostr << "\tHoles: " << holes.size()  << " regions(s)." << std::endl;
    for (std::vector<OGRGeometry*>::iterator g = holes.begin() ; g != holes.end(); g++) {
      OGRPolygon *tmp = static_cast<OGRPolygon*>(*g);
      std::cout << "\t\t- " << tmp->get_Area() << " unit^2" << std::endl;
    }
  }
  else {
    std::vector<OGRGeometry*> slivers;
    for (std::vector<OGRGeometry*>::iterator g = holes.begin() ; g != holes.end(); g++) {
      OGRPolygon *tmp = static_cast<OGRPolygon*>(*g);
      //-- use the magic formula from ELF project
      double thinness = 4 * CGAL_PI * tmp->get_Area() / pow(tmp->getExteriorRing()->get_Length(), 2);
      if ( (thinness < 0.30) && (tmp->get_Area() < minSliverArea) ) {
        slivers.push_back(*g);
      }
    }
    ostr << "\tHoles (valid slivers): " << slivers.size()  << " regions(s)." << std::endl;
    for (std::vector<OGRGeometry*>::iterator g = slivers.begin() ; g != slivers.end(); g++) {
      OGRPolygon *tmp = static_cast<OGRPolygon*>(*g);
      std::cout << "\t\t- " << tmp->get_Area() << " unit^2" << std::endl;
    }
  }
}


bool PlanarPartition::exportProblemRegionsAsSHP(std::string &file, double thinness, double minSliverArea) {
  std::cout << "Exporting problematic regions as a SHP..." << std::endl;
  std::vector<OGRGeometry*> holes;
  std::vector<OGRGeometry*> overlaps;
  getProblemRegionsAsOGR(holes, overlaps);
  
  const char *driverName = "ESRI Shapefile";
  #if GDAL_VERSION_MAJOR < 2
    OGRSFDriver *driver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(driverName);
  #else
  	GDALDriver *driver = GetGDALDriverManager()->GetDriverByName(driverName);
  #endif
  if (driver == NULL) {
    std::cout << "Driver not found." << std::endl;
    return false;
  }
  #if GDAL_VERSION_MAJOR < 2
  	OGRDataSource *dataSource = driver->Open(file.c_str(), false);
 #else
 	  GDALDataset *dataSource = (GDALDataset*) GDALOpenEx(file.c_str(), GDAL_OF_READONLY, NULL, NULL, NULL);
 #endif
  if (dataSource != NULL) {
    std::cout << "Erasing current file..." << std::endl;
    #if GDAL_VERSION_MAJOR < 2
      if (driver->DeleteDataSource(dataSource->GetName())!= OGRERR_NONE) {
    #else
      if (driver->Delete(file.c_str())!= CE_None) {
    #endif
        std::cout << "Couldn't erase current file." << std::endl;
        return false;
    }
    #if GDAL_VERSION_MAJOR < 2
      OGRDataSource::DestroyDataSource(dataSource);
    #else
      GDALClose(dataSource);
    #endif
  }
  std::cout << "Writing file " << file << "..." << std::endl;
  #if GDAL_VERSION_MAJOR < 2
  	  dataSource = driver->CreateDataSource(file.c_str(), NULL);
  #else
    dataSource = driver->Create(file.c_str(),0,0,0,GDT_Unknown,NULL);
  #endif
  if (dataSource == NULL) {
    std::cout << "Could not create file." << std::endl;
    return false;
  }
  OGRLayer *layer = dataSource->CreateLayer("triangles", NULL, wkbPolygon, NULL);
  if (layer == NULL) {
    std::cout << "Could not create layer." << std::endl;
    return false;
  }
  OGRFieldDefn numberOfTagsField("Tags", OFTInteger);
  if (layer->CreateField(&numberOfTagsField) != OGRERR_NONE) {
    std::cout << "Could not create field Tags." << std::endl;
    return false;
  }
  //-- overlaps
  for (std::vector<OGRGeometry*>::iterator g = overlaps.begin() ; g != overlaps.end(); g++) {
    OGRPolygon *tmp = static_cast<OGRPolygon*>(*g);
    OGRFeature *feature = OGRFeature::CreateFeature(layer->GetLayerDefn());
    feature->SetField("Tags", 1);
    feature->SetGeometry(tmp);
    if (layer->CreateFeature(feature) != OGRERR_NONE)
      std::cout << "Could not create feature." << std::endl;
    // Free OGR feature
    OGRFeature::DestroyFeature(feature);
  }
  //-- holes
  if (thinness < 0.0) {
    for (std::vector<OGRGeometry*>::iterator g = holes.begin() ; g != holes.end(); g++) {
      OGRPolygon *tmp = static_cast<OGRPolygon*>(*g);
      OGRFeature *feature = OGRFeature::CreateFeature(layer->GetLayerDefn());
      feature->SetField("Tags", 0);
      feature->SetGeometry(tmp);
      if (layer->CreateFeature(feature) != OGRERR_NONE)
        std::cout << "Could not create feature." << std::endl;
      // Free OGR feature
      OGRFeature::DestroyFeature(feature);
    }
  }
  else {
    std::cout << "(ignoring slivers having [thinness < 0.3] and [minarea < " << minSliverArea << "])" << std::endl;
    std::vector<OGRGeometry*> slivers;
    for (std::vector<OGRGeometry*>::iterator g = holes.begin() ; g != holes.end(); g++) {
      OGRPolygon *tmp = static_cast<OGRPolygon*>(*g);
      //-- use the magic formula from ELF project
      double thinness = 4 * CGAL_PI * tmp->get_Area() / pow(tmp->getExteriorRing()->get_Length(), 2);
      if ( (thinness < 0.30) && (tmp->get_Area() < minSliverArea) ) {
        slivers.push_back(*g);
      }
    }
    for (std::vector<OGRGeometry*>::iterator g = slivers.begin() ; g != slivers.end(); g++) {
      OGRPolygon *tmp = static_cast<OGRPolygon*>(*g);
      OGRFeature *feature = OGRFeature::CreateFeature(layer->GetLayerDefn());
      feature->SetField("Tags", 0);
      feature->SetGeometry(tmp);
      if (layer->CreateFeature(feature) != OGRERR_NONE)
        std::cout << "Could not create feature." << std::endl;
      // Free OGR feature
      OGRFeature::DestroyFeature(feature);
    }
  }
  #if GDAL_VERSION_MAJOR < 2
    OGRDataSource::DestroyDataSource(dataSource);
  #else
    GDALClose(dataSource);
  #endif
  return true;
}

void PlanarPartition::removeVertices() {
  // Remove unnecessary vertices completely surrounded by the same polygon
  // TODO: This can be optimised
  std::cout << "\tBefore: " << triangulation.number_of_faces() << " triangles in the triangulation" << std::endl;
  
  unsigned long long int surroundedVerticesRemoved = 0;
  Triangulation::Finite_vertices_iterator currentVertex = triangulation.finite_vertices_begin();
  while (currentVertex != triangulation.finite_vertices_end()) {
    if (triangulation.are_there_incident_constraints(currentVertex)) {
      ++currentVertex;
      continue;
    }
    
    Triangulation::Face_circulator firstFace = triangulation.incident_faces(currentVertex), currentFace = firstFace;
    ++currentFace;
    bool allEqual = true;
    while (currentFace != firstFace) {
      if (currentFace->info().getOneTag() != firstFace->info().getOneTag()) {
        allEqual = false;
        break;
      } ++currentFace;
    }
    
    if (allEqual) {
      Triangulation::Finite_vertices_iterator vertexToRemove = currentVertex;
      ++currentVertex;
      
      Point location = vertexToRemove->point();
      //Triangulation::Face_handle approximateLocation;
      PolygonHandle *tag = triangulation.incident_faces(vertexToRemove)->info().getOneTag();
      triangulation.remove(vertexToRemove);
      std::stack<Triangulation::Face_handle> stack;
      Triangulation::Face_handle emptyFace = triangulation.locate(location);
      stack.push(emptyFace);
      tagStack(stack, tag);
      
      ++surroundedVerticesRemoved;
    } else {
      ++currentVertex;
    }
  }
  std::cout << "\tRemoved " << surroundedVerticesRemoved << " surrounded vertices" << std::endl;
  std::cout << "\tAfter: " << triangulation.number_of_faces() << " triangles in the triangulation" << std::endl;
}


void PlanarPartition::removeConstraints() {
	// Remove constrained edges that have the same polygon on both sides
  unsigned long long int constrainedEdgesRemoved = 0;
  for (Triangulation::All_edges_iterator currentEdge = triangulation.all_edges_begin(); currentEdge != triangulation.all_edges_end(); ++currentEdge) {
    if (!triangulation.is_constrained(*currentEdge)) continue;
    if (currentEdge->first->info().getOneTag() == currentEdge->first->neighbor(currentEdge->second)->info().getOneTag()) {
      triangulation.remove_constrained_edge(currentEdge->first, currentEdge->second);
      ++constrainedEdgesRemoved;
    }
  }
  std::cout << "\tRemoved " << constrainedEdgesRemoved << " constrained edges" << std::endl;
}

std::list<Triangulation::Vertex_handle>*
PlanarPartition::getBoundary(Triangulation::Face_handle face, int edge, PolygonHandle *polygon) {
  std::list<Triangulation::Vertex_handle> *vertices = new std::list<Triangulation::Vertex_handle>();
	// Check clockwise edge
	if (!face->is_constrained(face->cw(edge)) && !face->neighbor(face->cw(edge))->info().hasNoTags()) {
		face->neighbor(face->cw(edge))->info().removeAllTags();
		std::list<Triangulation::Vertex_handle> *v1 = getBoundary(face->neighbor(face->cw(edge)), face->neighbor(face->cw(edge))->index(face), polygon);
    vertices->splice(vertices->end(), *v1);
		delete v1;
	}
	// Add central vertex
	vertices->push_back(face->vertex(edge));
	// Check counterclockwise edge
	if (!face->is_constrained(face->ccw(edge)) && !face->neighbor(face->ccw(edge))->info().hasNoTags()) {
		face->neighbor(face->ccw(edge))->info().removeAllTags();
		std::list<Triangulation::Vertex_handle> *v2 = getBoundary(face->neighbor(face->ccw(edge)), face->neighbor(face->ccw(edge))->index(face), polygon);
    vertices->splice(vertices->end(), *v2);
		delete v2;
	}
	return vertices;
}
