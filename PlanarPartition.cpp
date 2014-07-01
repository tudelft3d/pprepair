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

#include "PlanarPartition.h"


PlanarPartition::PlanarPartition() {
	// Registers drivers for all supported formats in OGR
	OGRRegisterAll();
	// Set internal states
	state = CREATED;
	// std::cout precision (for debugging)
	std::cout.setf(std::ios::fixed,std::ios::floatfield);
	std::cout.precision(6);
}


PlanarPartition::~PlanarPartition() {
	triangulation.clear();
}

int PlanarPartition::noPolygons() {
  return polygons.size();
}

bool PlanarPartition::addOGRdataset(std::string &file) {
  // Check if we have already made changes to the triangulation
  if (state > TRIANGULATED) {
    std::cerr << "Error: The triangulation has already been tagged. It cannot be modified!" << std::endl;
		return false;
	}
  std::cout << "Adding new dataset" << std::endl << "\t" << file << std::endl;
  std::vector<OGRFeature*> lsInputFeatures;
  getOGRFeatures(file, lsInputFeatures);
  addFeatures(lsInputFeatures);
  
  // std::cout << "# of features: " << lsInputFeatures.size() << std::endl;
  // std::cout << "deff: " << lsInputFeatures[0]->GetDefnRef()->GetName() << std::endl;
  // std::cout << "deff: " << lsInputFeatures[0]->GetFID() << std::endl;
  // std::cout << "deff: " << lsInputFeatures[1]->GetFID() << std::endl;
  // std::cout << "same ogrdefn? " << lsInputFeatures[0]->GetDefnRef()->IsSame(lsInputFeatures[1]->GetDefnRef()) << std::endl;

  return true;
}


bool PlanarPartition::addFeatures(std::vector<OGRFeature*> &lsOGRFeatures) {
  for (std::vector<OGRFeature*>::iterator f = lsOGRFeatures.begin() ; f != lsOGRFeatures.end(); f++) {
    std::vector<Polygon> polygonsVector;
    std::vector<std::list<Point> > outerRingsList;
    std::vector<std::list<Point> > innerRingsList;
    switch((*f)->GetGeometryRef()->getGeometryType()) {
      case wkbPolygon:
      case wkbPolygon25D: {
        OGRPolygon *geometry = static_cast<OGRPolygon *>((*f)->GetGeometryRef());
        outerRingsList.push_back(std::list<Point>());
        // Get outer ring
        for (int currentPoint = 0; currentPoint < (geometry->getExteriorRing()->getNumPoints() - 1); currentPoint++)
          outerRingsList.back().push_back(Point(geometry->getExteriorRing()->getX(currentPoint),
                                                geometry->getExteriorRing()->getY(currentPoint)));
        // Get inner rings
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
        polygonsVector.push_back(Polygon(oring, irings.begin(), irings.end()));
        break;
      }
      case wkbMultiPolygon: {
        OGRMultiPolygon *geometry = static_cast<OGRMultiPolygon *>((*f)->GetGeometryRef());
        // Check each polygon
        for (int currentPolygon = 0; currentPolygon < geometry->getNumGeometries(); currentPolygon++) {
          OGRPolygon *thisGeometry = static_cast<OGRPolygon *>(geometry->getGeometryRef(currentPolygon));
          outerRingsList.push_back(std::list<Point>());
          
          // Get outer ring
          for (int currentPoint = 0; currentPoint < (thisGeometry->getExteriorRing()->getNumPoints() - 1); currentPoint++)
            outerRingsList.back().push_back(Point(thisGeometry->getExteriorRing()->getX(currentPoint),
                                                  thisGeometry->getExteriorRing()->getY(currentPoint)));
          
          // Get inner rings
          innerRingsList.reserve(innerRingsList.size()+thisGeometry->getNumInteriorRings());
          for (int currentRing = 0; currentRing < thisGeometry->getNumInteriorRings(); currentRing++) {
            innerRingsList.push_back(std::list<Point>());
            for (int currentPoint = 0; currentPoint < (thisGeometry->getInteriorRing(currentRing)->getNumPoints() - 1); currentPoint++) {
              innerRingsList.back().push_back(Point(thisGeometry->getInteriorRing(currentRing)->getX(currentPoint),
                                                    thisGeometry->getInteriorRing(currentRing)->getY(currentPoint)));
            }
          }
          // TODO: fix here for multipolygon
          Ring oring(outerRingsList[0].begin(), outerRingsList[0].begin());
          outerRingsList.clear();
          std::vector<Ring> irings;
          for (unsigned int currentRing = 0; currentRing < innerRingsList.size(); currentRing++) {
            irings.push_back(Ring(innerRingsList[currentRing].begin(), innerRingsList[currentRing].end()));
            innerRingsList[currentRing].clear();
          }
          polygonsVector.push_back(Polygon(oring, irings.begin(), irings.end()));
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
      PolygonHandle *handle = new PolygonHandle(*f);
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



bool PlanarPartition::validateSingleGeom(std::vector<OGRFeature*> &lsOGRFeatures) {
  for (std::vector<OGRFeature*>::iterator it = lsOGRFeatures.begin() ; it != lsOGRFeatures.end(); ++it) {
    switch((*it)->GetGeometryRef()->getGeometryType()) {
      case wkbPolygon:
      case wkbPolygon25D: {
        OGRPolygon *geometry = (OGRPolygon *)(*it)->GetGeometryRef();
//        OGRPolygon* a = (OGRPolygon *)geometry->clone();
        break;
      }
      case wkbMultiPolygon:
      case wkbMultiPolygon25D: {
        OGRMultiPolygon *geometry = static_cast<OGRMultiPolygon *>((*it)->GetGeometryRef());
        for (int cur = 0; cur < geometry->getNumGeometries(); cur++) {
          OGRPolygon *thisGeometry = (OGRPolygon *)geometry->getGeometryRef(cur);
//            OGRPolygon* a = (OGRPolygon *)thisGeometry->clone();
//            lsOGRFeatures.push_back(a);
        }
        break;
      }
      default: {
        std::cout << "UNKNOWN GEOMETRY TYPE, skipping feature." << std::endl;
      }
    }
  }
  return true;
}

bool PlanarPartition::getOGRFeatures(std::string file, std::vector<OGRFeature*> &lsOGRFeatures) {
	OGRDataSource *dataSource = OGRSFDriverRegistrar::Open(file.c_str(), false);
	if (dataSource == NULL) {
		std::cerr << "Error: Could not open file." << std::endl;
		return false;
	}
	int numberOfLayers = dataSource->GetLayerCount();
  for (int currentLayer = 0; currentLayer < numberOfLayers; currentLayer++) {
    OGRLayer *dataLayer = dataSource->GetLayer(currentLayer);
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
  OGRDataSource::DestroyDataSource(dataSource);
  std::cout << "\tdone." << std::endl;
  return true;
}



bool PlanarPartition::addToTriangulation(const char *file, unsigned int schemaIndex) {
  // Check if we have already made changes to the triangulation
  if (state > TRIANGULATED) {
    std::cerr << "Error: The triangulation has already been tagged. It cannot be modified!" << std::endl;
		return false;
	}
  
  std::cout << "Adding a new set of polygons to the triangulation..." << std::endl;
	time_t thisTime = time(NULL);
  
  bool returnValue = io.addToTriangulation(triangulation, edgesToTag, file, schemaIndex);
  if (triangulation.number_of_faces() > 0) state = TRIANGULATED;
  
  // Triangulation stats
	std::cout << "Polygons added (" << time(NULL)-thisTime << " s). The triangulation has now:" << std::endl;
	std::cout << "\tVertices: " << triangulation.number_of_vertices() << std::endl;
	std::cout << "\tEdges: " << triangulation.tds().number_of_edges() << std::endl;
	std::cout << "\tTriangles: " << triangulation.number_of_faces() << std::endl;
  
  return returnValue;
}


bool PlanarPartition::buildPP() {
  if (state < TRIANGULATED) {
		std::cout << "No triangulation to tag!" << std::endl;
		return false;
	} if (state > TRIANGULATED) {
		std::cout << "Triangulation already tagged!" << std::endl;
		return false;
	}
  
  std::cout << "Building the PP (tagging the triangles)..." << std::endl;
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
			if ((*previousVertex)->point() == edgesToTag[currentPolygon].first[currentEdge]->point()) sameOrder = true;
			else sameOrder = false;
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
		edgesToTag[currentPolygon].second.clear();
		// Expand the tags
		tagStack(stack, polygons[currentPolygon]);
	}
	
	// Free remaining memory
	edgesToTag.clear();
	// Tag the universe
	currentFace = triangulation.infinite_face();
	stack.push(currentFace);
	tagStack(stack, &universe);
	
  state = TAGGED;
  return true;
}


void PlanarPartition::tagStack(std::stack<Triangulation::Face_handle> &stack, PolygonHandle *handle) {
	while (!stack.empty()) {
		Triangulation::Face_handle currentFace = stack.top();
		stack.pop();
		currentFace->info().addTag(handle);
		if (!currentFace->neighbor(0)->info().hasTag(handle) && !currentFace->is_constrained(0)) {
			currentFace->neighbor(0)->info().addTag(handle);
			stack.push(currentFace->neighbor(0));
		} if (!currentFace->neighbor(1)->info().hasTag(handle) && !currentFace->is_constrained(1)) {
			currentFace->neighbor(1)->info().addTag(handle);
			stack.push(currentFace->neighbor(1));
		} if (!currentFace->neighbor(2)->info().hasTag(handle) && !currentFace->is_constrained(2)) {
			currentFace->neighbor(2)->info().addTag(handle);
			stack.push(currentFace->neighbor(2));
		}
	}
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
					if (alsoUniverse || tagToAssign != &universe)
            break;
				}
        else {
					std::list<PolygonHandle *>::const_iterator randomTag = static_cast<MultiPolygonHandle *>((*randomFace)->neighbor(neighbourIndex)->info().getTags())->getHandles()->begin();
					std::advance(randomTag, rand()%numberOfTags);
					tagToAssign = *randomTag;
					if (alsoUniverse || tagToAssign != &universe)
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

bool PlanarPartition::repair(const std::string &method, bool alsoUniverse, const std::string &priority) {
	if (state < TAGGED) {
		std::cout << "Triangulation not yet tagged. Cannot repair!" << std::endl;
		return false;
	}
  if (state > TAGGED) {
		std::cout << "Triangulation already repaired!" << std::endl;
		return false;
	}
	
  time_t thisTime = time(NULL);
  bool repaired;
  if (method == "RN") {
    std::cout << "Repairing regions by random neighbour...";
    repaired = repairRN(alsoUniverse);
  }
  else if (method == "LB") {
    std::cout << "Repairing regions by longest boundary...";
    repaired = repairLB(alsoUniverse);
  }
  
  std::cout << " done" << std::endl;
	if (repaired) {
		std::cout << "Repair successful (" << time(NULL)-thisTime << " s). All polygons are now valid." << std::endl;
	}
  else {
		std::cout << "Repair of all polygons not possible (" << time(NULL)-thisTime << " s)." << std::endl;
	}
	if (repaired) state = REPAIRED;
	return repaired;
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
				if (currentLength->first != NULL && (alsoUniverse || currentLength->first != &universe)) {
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


bool PlanarPartition::tagTriangulation() {
	
	if (state < TRIANGULATED) {
		std::cout << "No triangulation to tag!" << std::endl;
		return false;
	} if (state > TRIANGULATED) {
		std::cout << "Triangulation already tagged!" << std::endl;
		return false;
	}
	
	std::cout << "Tagging..." << std::endl;
	time_t thisTime = time(NULL);
  
  bool returnValue = io.tagTriangulation(triangulation, edgesToTag);
	
	// Mark as tagged (for export)
	if (returnValue) state = TAGGED;
	
	std::cout << "Tagging done (" << time(NULL)-thisTime << " s)." << std::endl;
	
	return true;
}

bool PlanarPartition::makeAllHolesValid() {
  return io.makeAllHolesValid(triangulation);
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


bool PlanarPartition::checkValidity() {
	if (state < TAGGED) {
		std::cout << "Triangulation not yet tagged. Cannot check!" << std::endl;
		return false;
	} 
  if (state >= REPAIRED) 
    return true;

	for (Triangulation::Finite_faces_iterator currentFace = triangulation.finite_faces_begin(); currentFace != triangulation.finite_faces_end(); ++currentFace) {
		if (!(*currentFace).info().hasOneTag()) {
      return true;	// true means successful operation, not that everything is valid
    }
	}
	state = REPAIRED;
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
	
	if (ratio <= 1.0) return false;
	
	bool returnValue = io.splitRegions(triangulation, ratio);
	
	std::cout << "\tRegions split (" << time(NULL)-thisTime << " s)." << std::endl;
	
	return returnValue;
}

bool PlanarPartition::repairTrianglesByNumberOfNeighbours(bool alsoUniverse) {
	
	if (state < TAGGED) {
		std::cout << "Triangulation not yet tagged. Cannot repair!" << std::endl;
		return false;
	} if (state > TAGGED) {
		std::cout << "Triangulation already repaired!" << std::endl;
		return false;
	}
	
	std::cout << "Repairing triangles by number of neighbours..." << std::endl;
	time_t thisTime = time(NULL);
	
	bool repaired = io.repairTrianglesByNumberOfNeighbours(triangulation, alsoUniverse);
	
	if (repaired) {
		std::cout << "Repair successful (" << time(NULL)-thisTime << " s). All polygons are now valid." << std::endl;
	} else {
		std::cout << "Repair of all polygons not possible (" << time(NULL)-thisTime << " s)." << std::endl;
	}
	
	// Return whether the planar partition is now valid
	if (repaired) state = REPAIRED;
	return repaired;
}

bool PlanarPartition::repairTrianglesByAbsoluteMajority(bool alsoUniverse) {
	
	if (state < TAGGED) {
		std::cout << "Triangulation not yet tagged. Cannot repair!" << std::endl;
		return false;
	} if (state > TAGGED) {
		std::cout << "Triangulation already repaired!" << std::endl;
		return false;
	}
	
	std::cout << "Repairing triangles by absolute majority..." << std::endl;
	time_t thisTime = time(NULL);
	
	bool repaired = io.repairTrianglesByAbsoluteMajority(triangulation, alsoUniverse);
	
	if (repaired) {
		std::cout << "Repair successful (" << time(NULL)-thisTime << " s). All polygons are now valid." << std::endl;
	} else {
		std::cout << "Repair of all polygons not possible (" << time(NULL)-thisTime << " s)." << std::endl;
	}
	
	// Return whether the planar partition is now valid
	if (repaired) state = REPAIRED;
	return repaired;
}

bool PlanarPartition::repairTrianglesByLongestBoundary(bool alsoUniverse) {
	
	if (state < TAGGED) {
		std::cout << "Triangulation not yet tagged. Cannot repair!" << std::endl;
		return false;
	} if (state > TAGGED) {
		std::cout << "Triangulation already repaired!" << std::endl;
		return false;
	}
	
	std::cout << "Repairing triangles by longest boundary..." << std::endl;
	time_t thisTime = time(NULL);
	
	bool repaired = io.repairTrianglesByLongestBoundary(triangulation, alsoUniverse);
	
	if (repaired) {
		std::cout << "Repair successful (" << time(NULL)-thisTime << " s). All polygons are now valid." << std::endl;
	} else {
		std::cout << "Repair of all polygons not possible (" << time(NULL)-thisTime << " s)." << std::endl;
	}
	
	// Return whether the planar partition is now valid
	if (repaired) state = REPAIRED;
	return repaired;
}

bool PlanarPartition::repairRegionsByLongestBoundary(bool alsoUniverse) {
	
	if (state < TAGGED) {
		std::cout << "Triangulation not yet tagged. Cannot repair!" << std::endl;
		return false;
	} if (state > TAGGED) {
		std::cout << "Triangulation already repaired!" << std::endl;
		return false;
	}
	
	std::cout << "Repairing regions by longest boundary..." << std::endl;
	time_t thisTime = time(NULL);
	
	bool repaired = io.repairRegionsByLongestBoundary(triangulation, alsoUniverse);
	
	if (repaired) {
		std::cout << "Repair successful (" << time(NULL)-thisTime << " s). All polygons are now valid." << std::endl;
	} else {
		std::cout << "Repair of all polygons not possible (" << time(NULL)-thisTime << " s)." << std::endl;
	}
	
	// Return whether the planar partition is now valid
	if (repaired) state = REPAIRED;
	return repaired;
}


bool PlanarPartition::repairRegionsByRandomNeighbour(bool alsoUniverse) {
	if (state < TAGGED) {
		std::cout << "Triangulation not yet tagged. Cannot repair!" << std::endl;
		return false;
	}
  if (state > TAGGED) {
		std::cout << "Triangulation already repaired!" << std::endl;
		return false;
	}
	
	std::cout << "Repairing regions by random neighbour..." << std::endl;
	time_t thisTime = time(NULL);
	
	bool repaired = io.repairRegionsByRandomNeighbour(triangulation, alsoUniverse);
	
	if (repaired) {
		std::cout << "Repair successful (" << time(NULL)-thisTime << " s). All polygons are now valid." << std::endl;
	} else {
		std::cout << "Repair of all polygons not possible (" << time(NULL)-thisTime << " s)." << std::endl;
	}
	
	// Return whether the planar partition is now valid
	if (repaired) state = REPAIRED;
	return repaired;
}

bool PlanarPartition::repairByPriorityList(const char *file) {
	
	if (state < TAGGED) {
		std::cout << "Triangulation not yet tagged. Cannot repair!" << std::endl;
		return false;
	} if (state > TAGGED) {
		std::cout << "Triangulation already repaired!" << std::endl;
		return false;
	}
	
	std::cout << "Repairing by priority list..." << std::endl;
	time_t thisTime = time(NULL);
  
  bool repaired = io.repairByPriorityList(triangulation, file);
	
	if (repaired) {
		std::cout << "Repair successful (" << time(NULL)-thisTime << " s). All polygons are now valid." << std::endl;
	} else {
		std::cout << "Repair of all polygons not possible (" << time(NULL)-thisTime << " s)." << std::endl;
	}
	
	// Return whether the planar partition is now valid
	if (repaired) state = REPAIRED;
	return repaired;
}

bool PlanarPartition::repairEdgeMatching(const char *file) {
	
	if (state < TAGGED) {
		std::cout << "Triangulation not yet tagged. Cannot repair!" << std::endl;
		return false;
	} if (state > TAGGED) {
		std::cout << "Triangulation already repaired!" << std::endl;
		return false;
	}
	
	std::cout << "Repairing by priority list..." << std::endl;
	time_t thisTime = time(NULL);
  
  bool repaired = io.repairEdgeMatching(triangulation, file);
	
	if (repaired) {
		std::cout << "Repair successful (" << time(NULL)-thisTime << " s). All polygons are now valid." << std::endl;
	} else {
		std::cout << "Repair of all polygons not possible (" << time(NULL)-thisTime << " s)." << std::endl;
	}
	
	// Return whether the planar partition is now valid
	if (repaired) state = REPAIRED;
	return repaired;
}

bool PlanarPartition::matchSchemata() {
	
	if (state < TAGGED) {
		std::cout << "Triangulation not tagged. Cannot match schemata!" << std::endl;
		return false;
	} if (state < REPAIRED) std::cout << "Warning: Triangulation not yet repaired. There could be errors..." << std::endl;
	if (state > REPAIRED) {
		std::cout << "Polygons already reconstructed. Cannot match schemata!" << std::endl;
		return false;
	}
	
	bool returnValue = io.matchSchemata(triangulation);
	
	std::cout << "Schemata matched." << std::endl;
	
	return returnValue;
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
    if (currentTag == &universe) {
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


bool PlanarPartition::exportPolygonsSHP() {
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
  std::map<OGRFeatureDefn*,OGRDataSource*> allshps;
  std::map<OGRFeatureDefn*,OGRDataSource*>::iterator allshpsit;
  for (std::set<OGRFeatureDefn*>::iterator it = allFDefs.begin(); it != allFDefs.end(); ++it) {
    const char *driverName = "ESRI Shapefile";
    OGRSFDriver *driver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(driverName);
    if (driver == NULL) {
      std::cout << "\tError: OGR Shapefile driver not found." << std::endl;
      return false;
    }
    std::string tmp = (*it)->GetName();
    std::string outname = "/Users/hugo/temp/" + tmp + "-out.shp";
    OGRDataSource *dataSource = driver->Open(outname.c_str(), false);
    allshps[*it] = dataSource;
    if (dataSource != NULL) {
      std::cout << "\tOverwriting file..." << std::endl;
      if (driver->DeleteDataSource(dataSource->GetName())!= OGRERR_NONE) {
        std::cout << "\tError: Couldn't erase file with same name." << std::endl;
        return false;
      } OGRDataSource::DestroyDataSource(dataSource);
    }
    std::cout << "\tWriting file... " << std::endl;
    dataSource = driver->CreateDataSource(outname.c_str(), NULL);
    if (dataSource == NULL) {
      std::cout << "\tError: Could not create file." << std::endl;
      return false;
    }
    //	OGRLayer *layer = dataSource->CreateLayer("polygons", spatialReference, wkbPolygon, NULL);
    OGRLayer *layer = dataSource->CreateLayer("polygons", NULL, wkbPolygon, NULL);
    if (layer == NULL) {
      std::cout << "\tError: Could not create layer." << std::endl;
      return false;
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
    
//    OGRFeature *feature = OGRFeature::CreateFeature(layer->GetLayerDefn());
		
		f->SetGeometry(&polygon);
		// Create OGR feature
		if (layer->CreateFeature(f) != OGRERR_NONE)
      std::cout << "\tError: Could not create feature." << std::endl;
		OGRFeature::DestroyFeature(f); //-- TODO: free the features? hmmmm...
	}
  
  //-- clear memory for all the created SHP
  for (allshpsit = allshps.begin(); allshpsit != allshps.end(); ++allshpsit) {
    OGRDataSource::DestroyDataSource(allshpsit->second);
  }
  
	std::cout << "Polygons exported (" << time(NULL)-thisTime << " s)." << std::endl;
  return true;
}


bool PlanarPartition::exportTriangulation(const char *file, bool withNumberOfTags, bool withFields, bool withProvenance) {
	
	// To accept external triangulations for debugging, comment...
	if (state < TRIANGULATED || state > REPAIRED) {
		std::cout << "No triangulation to export!" << std::endl;
		return false;
	}
	
	io.exportTriangulation(triangulation, file, withNumberOfTags, withFields, withProvenance);
	
	return true;
}


void PlanarPartition::printInfo(std::ostream &ostr) {
// TODO: change so that it prints for regions not triangles?
  std::cout << "*** Triangulation ***" << std::endl;
  std::cout << "\tVertices: " << triangulation.number_of_vertices() << std::endl;
	std::cout << "\tEdges: " << triangulation.tds().number_of_edges() << std::endl;
	std::cout << "\tTriangles: " << triangulation.number_of_faces() << std::endl;

	// Number of tags
	unsigned int untagged = 0, onetag = 0, multipletags = 0, total;
	for (Triangulation::Finite_faces_iterator currentFace = triangulation.finite_faces_begin(); currentFace != triangulation.finite_faces_end(); ++currentFace) {
		if ((*currentFace).info().hasNoTags()) untagged++;
		else if ((*currentFace).info().hasOneTag()) onetag++;
		else multipletags++;
	}
  total = onetag + multipletags + untagged;
  ostr << "\tHoles:    " << untagged << " triangles (" << 100.0*untagged/total << " %)" << std::endl <<
  "\tOk:       " << onetag << " triangles (" << 100.0*onetag/total << " %)" << std::endl <<
  "\tOverlaps: " << multipletags << " triangles (" << 100.0*multipletags/total << " %)" << std::endl;
  std::cout << "*********************" << std::endl;
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
