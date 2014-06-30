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


bool PlanarPartition::addOGRdataset(std::string &file) {
  // Check if we have already made changes to the triangulation
  if (state > TRIANGULATED) {
    std::cerr << "Error: The triangulation has already been tagged. It cannot be modified!" << std::endl;
		return false;
	}
  std::cout << "Adding a new dataset to the PP: " << std::endl << file << std::endl;
  std::vector<OGRFeature*> lsInputFeatures;
  getOGRFeatures(file, lsInputFeatures);
//  validateSingleGeom(lsInputFeatures);
  addFeatures(lsInputFeatures);
  
  std::cout << "# of features: " << lsInputFeatures.size() << std::endl;
  
  std::cout << "deff: " << lsInputFeatures[0]->GetDefnRef()->GetName() << std::endl;
  std::cout << "deff: " << lsInputFeatures[0]->GetFID() << std::endl;
  std::cout << "deff: " << lsInputFeatures[1]->GetFID() << std::endl;
  std::cout << "same ogrdefn? " << lsInputFeatures[0]->GetDefnRef()->IsSame(lsInputFeatures[1]->GetDefnRef()) << std::endl;

  return true;
}


bool PlanarPartition::addFeatures(std::vector<OGRFeature*> &lsOGRFeatures) {
  for (std::vector<OGRFeature*>::iterator f = lsOGRFeatures.begin() ; f != lsOGRFeatures.end(); ++f) {
    std::vector<Polygon> polygonsVector;
    std::vector<std::list<Point> > outerRingsList;
    std::vector<std::list<Point> > innerRingsList;
    switch((*f)->GetGeometryRef()->getGeometryType()) {
      case wkbPolygon:
      case wkbPolygon25D: {
        OGRPolygon *geometry = static_cast<OGRPolygon *>((*f)->GetGeometryRef());
        outerRingsList.push_back(std::list<Point>());
        // Get outer ring
        for (int currentPoint = 0; currentPoint < geometry->getExteriorRing()->getNumPoints(); currentPoint++)
          outerRingsList.back().push_back(Point(geometry->getExteriorRing()->getX(currentPoint),
                                                geometry->getExteriorRing()->getY(currentPoint)));
        // Get inner rings
        innerRingsList.reserve(geometry->getNumInteriorRings());
        for (int currentRing = 0; currentRing < geometry->getNumInteriorRings(); currentRing++) {
          innerRingsList.push_back(std::list<Point>());
          for (int currentPoint = 0; currentPoint < geometry->getInteriorRing(currentRing)->getNumPoints(); currentPoint++) {
            innerRingsList.back().push_back(Point(geometry->getInteriorRing(currentRing)->getX(currentPoint),
                                                  geometry->getInteriorRing(currentRing)->getY(currentPoint)));
          }
        }
        
        Ring oring(outerRingsList[0].begin(), outerRingsList[0].begin());
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
          for (int currentPoint = 0; currentPoint < thisGeometry->getExteriorRing()->getNumPoints(); currentPoint++)
            outerRingsList.back().push_back(Point(thisGeometry->getExteriorRing()->getX(currentPoint),
                                                  thisGeometry->getExteriorRing()->getY(currentPoint)));
          
          // Get inner rings
          innerRingsList.reserve(innerRingsList.size()+thisGeometry->getNumInteriorRings());
          for (int currentRing = 0; currentRing < thisGeometry->getNumInteriorRings(); currentRing++) {
            innerRingsList.push_back(std::list<Point>());
            for (int currentPoint = 0; currentPoint < thisGeometry->getInteriorRing(currentRing)->getNumPoints(); currentPoint++) {
              innerRingsList.back().push_back(Point(thisGeometry->getInteriorRing(currentRing)->getX(currentPoint),
                                                    thisGeometry->getInteriorRing(currentRing)->getY(currentPoint)));
            }
          }
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
    
    for (std::vector<Polygon>::iterator currentPolygon = polygonsVector.begin(); currentPolygon != polygonsVector.end(); ++currentPolygon) {
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

bool PlanarPartition::getOGRFeatures(std::string &file, std::vector<OGRFeature*> &lsOGRFeatures) {
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
    std::cout << "\tReading layer #" << currentLayer+1 << " (" << numberOfPolygons << " polygons)..." << std::endl;
    
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

bool PlanarPartition::checkValidity() {
	
	if (state < TAGGED) {
		std::cout << "Triangulation not yet tagged. Cannot check!" << std::endl;
		return false;
	} if (state >= REPAIRED) return true;
	
	for (Triangulation::Finite_faces_iterator currentFace = triangulation.finite_faces_begin(); currentFace != triangulation.finite_faces_end(); ++currentFace) {
		if (!(*currentFace).info().hasOneTag()) return true;	// true means successful operation, not that everything is valid
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
	} if (state > TAGGED) {
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
  
  io.removeConstraints(triangulation);
  if (removeVertices) io.removeVertices(triangulation);
  bool returnValue = io.reconstructPolygons(triangulation, outputPolygons);
  
  // Mark as reconstructed
	if (returnValue) state = RECONSTRUCTED;
	
	std::cout << "Polygons reconstructed (" << time(NULL)-thisTime << " s)." << std::endl;
  return returnValue;
}

bool PlanarPartition::exportPolygons(const char *file, bool withProvenance) {
	
	if (state < RECONSTRUCTED) {
		std::cout << "Polygons have not been reconstructed. Nothing to export!" << std::endl;
		return false;
	}
	
	std::cout << "Exporting polygons..." << std::endl;
	time_t thisTime = time(NULL);
  
  bool returnValue = io.exportPolygons(outputPolygons, file, withProvenance);
	
	std::cout << "Polygons exported (" << time(NULL)-thisTime << " s)." << std::endl;
  return returnValue;
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

}
