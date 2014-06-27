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

#include "IOWorker.h"

IOWorker::IOWorker() {
  startingSearchFace = Triangulation::Face_handle();
}

bool IOWorker::addToTriangulation(Triangulation &triangulation, TaggingVector &edgesToTag, const char *file, unsigned int schemaIndex) {
  // Open file
	OGRDataSource *dataSource = OGRSFDriverRegistrar::Open(file, false);
	if (dataSource == NULL) {
		std::cerr << "Error: Could not open file." << std::endl;
		return false;
	}
  
  char *name = new char[strlen(dataSource->GetName())+1];
	strcpy(name, dataSource->GetName());
	fileNames.push_back(name);
	std::cout << "\tPath: " << name << std::endl;
	std::cout << "\tType: " << dataSource->GetDriver()->GetName() << std::endl;
	int numberOfLayers = dataSource->GetLayerCount();
	std::cout << "\tLayers: " << numberOfLayers << std::endl;
  
  // Read layer by layer
  for (int currentLayer = 0; currentLayer < numberOfLayers; currentLayer++) {
    OGRLayer *dataLayer = dataSource->GetLayer(currentLayer);
    dataLayer->ResetReading();
    OGRSpatialReference* tmp = dataLayer->GetSpatialRef();
    if ( (tmp != NULL) && (spatialReference != NULL) )
      spatialReference = tmp->CloneGeogCS();
		
		unsigned int numberOfPolygons = dataLayer->GetFeatureCount(true);
		std::cout << "\tReading layer #" << currentLayer+1 << " (" << numberOfPolygons << " polygons)...";
		polygons.reserve(polygons.size()+numberOfPolygons);
    
    // Check fields and the schema type
    OGRFeatureDefn *layerDefinition = dataLayer->GetLayerDefn();
    insertToStream(std::cout, layerDefinition, 1, schemaIndex);
    
    // If it's the first input file, assign the schema type of it
    if (triangulation.number_of_faces() == 0) {
      schemaFieldType = layerDefinition->GetFieldDefn(schemaIndex)->GetType();
    } // Otherwise, check if it matches the previous one
    else {
      if (layerDefinition->GetFieldDefn(schemaIndex)->GetType() != schemaFieldType) {
        std::cerr << "\tError: The schema field type in this layer is incompatible with the previous one. Skipped." << std::endl;
        continue;
      }
    }
    
    // Save the field names and types
		for (int currentField = 0; currentField < layerDefinition->GetFieldCount(); currentField++) {
			OGRFieldDefn *fieldDefinition = layerDefinition->GetFieldDefn(currentField);
			FieldDefinition *newField = new FieldDefinition(fieldDefinition->GetNameRef(), fieldDefinition->GetType(), fieldDefinition->GetJustify(), fieldDefinition->GetWidth(), fieldDefinition->GetPrecision());
			unsigned int currentCheck;
			for (currentCheck = 0; currentCheck < fields.size(); currentCheck++) {
				if (newField->matches(fields[currentCheck])) break;
			} if (currentCheck == (unsigned int)fields.size()) {
				// It's a new field
				fields.push_back(newField);
				fieldEquivalencies[FieldDescriptor(name, currentLayer, currentField)] = ((unsigned int)fields.size())-1;
			} else {
				// The field matches an older one, don't add
				delete newField;
				fieldEquivalencies[FieldDescriptor(name, currentLayer, currentField)] = currentCheck;
			}
		}
    
    // Reads all features in this layer
		OGRFeature *feature;
		while ((feature = dataLayer->GetNextFeature()) != NULL) {
			
			// STEP 1: Get polygons from input
			std::vector<std::list<Point> > outerRingsList;
			std::vector<std::list<Point> > innerRingsList;
			switch(feature->GetGeometryRef()->getGeometryType()) {
          
          // Most typical case, receiving polygons
        case wkbPolygon:
        case wkbPolygon25D: {
					OGRPolygon *geometry = static_cast<OGRPolygon *>(feature->GetGeometryRef());
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
					} break;
				}
					
          // Receiving multi polygons
				case wkbMultiPolygon: {
					OGRMultiPolygon *geometry = static_cast<OGRMultiPolygon *>(feature->GetGeometryRef());
					
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
					} break;
				}
          
				default:
					std::cerr << "\tFeature #" << feature->GetFID() << ": unsupported type (";
					insertToStream(std::cout, feature->GetGeometryRef()->getGeometryType());
					std::cerr << "). Skipped." << std::endl;
					continue;
					break;
          
          // TODO: Implement other cases: points, lines, containers with multiple features, etc.
			}
			
			// STEP 2: Check validity of individual polygons
			//  it's more efficient doing this during creation, but this is more readable and maintainable (check SVN v59).
			//  After all, this is not the main focus here.
			
			std::vector<Polygon> polygonsVector;
			// CHECKS ON VERTICES
      
      // Remove repeated vertices. One per ring is considered normal, since some specifications allow or require it (first == last).
      for (unsigned int currentRing = 0; currentRing < outerRingsList.size(); ++currentRing) {
        if (removeDuplicateVertices(outerRingsList[currentRing]) > 1)
          std::cout << "\tFeature #" << feature->GetFID() << ": duplicate vertices in outer boundary #" << currentRing << ". Removed duplicates." << std::endl;
      } for (unsigned int currentRing = 0; currentRing < innerRingsList.size(); ++currentRing) {
        if (removeDuplicateVertices(innerRingsList[currentRing]) > 1)
          std::cout << "\tFeature #" << feature->GetFID() << ": duplicate vertices in inner boundary #" << currentRing << ". Removed duplicates." << std::endl;
      }
      
      // Trivial check for rings with less than 3 vertices. The ones with 3 or more vertices will be done in the triangulation
      for (int currentRing = 0; currentRing < (int)outerRingsList.size(); ++currentRing) {
        if (outerRingsList[currentRing].size() < 3) {
          std::cout << "\tFeature #" << feature->GetFID() << ": less than 3 vertices in outer boundary #" << currentRing << ". Removed." << std::endl;
          outerRingsList.erase(outerRingsList.begin()+currentRing);
          --currentRing;
        }
      } for (int currentRing = 0; currentRing < (int)innerRingsList.size(); ++currentRing) {
        if (innerRingsList[currentRing].size() < 3) {
          std::cout << "\tFeature #" << feature->GetFID() << ": less than 3 vertices in inner boundary #" << currentRing << ". Removed." << std::endl;
          innerRingsList.erase(innerRingsList.begin()+currentRing);
          --currentRing;
        }
      }
      
      // CHECKS ON RINGS
      
      // Let's move on to the CGAL data structures for Rings
      std::vector<Ring> outerRings;
      std::vector<Ring> innerRings;
      outerRings.reserve(outerRingsList.size());
      innerRings.reserve(innerRingsList.size());
      for (unsigned int currentRing = 0; currentRing < outerRingsList.size(); currentRing++) {
        outerRings.push_back(Ring(outerRingsList[currentRing].begin(), outerRingsList[currentRing].end()));
        outerRingsList[currentRing].clear();
      } for (unsigned int currentRing = 0; currentRing < innerRingsList.size(); currentRing++) {
        innerRings.push_back(Ring(innerRingsList[currentRing].begin(), innerRingsList[currentRing].end()));
        innerRingsList[currentRing].clear();
      }
      
      // Split self touching rings and correct winding
      std::vector<Ring *> outerRingsToBuild;
      std::vector<Ring *> innerRingsToClassify;
      std::vector<std::vector<Ring> > innerRingsToBuild;
      
      // Get outer rings
      for (unsigned int currentRings = 0; currentRings < outerRings.size(); currentRings++) {
        if (!outerRings[currentRings].is_simple()) {
          std::cout << "\tFeature #" << feature->GetFID() << " (" << outerRings[currentRings].size() << " vertices): self intersecting outer boundary #" << currentRings << ". Split." << std::endl;
          std::vector<Ring *> receivedRings = splitRing(outerRings[currentRings]);
          for (std::vector<Ring *>::iterator currentRing = receivedRings.begin(); currentRing != receivedRings.end(); ++currentRing) {
            if ((*currentRing)->is_clockwise_oriented()) {
              outerRingsToBuild.push_back(*currentRing);
            } else {
              innerRingsToClassify.push_back(*currentRing);
            }
          }
        } else {
          if (outerRings[currentRings].is_counterclockwise_oriented()) {
            std::cout << "\tFeature #" << feature->GetFID() << ": incorrect winding in outer boundary #" << currentRings << ". Reversed." << std::endl;
            outerRings[currentRings].reverse_orientation();
          } outerRingsToBuild.push_back(new Ring(outerRings[currentRings]));
          outerRings[currentRings].clear();
        }
      }
      
      // Get inner rings
      for (unsigned int currentRings = 0; currentRings < innerRings.size(); currentRings++) {
        if (!innerRings[currentRings].is_simple()) {
          std::cout << "\tFeature #" << feature->GetFID() << " (" << innerRings[currentRings].size() << " vertices): self intersecting inner boundary #" << currentRings << ". Split." << std::endl;
          std::vector<Ring *> receivedRings = splitRing(innerRings[currentRings]);
          for (std::vector<Ring *>::iterator currentRing = receivedRings.begin(); currentRing != receivedRings.end(); ++currentRing) {
            if ((*currentRing)->is_clockwise_oriented()) {
              innerRingsToClassify.push_back(*currentRing);
            } else {
              outerRingsToBuild.push_back(*currentRing);
            }
          }
        } else {
          if (innerRings[currentRings].is_clockwise_oriented()) {
            std::cout << "\tFeature #" << feature->GetFID() << ": incorrect winding in inner boundary #" << currentRings << ". Reversed." << std::endl;
            innerRings[currentRings].reverse_orientation();
          } innerRingsToClassify.push_back(new Ring(innerRings[currentRings]));
          innerRings[currentRings].clear();
        }
      }
      
      // Make space for inner rings
      for (std::vector<Ring *>::iterator currentRing = outerRingsToBuild.begin(); currentRing != outerRingsToBuild.end(); ++currentRing) {
        innerRingsToBuild.push_back(std::vector<Ring>());
      }
      
      // Put inner rings into the correct outer ring (and likely other ones). Incorrectly nested rings are found here.
      if (outerRingsToBuild.size() == 0) {
        // Outer ring had no area or there wasn't any. Delete all inner rings
        std::cout << "\tFeature #" << feature->GetFID() << ": zero area outer boundary. Inner boundaries removed." << std::endl;
        for (std::vector<Ring *>::iterator currentRing = innerRingsToClassify.begin(); currentRing != innerRingsToClassify.end(); ++currentRing) {
          delete *currentRing;
        }
      }
      
      // Now check them and put them in place
      else if (innerRingsToClassify.size() > 0) {
        testRings(outerRingsToBuild, innerRingsToClassify, innerRingsToBuild, feature->GetFID());
      }
      
      // Let's move on to CGAL data structures for Polygons
      for (unsigned int currentPolygon = 0; currentPolygon < outerRingsToBuild.size(); ++currentPolygon) {
        polygonsVector.push_back(Polygon(*outerRingsToBuild[currentPolygon], innerRingsToBuild[currentPolygon].begin(), innerRingsToBuild[currentPolygon].end()));
      } outerRingsToBuild.clear();
      innerRingsToBuild.clear();
			
			// STEP 3: Introduce edges as constraints in the triangulation
      for (std::vector<Polygon>::iterator currentPolygon = polygonsVector.begin(); currentPolygon != polygonsVector.end(); ++currentPolygon) {
				
				// Create and save polygon handle
//				PolygonHandle *handle = new PolygonHandle(schemaIndex, fileNames.back(), currentLayer, feature->GetFID());
//				polygons.push_back(handle);
				
				// Save other attributes to put back later
//				copyFields(feature, handle);
				
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
				} for (Polygon::Hole_const_iterator currentRing = currentPolygon->holes_begin(); currentRing != currentPolygon->holes_end(); ++currentRing) {
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
			
			// Free memory
			polygonsVector.clear();
			
			// Free OGR feature
			OGRFeature::DestroyFeature(feature);
		}
  }
  
  // Free OGR data source
	OGRDataSource::DestroyDataSource(dataSource);
  
  return true;
}

bool IOWorker::tagTriangulation(Triangulation &triangulation, TaggingVector &edgesToTag) {
	
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
	
	return true;
}

bool IOWorker::makeAllHolesValid(Triangulation &triangulation) {
	
	for (Triangulation::Finite_faces_iterator currentFace = triangulation.finite_faces_begin(); currentFace != triangulation.finite_faces_end(); ++currentFace) {
		if (currentFace->info().hasNoTags()) {
			currentFace->info().addTag(&universe);
		}
	}
	
	return true;
}

bool IOWorker::splitRegions(Triangulation &triangulation, double ratio) {
	
	double shortSide, longSide, thisSide;
	unsigned int whichSide, splits = 0;
	
	for (Triangulation::Finite_faces_iterator currentFace = triangulation.finite_faces_begin(); currentFace != triangulation.finite_faces_end(); ++currentFace) {
		
		// Check for the longest and shortest sides
		shortSide = longSide = sqrt(CGAL::to_double(triangulation.segment(currentFace, 0).squared_length()));
		whichSide = 0;
		thisSide = sqrt(CGAL::to_double(triangulation.segment(currentFace, 1).squared_length()));
		if (thisSide > longSide) longSide = thisSide;
		else if (thisSide < shortSide) {
			shortSide = thisSide;
			whichSide = 1;
		} thisSide = sqrt(CGAL::to_double(triangulation.segment(currentFace, 2).squared_length()));
		if (thisSide > longSide) longSide = thisSide;
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
	
	return true;
}

bool IOWorker::repairTrianglesByNumberOfNeighbours(Triangulation &triangulation, bool alsoUniverse) {
	
	bool repaired = true;
	
	// Use a temporary vector to make it deterministic and order independent
	std::vector<std::pair<Triangulation::Face_handle, PolygonHandle *> > facesToRepair;
	for (Triangulation::Finite_faces_iterator currentFace = triangulation.finite_faces_begin(); currentFace != triangulation.finite_faces_end(); ++currentFace) {
		std::map<PolygonHandle *, unsigned int> tagCount;
		if (!currentFace->info().hasOneTag()) {
			
			// Count the number of times each tag appears
			addtoCount(tagCount, currentFace->neighbor(0)->info().getTags());
			addtoCount(tagCount, currentFace->neighbor(1)->info().getTags());
			addtoCount(tagCount, currentFace->neighbor(2)->info().getTags());
      
			// Find the tag with highest count
			unsigned int maxCount = 0;
			std::map<PolygonHandle *, unsigned int>::iterator mostTimesAppeared = tagCount.end();
			for (std::map<PolygonHandle *, unsigned int>::iterator currentCount = tagCount.begin(); currentCount != tagCount.end(); ++currentCount) {
				if (currentCount->first != NULL && (alsoUniverse || currentCount->first != &universe)) {
					if (currentCount->second > maxCount && (currentFace->info().hasTag(currentCount->first) || currentFace->info().hasNoTags())) {
						currentCount->second = maxCount;
						mostTimesAppeared = currentCount;
					} else if (currentCount->second == maxCount) {
						mostTimesAppeared = tagCount.end();
					}
				}
			}
			
			// Assign the triangle to the tag with the highest count (if there is one)
			if (mostTimesAppeared == tagCount.end()) repaired = false;
			else facesToRepair.push_back(std::pair<Triangulation::Face_handle, PolygonHandle *>(currentFace, mostTimesAppeared->first));
		}
	}
	
	// Re-tag faces in the vector
	for (std::vector<std::pair<Triangulation::Face_handle, PolygonHandle *> >::iterator currentFace = facesToRepair.begin(); currentFace != facesToRepair.end(); ++currentFace) {
		currentFace->first->info().removeAllTags();
		currentFace->first->info().addTag(currentFace->second);
	}
	
	return repaired;
}

bool IOWorker::repairTrianglesByAbsoluteMajority(Triangulation &triangulation, bool alsoUniverse) {
	
	bool repaired = true;
	
	// Put faces to repair in the vector
	// Use a temporary vector to make it deterministic and order independent
	std::vector<std::pair<Triangulation::Face_handle, Triangulation::Face_handle> > facesToRepair;
	for (Triangulation::Finite_faces_iterator currentFace = triangulation.finite_faces_begin(); currentFace != triangulation.finite_faces_end(); ++currentFace) {
		if (!currentFace->info().hasOneTag()) {
			if (currentFace->neighbor(0)->info().hasOneTag() && currentFace->neighbor(1)->info().hasOneTag() &&
          currentFace->neighbor(0)->info().getTags() == currentFace->neighbor(1)->info().getTags() &&
          (currentFace->info().hasTag(currentFace->neighbor(0)->info().getTags()) || currentFace->info().hasNoTags())) {
				if ((currentFace->neighbor(0)->info().getTags() != &universe &&
             currentFace->neighbor(0)->info().getTags() != NULL) ||
            alsoUniverse) {
					facesToRepair.push_back(std::pair<Triangulation::Face_handle, Triangulation::Face_handle>(currentFace, currentFace->neighbor(0)));
				}
			} else if (currentFace->neighbor(0)->info().hasOneTag() && currentFace->neighbor(2)->info().hasOneTag() &&
                 currentFace->neighbor(0)->info().getTags() == currentFace->neighbor(2)->info().getTags() &&
                 (currentFace->info().hasTag(currentFace->neighbor(2)->info().getTags()) || currentFace->info().hasNoTags())) {
				if ((currentFace->neighbor(2)->info().getTags() != &universe &&
             currentFace->neighbor(2)->info().getTags() != NULL) ||
            alsoUniverse) {
					facesToRepair.push_back(std::pair<Triangulation::Face_handle, Triangulation::Face_handle>(currentFace, currentFace->neighbor(2)));
				}
			} else if (currentFace->neighbor(1)->info().hasOneTag() && currentFace->neighbor(2)->info().hasOneTag() &&
                 currentFace->neighbor(1)->info().getTags() == currentFace->neighbor(2)->info().getTags() &&
                 (currentFace->info().hasTag(currentFace->neighbor(1)->info().getTags()) || currentFace->info().hasNoTags())) {
				if ((currentFace->neighbor(1)->info().getTags() != &universe &&
             currentFace->neighbor(1)->info().getTags() != NULL) ||
            alsoUniverse) {
					facesToRepair.push_back(std::pair<Triangulation::Face_handle, Triangulation::Face_handle>(currentFace, currentFace->neighbor(1)));
				}
			} else {
				repaired = false;
			}
		}
	}
	
	// Re-tag faces in the vector
	for (std::vector<std::pair<Triangulation::Face_handle, Triangulation::Face_handle> >::iterator currentFace = facesToRepair.begin();
       currentFace != facesToRepair.end();
       ++currentFace) {
		currentFace->first->info().removeAllTags();
		currentFace->first->info().addTag(currentFace->second->info().getTags());
	}
	
	return repaired;
}

bool IOWorker::repairTrianglesByLongestBoundary(Triangulation &triangulation, bool alsoUniverse) {
	
	bool repaired = true;
	
	// Use a temporary vector to make it deterministic and order independent
	std::vector<std::pair<Triangulation::Face_handle, PolygonHandle *> > facesToRepair;
	for (Triangulation::Finite_faces_iterator currentFace = triangulation.finite_faces_begin(); currentFace != triangulation.finite_faces_end(); ++currentFace) {
		std::map<PolygonHandle *, double> tagBoundaryLength;
		if (!currentFace->info().hasOneTag()) {
			
			// Add up the boundary for each tag
			addToLength(tagBoundaryLength, currentFace->neighbor(0)->info().getTags(), sqrt(CGAL::to_double(triangulation.segment(currentFace, 0).squared_length())));
			addToLength(tagBoundaryLength, currentFace->neighbor(1)->info().getTags(), sqrt(CGAL::to_double(triangulation.segment(currentFace, 1).squared_length())));
			addToLength(tagBoundaryLength, currentFace->neighbor(2)->info().getTags(), sqrt(CGAL::to_double(triangulation.segment(currentFace, 2).squared_length())));
      
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
			
			// Assign the triangle to the tag with the longest boundary (if there is one)
			if (longest == tagBoundaryLength.end()) repaired = false;
			else facesToRepair.push_back(std::pair<Triangulation::Face_handle, PolygonHandle *>(currentFace, longest->first));
		}
	}
	
	// Re-tag faces in the vector
	for (std::vector<std::pair<Triangulation::Face_handle, PolygonHandle *> >::iterator currentFace = facesToRepair.begin(); currentFace != facesToRepair.end(); ++currentFace) {
		currentFace->first->info().removeAllTags();
		currentFace->first->info().addTag(currentFace->second);
	}
	
	return repaired;
}

bool IOWorker::repairRegionsByLongestBoundary(Triangulation &triangulation, bool alsoUniverse) {
	
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

bool IOWorker::repairRegionsByRandomNeighbour(Triangulation &triangulation, bool alsoUniverse) {
	
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
				std::advance(randomFace, rand()%facesInRegion.size());
				int neighbourIndex = rand()%3;
				unsigned int numberOfTags = (*randomFace)->neighbor(neighbourIndex)->info().numberOfTags();
				if (numberOfTags == 0) continue;
				if (numberOfTags == 1) {
					tagToAssign = (*randomFace)->neighbor(neighbourIndex)->info().getTags();
					if (alsoUniverse || tagToAssign != &universe) break;
				} else {
					std::list<PolygonHandle *>::const_iterator randomTag = static_cast<MultiPolygonHandle *>((*randomFace)->neighbor(neighbourIndex)->info().getTags())->getHandles()->begin();
					std::advance(randomTag, rand()%numberOfTags);
					tagToAssign = *randomTag;
					if (alsoUniverse || tagToAssign != &universe) break;
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

bool IOWorker::repairByPriorityList(Triangulation &triangulation, const char *file) {
	
	// Process priority file
	std::ifstream priorityFile;
	priorityFile.open(file, std::ios::in);
	if (!priorityFile.is_open()) {
		std::cout << "Priority file could not be opened." << std::endl;
		return false;
	} std::map<Field *, unsigned int, FieldComparator> priorityMap;
	unsigned int currentPriority = 0;
	while (!priorityFile.eof()) {
		switch (schemaFieldType) {
			case OFTString: {
				std::string fieldAsString;
				std::getline(priorityFile, fieldAsString);		// If we deal with strings take a whole line (since spaces could be valid)
				StringField *newField = new StringField(fieldAsString.c_str());
				priorityMap[newField] = currentPriority;
				break;
			} case OFTReal: {
				double fieldAsDouble;
				priorityFile >> fieldAsDouble;
				DoubleField *newField = new DoubleField(fieldAsDouble);
				priorityMap[newField] = currentPriority;
			} case OFTInteger: {
				int fieldAsInt;
				priorityFile >> fieldAsInt;
				IntField *newField = new IntField(fieldAsInt);
				priorityMap[newField] = currentPriority;
			} default: {
				std::cout << "Field type not supported." << std::endl;
				std::string fieldAsString;
				std::getline(priorityFile, fieldAsString);
				break;
			}
		} ++currentPriority;
	} priorityFile.close();
	
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
			for (std::set<Triangulation::Face_handle>::iterator currentFaceInRegion = facesInRegion.begin(); currentFaceInRegion != facesInRegion.end(); ++currentFaceInRegion) {
				// Gap, check neighbours
				if ((*currentFaceInRegion)->info().hasNoTags()) {
					if (!(*currentFaceInRegion)->neighbor(0)->info().hasNoTags()) {
						if ((*currentFaceInRegion)->neighbor(0)->info().hasOneTag() && (*currentFaceInRegion)->neighbor(0)->info().getTags() != &universe) {
							if (priorityMap[(*currentFaceInRegion)->neighbor(0)->info().getTags()->getSchemaField()] < priorityOfTag) {
								priorityOfTag = priorityMap[(*currentFaceInRegion)->neighbor(0)->info().getTags()->getSchemaField()];
								tagToAssign = (*currentFaceInRegion)->neighbor(0)->info().getTags();
							}
						} else {
							MultiPolygonHandle *handle = static_cast<MultiPolygonHandle *>((*currentFaceInRegion)->neighbor(0)->info().getTags());
							for (std::list<PolygonHandle *>::const_iterator currentTag = handle->getHandles()->begin(); currentTag != handle->getHandles()->end(); ++currentTag) {
								if (*currentTag == &universe) continue;
								if (priorityMap[(*currentTag)->getSchemaField()] < priorityOfTag) {
									priorityOfTag = priorityMap[(*currentTag)->getSchemaField()];
									tagToAssign = *currentTag;
								}
							}
						}
					} if (!(*currentFaceInRegion)->neighbor(1)->info().hasNoTags()) {
						if ((*currentFaceInRegion)->neighbor(1)->info().hasOneTag() && (*currentFaceInRegion)->neighbor(1)->info().getTags() != &universe) {
							if (priorityMap[(*currentFaceInRegion)->neighbor(1)->info().getTags()->getSchemaField()] < priorityOfTag) {
								priorityOfTag = priorityMap[(*currentFaceInRegion)->neighbor(1)->info().getTags()->getSchemaField()];
								tagToAssign = (*currentFaceInRegion)->neighbor(1)->info().getTags();
							}
						} else {
							MultiPolygonHandle *handle = static_cast<MultiPolygonHandle *>((*currentFaceInRegion)->neighbor(1)->info().getTags());
							for (std::list<PolygonHandle *>::const_iterator currentTag = handle->getHandles()->begin(); currentTag != handle->getHandles()->end(); ++currentTag) {
								if (*currentTag == &universe) continue;
								if (priorityMap[(*currentTag)->getSchemaField()] < priorityOfTag) {
									priorityOfTag = priorityMap[(*currentTag)->getSchemaField()];
									tagToAssign = *currentTag;
								}
							}
						}
					} if (!(*currentFaceInRegion)->neighbor(2)->info().hasNoTags()) {
						if ((*currentFaceInRegion)->neighbor(2)->info().hasOneTag() && (*currentFaceInRegion)->neighbor(2)->info().getTags() != &universe) {
							if (priorityMap[(*currentFaceInRegion)->neighbor(2)->info().getTags()->getSchemaField()] < priorityOfTag) {
								priorityOfTag = priorityMap[(*currentFaceInRegion)->neighbor(2)->info().getTags()->getSchemaField()];
								tagToAssign = (*currentFaceInRegion)->neighbor(2)->info().getTags();
							}
						} else {
							MultiPolygonHandle *handle = static_cast<MultiPolygonHandle *>((*currentFaceInRegion)->neighbor(2)->info().getTags());
							for (std::list<PolygonHandle *>::const_iterator currentTag = handle->getHandles()->begin(); currentTag != handle->getHandles()->end(); ++currentTag) {
								if (*currentTag == &universe) continue;
								if (priorityMap[(*currentTag)->getSchemaField()] < priorityOfTag) {
									priorityOfTag = priorityMap[(*currentTag)->getSchemaField()];
									tagToAssign = *currentTag;
								}
							}
						}
					}
				}
				
				// Overlap, check this one
				else {
					if ((*currentFaceInRegion)->info().hasOneTag() && (*currentFaceInRegion)->info().getTags() != &universe) {
						if (priorityMap[(*currentFaceInRegion)->info().getTags()->getSchemaField()] < priorityOfTag) {
							priorityOfTag = priorityMap[(*currentFaceInRegion)->info().getTags()->getSchemaField()];
							tagToAssign = (*currentFaceInRegion)->info().getTags();
						}
					} else {
						MultiPolygonHandle *handle = static_cast<MultiPolygonHandle *>((*currentFaceInRegion)->info().getTags());
						for (std::list<PolygonHandle *>::const_iterator currentTag = handle->getHandles()->begin(); currentTag != handle->getHandles()->end(); ++currentTag) {
							if (*currentTag == &universe) continue;
							if (priorityMap[(*currentTag)->getSchemaField()] < priorityOfTag) {
								priorityOfTag = priorityMap[(*currentTag)->getSchemaField()];
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
	for (std::vector<std::pair<Triangulation::Face_handle, PolygonHandle *> >::iterator currentFace = facesToRepair.begin(); currentFace != facesToRepair.end(); ++currentFace) {
		currentFace->first->info().removeAllTags();
		currentFace->first->info().addTag(currentFace->second);
	}
	
	return true;
}

bool IOWorker::repairEdgeMatching(Triangulation &triangulation, const char *file) {
	
	// Process priority file
	std::ifstream priorityFile;
	priorityFile.open(file, std::ios::in);
	if (!priorityFile.is_open()) {
		std::cout << "Priority file could not be opened." << std::endl;
		return false;
	} std::map<Field *, unsigned int, FieldComparator> priorityMap;
	unsigned int currentPriority = 0;
	while (!priorityFile.eof()) {
		switch (schemaFieldType) {
			case OFTString: {
				std::string fieldAsString;
				std::getline(priorityFile, fieldAsString);		// If we deal with strings take a whole line (since spaces could be valid)
				StringField *newField = new StringField(fieldAsString.c_str());
				priorityMap[newField] = currentPriority;
				break;
			} case OFTReal: {
				double fieldAsDouble;
				priorityFile >> fieldAsDouble;
				DoubleField *newField = new DoubleField(fieldAsDouble);
				priorityMap[newField] = currentPriority;
			} case OFTInteger: {
				int fieldAsInt;
				priorityFile >> fieldAsInt;
				IntField *newField = new IntField(fieldAsInt);
				priorityMap[newField] = currentPriority;
			} default: {
				std::cout << "Field type not supported." << std::endl;
				std::string fieldAsString;
				std::getline(priorityFile, fieldAsString);
				break;
			}
		} ++currentPriority;
	} priorityFile.close();
	
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
			unsigned int priorityOfTagh = UINT_MAX;
      unsigned int priorityOfTagg = 0;
			for (std::set<Triangulation::Face_handle>::iterator currentFaceInRegion = facesInRegion.begin(); currentFaceInRegion != facesInRegion.end(); ++currentFaceInRegion) {
				// Gap, check neighbours
				if ((*currentFaceInRegion)->info().hasNoTags()) {
					if (!(*currentFaceInRegion)->neighbor(0)->info().hasNoTags()) {
						if ((*currentFaceInRegion)->neighbor(0)->info().hasOneTag() && (*currentFaceInRegion)->neighbor(0)->info().getTags() != &universe) {
							if (priorityMap[(*currentFaceInRegion)->neighbor(0)->info().getTags()->getSchemaField()] >= priorityOfTagg) {
								priorityOfTagg = priorityMap[(*currentFaceInRegion)->neighbor(0)->info().getTags()->getSchemaField()];
								tagToAssign = (*currentFaceInRegion)->neighbor(0)->info().getTags();
							}
						} else {
							MultiPolygonHandle *handle = static_cast<MultiPolygonHandle *>((*currentFaceInRegion)->neighbor(0)->info().getTags());
							for (std::list<PolygonHandle *>::const_iterator currentTag = handle->getHandles()->begin(); currentTag != handle->getHandles()->end(); ++currentTag) {
								if (*currentTag == &universe) continue;
								if (priorityMap[(*currentTag)->getSchemaField()] < priorityOfTagg) {
									priorityOfTagg = priorityMap[(*currentTag)->getSchemaField()];
									tagToAssign = *currentTag;
								}
							}
						}
					} if (!(*currentFaceInRegion)->neighbor(1)->info().hasNoTags()) {
						if ((*currentFaceInRegion)->neighbor(1)->info().hasOneTag() && (*currentFaceInRegion)->neighbor(1)->info().getTags() != &universe) {
							if (priorityMap[(*currentFaceInRegion)->neighbor(1)->info().getTags()->getSchemaField()] >= priorityOfTagg) {
								priorityOfTagg = priorityMap[(*currentFaceInRegion)->neighbor(1)->info().getTags()->getSchemaField()];
								tagToAssign = (*currentFaceInRegion)->neighbor(1)->info().getTags();
							}
						} else {
							MultiPolygonHandle *handle = static_cast<MultiPolygonHandle *>((*currentFaceInRegion)->neighbor(1)->info().getTags());
							for (std::list<PolygonHandle *>::const_iterator currentTag = handle->getHandles()->begin(); currentTag != handle->getHandles()->end(); ++currentTag) {
								if (*currentTag == &universe) continue;
								if (priorityMap[(*currentTag)->getSchemaField()] < priorityOfTagg) {
									priorityOfTagg = priorityMap[(*currentTag)->getSchemaField()];
									tagToAssign = *currentTag;
								}
							}
						}
					} if (!(*currentFaceInRegion)->neighbor(2)->info().hasNoTags()) {
						if ((*currentFaceInRegion)->neighbor(2)->info().hasOneTag() && (*currentFaceInRegion)->neighbor(2)->info().getTags() != &universe) {
							if (priorityMap[(*currentFaceInRegion)->neighbor(2)->info().getTags()->getSchemaField()] < priorityOfTagg) {
								priorityOfTagg = priorityMap[(*currentFaceInRegion)->neighbor(2)->info().getTags()->getSchemaField()];
								tagToAssign = (*currentFaceInRegion)->neighbor(2)->info().getTags();
							}
						} else {
							MultiPolygonHandle *handle = static_cast<MultiPolygonHandle *>((*currentFaceInRegion)->neighbor(2)->info().getTags());
							for (std::list<PolygonHandle *>::const_iterator currentTag = handle->getHandles()->begin(); currentTag != handle->getHandles()->end(); ++currentTag) {
								if (*currentTag == &universe) continue;
								if (priorityMap[(*currentTag)->getSchemaField()] < priorityOfTagg) {
									priorityOfTagg = priorityMap[(*currentTag)->getSchemaField()];
									tagToAssign = *currentTag;
								}
							}
						}
					}
				}
				
				// Overlap, check this one
				else {
					if ((*currentFaceInRegion)->info().hasOneTag() && (*currentFaceInRegion)->info().getTags() != &universe) {
						if (priorityMap[(*currentFaceInRegion)->info().getTags()->getSchemaField()] < priorityOfTagh) {
							priorityOfTagh = priorityMap[(*currentFaceInRegion)->info().getTags()->getSchemaField()];
							tagToAssign = (*currentFaceInRegion)->info().getTags();
						}
					} else {
						MultiPolygonHandle *handle = static_cast<MultiPolygonHandle *>((*currentFaceInRegion)->info().getTags());
						for (std::list<PolygonHandle *>::const_iterator currentTag = handle->getHandles()->begin(); currentTag != handle->getHandles()->end(); ++currentTag) {
							if (*currentTag == &universe) continue;
							if (priorityMap[(*currentTag)->getSchemaField()] < priorityOfTagh) {
								priorityOfTagh = priorityMap[(*currentTag)->getSchemaField()];
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
	for (std::vector<std::pair<Triangulation::Face_handle, PolygonHandle *> >::iterator currentFace = facesToRepair.begin(); currentFace != facesToRepair.end(); ++currentFace) {
		currentFace->first->info().removeAllTags();
		currentFace->first->info().addTag(currentFace->second);
	}
	
	return true;
}

bool IOWorker::matchSchemata(Triangulation &triangulation) {
	
	std::map<Field *, PolygonHandle *, FieldComparator> fieldMatch;
	std::map<PolygonHandle *, PolygonHandle *> equivalencies;
	
	// Find equivalencies
	for (std::vector<PolygonHandle *>::iterator currentPolygon = polygons.begin(); currentPolygon != polygons.end(); ++currentPolygon) {
		if (fieldMatch.count((*currentPolygon)->getSchemaField()) == 0) {
			fieldMatch[(*currentPolygon)->getSchemaField()] = *currentPolygon;
			equivalencies[*currentPolygon] = *currentPolygon;
		} else equivalencies[*currentPolygon] = fieldMatch[(*currentPolygon)->getSchemaField()];
	} std::cout << fieldMatch.size() << " classes found." << std::endl;
	
	// Re-tag
	for (Triangulation::Finite_faces_iterator currentFace = triangulation.finite_faces_begin(); currentFace != triangulation.finite_faces_end(); ++currentFace) {
		currentFace->info().substituteTagsWith(equivalencies[currentFace->info().getTags()]);
	}
	
	return true;
}

void IOWorker::removeConstraints(Triangulation &triangulation) {
	// Remove constrained edges that have the same polygon on both sides
  unsigned long long int constrainedEdgesRemoved = 0;
  for (Triangulation::All_edges_iterator currentEdge = triangulation.all_edges_begin(); currentEdge != triangulation.all_edges_end(); ++currentEdge) {
    if (!triangulation.is_constrained(*currentEdge)) continue;
    if (currentEdge->first->info().getOneTag() == currentEdge->first->neighbor(currentEdge->second)->info().getOneTag()) {
      triangulation.remove_constrained_edge(currentEdge->first, currentEdge->second);
      ++constrainedEdgesRemoved;
    }
  } std::cout << "\tRemoved " << constrainedEdgesRemoved << " constrained edges" << std::endl;
}

void IOWorker::removeVertices(Triangulation &triangulation) {
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
  } std::cout << "\tRemoved " << surroundedVerticesRemoved << " surrounded vertices" << std::endl;
  
  std::cout << "\tAfter: " << triangulation.number_of_faces() << " triangles in the triangulation" << std::endl;
}

bool IOWorker::reconstructPolygons(Triangulation &triangulation, std::vector<std::pair<PolygonHandle *, Polygon> > &outputPolygons) {
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
            if (repeatedVertices.count(newChain->front()) > 0) vertexChainMap[newChain->front()] = newChain;
            chainsStack.push(newChain);
          }
        } newChain = new std::list<Triangulation::Vertex_handle>();
      } newChain->push_back(*currentVertex);
    }
    
    // Final ring
    while (chainsStack.size() > 0) {
      newChain->splice(newChain->begin(), *chainsStack.top());
      chainsStack.pop();
    }
    
    // Degenerate (insufficient vertices to be valid)
    if (newChain->size() < 3) {
      delete newChain;
    } else {
      std::list<Triangulation::Vertex_handle>::iterator secondElement = newChain->begin();
      ++secondElement;
      // Degenerate (zero area)
      if (newChain->back() == *secondElement) delete newChain;
      // Valid
      else rings.push_back(newChain);
    }
    
    if (chainsStack.size() > 0) std::cout << "Error: Stack has " << chainsStack.size() << " elements. Should be empty." << std::endl;
    
    // STEP 5: Make a polygon from this list and save it
    std::vector<Ring> innerRings;
    Ring outerRing;
    for (std::list<std::list<Triangulation::Vertex_handle> *>::iterator currentRing = rings.begin(); currentRing != rings.end(); ++currentRing) {
      Ring newRing;
      for (std::list<Triangulation::Vertex_handle>::iterator currentPoint = (*currentRing)->begin(); currentPoint != (*currentRing)->end(); ++currentPoint) {
        newRing.push_back((*currentPoint)->point());
      } if (newRing.is_clockwise_oriented()) outerRing = newRing;
      else innerRings.push_back(newRing);
    } outputPolygons.push_back(std::pair<PolygonHandle *, Polygon>(currentTag, Polygon(outerRing, innerRings.begin(), innerRings.end())));
    // Free memory from the chains
    for (std::list<std::list<Triangulation::Vertex_handle> *>::iterator currentRing = rings.begin(); currentRing != rings.end(); ++currentRing) {
      delete *currentRing;
    }
  }
	
	return true;
}

bool IOWorker::exportPolygons(std::vector<std::pair<PolygonHandle *, Polygon> > &outputPolygons, const char *file, bool withProvenance) {
	
	// Prepare file
	const char *driverName = "ESRI Shapefile";
	OGRSFDriver *driver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(driverName);
	if (driver == NULL) {
		std::cout << "\tError: OGR Shapefile driver not found." << std::endl;
		return false;
	}
	
	OGRDataSource *dataSource = driver->Open(file, false);
	if (dataSource != NULL) {
		std::cout << "\tOverwriting file..." << std::endl;
		if (driver->DeleteDataSource(dataSource->GetName())!= OGRERR_NONE) {
			std::cout << "\tError: Couldn't erase file with same name." << std::endl;
			return false;
		} OGRDataSource::DestroyDataSource(dataSource);
	}
	
	std::cout << "\tWriting file... " << std::endl;
	dataSource = driver->CreateDataSource(file, NULL);
	if (dataSource == NULL) {
		std::cout << "\tError: Could not create file." << std::endl;
		return false;
	}
  
	OGRLayer *layer = dataSource->CreateLayer("polygons", spatialReference, wkbPolygon, NULL);
	if (layer == NULL) {
		std::cout << "\tError: Could not create layer." << std::endl;
		return false;
	}
	
	// Set up the fields that there will be
	if (withProvenance) {
		unsigned int longest = 0;
		for (std::vector<char *>::iterator currentFileName = fileNames.begin(); currentFileName != fileNames.end(); ++currentFileName) {
			if (strlen(*currentFileName) > longest) longest = strlen(*currentFileName);
		} OGRFieldDefn filenameField("File", OFTString);
		filenameField.SetWidth(longest);
		if (layer->CreateField(&filenameField) != OGRERR_NONE) {
			std::cout << "\tError: Could not create field File." << std::endl;
			return false;
		} OGRFieldDefn layerField("Layer", OFTInteger);
		if (layer->CreateField(&layerField) != OGRERR_NONE) {
			std::cout << "\tError: Could not create field Layer." << std::endl;
			return false;
		}
	}
	
	for (std::vector<FieldDefinition *>::iterator currentField = fields.begin(); currentField != fields.end(); ++currentField) {
		OGRFieldDefn newField((*currentField)->name, (*currentField)->type);
		newField.SetJustify((*currentField)->justification);
		newField.SetWidth((*currentField)->width);
		newField.SetPrecision((*currentField)->precision);
		if (layer->CreateField(&newField) != OGRERR_NONE) {
			std::cout << "\tError: Could not create field " << (*currentField)->name << "." << std::endl;
			return false;
		}
	}
	
	// Put fields in
	for (std::vector<std::pair<PolygonHandle *, Polygon> >::iterator currentPolygon = outputPolygons.begin(); currentPolygon != outputPolygons.end(); ++currentPolygon) {
		OGRPolygon polygon;
		OGRLinearRing outerRing;
    if (currentPolygon->second.outer_boundary().size() < 1) continue;
		for (Ring::Vertex_iterator currentVertex = currentPolygon->second.outer_boundary().vertices_begin();
         currentVertex != currentPolygon->second.outer_boundary().vertices_end();
         ++currentVertex) {
			outerRing.addPoint(CGAL::to_double(currentVertex->x()), CGAL::to_double(currentVertex->y()));
		} outerRing.addPoint(CGAL::to_double(currentPolygon->second.outer_boundary().vertex(0).x()), CGAL::to_double(currentPolygon->second.outer_boundary().vertex(0).y()));
		polygon.addRing(&outerRing);
		for (Polygon::Hole_const_iterator currentRing = currentPolygon->second.holes_begin(); currentRing != currentPolygon->second.holes_end(); ++currentRing) {
			OGRLinearRing innerRing;
			for (Ring::Vertex_iterator currentVertex = currentRing->vertices_begin(); currentVertex != currentRing->vertices_end(); ++currentVertex) {
				innerRing.addPoint(CGAL::to_double(currentVertex->x()), CGAL::to_double(currentVertex->y()));
			} innerRing.addPoint(CGAL::to_double(currentRing->vertex(0).x()), CGAL::to_double(currentRing->vertex(0).y()));
			polygon.addRing(&innerRing);
		} OGRFeature *feature = OGRFeature::CreateFeature(layer->GetLayerDefn());
		if (withProvenance) {
			feature->SetField("File", currentPolygon->first->getOriginalFile());
			feature->SetField("Layer", (int)currentPolygon->first->getLayer());
		} for (unsigned int currentField = 0; currentField < currentPolygon->first->getNumberOfFields(); currentField++) {
			switch (currentPolygon->first->getField(currentField)->getType()) {
				case OFTString:
					feature->SetField(fields[fieldEquivalencies[FieldDescriptor(currentPolygon->first->getOriginalFile(), currentPolygon->first->getLayer(), currentField)]]->name,
                            currentPolygon->first->getField(currentField)->getValueAsString());
					break;
				case OFTReal:
					feature->SetField(fields[fieldEquivalencies[FieldDescriptor(currentPolygon->first->getOriginalFile(), currentPolygon->first->getLayer(), currentField)]]->name,
                            currentPolygon->first->getField(currentField)->getValueAsDouble());
					break;
				case OFTInteger:
					feature->SetField(fields[fieldEquivalencies[FieldDescriptor(currentPolygon->first->getOriginalFile(), currentPolygon->first->getLayer(), currentField)]]->name,
                            currentPolygon->first->getField(currentField)->getValueAsInt());
					break;
				default:
					std::cout << "\tError: Type not implemented." << std::endl;
					break;
			}
		}
		
		// Put geometry in
		feature->SetGeometry(&polygon);
		
		// Create OGR feature
		if (layer->CreateFeature(feature) != OGRERR_NONE) std::cout << "\tError: Could not create feature." << std::endl;
		
		// Free OGR feature
		OGRFeature::DestroyFeature(feature);
	}
	
	// Free OGR data source
	OGRDataSource::DestroyDataSource(dataSource);
	
	return true;
}

bool IOWorker::exportTriangulation(Triangulation &t, const char *file, bool withNumberOfTags, bool withFields, bool withProvenance) {
	
	// Prepare file
	const char *driverName = "ESRI Shapefile";
	OGRSFDriver *driver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(driverName);
	if (driver == NULL) {
		std::cout << "Driver not found." << std::endl;
		return false;
	}
	
	OGRDataSource *dataSource = driver->Open(file, false);
	if (dataSource != NULL) {
		std::cout << "Erasing current file..." << std::endl;
		if (driver->DeleteDataSource(dataSource->GetName())!= OGRERR_NONE) {
			std::cout << "Couldn't erase current file." << std::endl;
			return false;
		} OGRDataSource::DestroyDataSource(dataSource);
	}
	
	std::cout << "Writing file... " << std::endl;
	dataSource = driver->CreateDataSource(file, NULL);
	if (dataSource == NULL) {
		std::cout << "Could not create file." << std::endl;
		return false;
	}
	
	OGRLayer *layer = dataSource->CreateLayer("triangles", spatialReference, wkbPolygon, NULL);
	if (layer == NULL) {
		std::cout << "Could not create layer." << std::endl;
		return false;
	}
	
	// Set up the fields that there will be
	if (withNumberOfTags) {
		OGRFieldDefn numberOfTagsField("Tags", OFTInteger);
		if (layer->CreateField(&numberOfTagsField) != OGRERR_NONE) {
			std::cout << "Could not create field Tags." << std::endl;
			return false;
		}
	}
	
	if (withProvenance) {
		unsigned int longest = 0;
		for (std::vector<char *>::iterator currentFileName = fileNames.begin(); currentFileName != fileNames.end(); ++currentFileName) {
			if (strlen(*currentFileName) > longest) longest = strlen(*currentFileName);
		} OGRFieldDefn filenameField("File", OFTString);
		filenameField.SetWidth(longest);
		if (layer->CreateField(&filenameField) != OGRERR_NONE) {
			std::cout << "Could not create field Filename." << std::endl;
			return false;
		} OGRFieldDefn layerField("Layer", OFTInteger);
		if (layer->CreateField(&layerField) != OGRERR_NONE) {
			std::cout << "Could not create field Layer." << std::endl;
			return false;
		}
	}
	
	if (withFields) {
		for (std::vector<FieldDefinition *>::iterator currentField = fields.begin(); currentField != fields.end(); ++currentField) {
			OGRFieldDefn newField((*currentField)->name, (*currentField)->type);
			newField.SetJustify((*currentField)->justification);
			newField.SetWidth((*currentField)->width);
			newField.SetPrecision((*currentField)->precision);
			if (layer->CreateField(&newField) != OGRERR_NONE) {
				std::cout << "Could not create field " << (*currentField)->name << "." << std::endl;
				return false;
			}
		}
	}
	
	// Put fields in
	for (CDT::Finite_faces_iterator currentFace = t.finite_faces_begin(); currentFace != t.finite_faces_end(); ++currentFace) {
		OGRLinearRing ring;
		ring.addPoint(CGAL::to_double((*(*currentFace).vertex(0)).point().x()), CGAL::to_double((*(*currentFace).vertex(0)).point().y()), 0.0);
		ring.addPoint(CGAL::to_double((*(*currentFace).vertex(1)).point().x()), CGAL::to_double((*(*currentFace).vertex(1)).point().y()), 0.0);
		ring.addPoint(CGAL::to_double((*(*currentFace).vertex(2)).point().x()), CGAL::to_double((*(*currentFace).vertex(2)).point().y()), 0.0);
		ring.addPoint(CGAL::to_double((*(*currentFace).vertex(0)).point().x()), CGAL::to_double((*(*currentFace).vertex(0)).point().y()), 0.0);
		OGRPolygon polygon;
		polygon.addRing(&ring);
		OGRFeature *feature = OGRFeature::CreateFeature(layer->GetLayerDefn());
		if (withNumberOfTags) {
      
      // Put number of tags
			if ((*currentFace).info().getTags() == NULL) {
				feature->SetField("Tags", 0);
			} else {
				feature->SetField("Tags", (int)(*currentFace).info().numberOfTags());
			}
      
      // Put number of tags, the universe doesn't count
      if ((*currentFace).info().getTags() == NULL) {
				feature->SetField("Tags", 0);
			} else if ((*currentFace).info().getTags() != &universe) {
        feature->SetField("Tags", (int)(*currentFace).info().numberOfTags());
			} else {
        feature->SetField("Tags", 0);
      }
		} if (withProvenance) {
      if ((*currentFace).info().getTags() == NULL) {
        feature->SetField("File", "");
        feature->SetField("Layer", -1);
      } else {
        feature->SetField("File", (*currentFace).info().getTags()->getOriginalFile());
        feature->SetField("Layer", (int)(*currentFace).info().getTags()->getLayer());
      }
		} if (withFields && (*currentFace).info().getTags() != NULL) for (unsigned int currentField = 0; currentField < (*currentFace).info().getTags()->getNumberOfFields(); currentField++) {
			switch ((*currentFace).info().getTags()->getField(currentField)->getType()) {
				case OFTString:
          feature->SetField(fields[fieldEquivalencies[FieldDescriptor((*currentFace).info().getTags()->getOriginalFile(), (*currentFace).info().getTags()->getLayer(), currentField)]]->name,
                            (*currentFace).info().getTags()->getField(currentField)->getValueAsString());
				case OFTReal:
          feature->SetField(fields[fieldEquivalencies[FieldDescriptor((*currentFace).info().getTags()->getOriginalFile(), (*currentFace).info().getTags()->getLayer(), currentField)]]->name,
                            (*currentFace).info().getTags()->getField(currentField)->getValueAsDouble());
				case OFTInteger:
          feature->SetField(fields[fieldEquivalencies[FieldDescriptor((*currentFace).info().getTags()->getOriginalFile(), (*currentFace).info().getTags()->getLayer(), currentField)]]->name,
                            (*currentFace).info().getTags()->getField(currentField)->getValueAsInt());
					break;
				default:
					std::cout << "Error: Type not implemented." << std::endl;
					break;
			}
      
		}
		
		// Put geometry in
		feature->SetGeometry(&polygon);
		
		// Create OGR feature
		if (layer->CreateFeature(feature) != OGRERR_NONE) std::cout << "Could not create feature." << std::endl;
		
		// Free OGR feature
		OGRFeature::DestroyFeature(feature);
	}
	
	// Free OGR data source
	OGRDataSource::DestroyDataSource(dataSource);
	
	return true;
}

unsigned int IOWorker::removeDuplicateVertices(std::list<Point> &ring) {
	unsigned int removed = 0;
	ring.push_back(ring.front());
	std::list<Point>::iterator previousVertex = ring.begin();
	std::list<Point>::iterator nextVertex = previousVertex;
	++nextVertex;
	while (nextVertex != ring.end()) {
		if (*previousVertex == *nextVertex) {
			nextVertex = ring.erase(nextVertex);
			++removed;
		} else {
			++previousVertex;
			++nextVertex;
		}
	} ring.pop_back();
	return removed;
}

std::vector<Ring *> IOWorker::splitRing(Ring &ring) {
	std::vector<Ring *> outputRings;
	
	// STEP 1: Put the edges in a triangulation
	Triangulation ringTriangulation;
  startingSearchFaceInRing = Triangulation::Face_handle();
	for (Ring::Edge_const_iterator currentEdge = ring.edges_begin(); currentEdge != ring.edges_end(); ++currentEdge) {
		Triangulation::Vertex_handle source = ringTriangulation.insert(currentEdge->source(), startingSearchFaceInRing);
    startingSearchFaceInRing = ringTriangulation.incident_faces(source);
		Triangulation::Vertex_handle target = ringTriangulation.insert(currentEdge->target(), startingSearchFaceInRing);
		Triangulation::Face_handle correspondingFace;
		int correspondingVertex;
    // Remove identical degenerate edges
		if (ringTriangulation.is_edge(source, target, correspondingFace, correspondingVertex)) {
			if (ringTriangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(correspondingFace, correspondingVertex))) {
        //std::cout << "Removing duplicate constraint <" << *source << ", " << *target << ">" << std::endl;
				ringTriangulation.remove_constraint(source, target);
				continue;
			}
		} //std::cout << "Inserting constraint <" << *source << ", " << *target << ">" << std::endl;
    ringTriangulation.insert_constraint(source, target);
    startingSearchFaceInRing = ringTriangulation.incident_faces(target);
	}
	
	// Free space of the old ring
	ring.clear();
	
	// STEP 2: Remove degenerate edges (not identical, so not caught during creation)
	for (Triangulation::Subconstraint_iterator currentEdge = ringTriangulation.subconstraints_begin();
       currentEdge != ringTriangulation.subconstraints_end();
       ++currentEdge) {
    //std::cout << "Checking subconstraint: <" << *(currentEdge->first.first) << ", " << *(currentEdge->first.second) << ">: " << ringTriangulation.number_of_enclosing_constraints(currentEdge->first.first, currentEdge->first.second) << " enclosing constraints." << std::endl;
		// Subconstraint_iterator has a weird return value...
		if (ringTriangulation.number_of_enclosing_constraints(currentEdge->first.first, currentEdge->first.second) % 2 == 0) {
			Triangulation::Face_handle f;
			int i;
			ringTriangulation.is_edge(currentEdge->first.first, currentEdge->first.second, f, i);
      if (ringTriangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(f, i))) {
        //std::cout << "Removing constraint..." << std::endl;
        ringTriangulation.remove_constraint(currentEdge->first.first, currentEdge->first.second);
      } else {
        //std::cout << "Adding constraint..." << std::endl;
        ringTriangulation.insert_constraint(currentEdge->first.first, currentEdge->first.second);
      }
		}
	}
	
	// STEP 3: Tag triangles
	PolygonHandle interior, exterior;
	std::stack<Triangulation::Face_handle> interiorStack, exteriorStack;
	exteriorStack.push(ringTriangulation.infinite_face());
	tagStack(interiorStack, exteriorStack, &interior, &exterior);
	
	// STEP 4: Get chains representing boundaries
	Triangulation::Face_handle seedingFace;
	std::vector<std::list<Triangulation::Vertex_handle> *> verticesList;
  for (Triangulation::Finite_faces_iterator seedingFace = ringTriangulation.finite_faces_begin(); seedingFace != ringTriangulation.finite_faces_end(); ++seedingFace) {
		
    if (seedingFace->info().getTags() != &interior) continue;
    
    std::list<Triangulation::Vertex_handle> *vertices = new std::list<Triangulation::Vertex_handle>();
		seedingFace->info().setTags(NULL);
		
		if (seedingFace->neighbor(2)->info().getTags() == &interior) {
			seedingFace->neighbor(2)->info().setTags(NULL);
			std::list<Triangulation::Vertex_handle> *l2 = getBoundary(seedingFace->neighbor(2), seedingFace->neighbor(2)->index(seedingFace), &interior);
			vertices->splice(vertices->end(), *l2);
			delete l2;
		} vertices->push_back(seedingFace->vertex(0));
		if (seedingFace->neighbor(1)->info().getTags() == &interior) {
			seedingFace->neighbor(1)->info().setTags(NULL);
			std::list<Triangulation::Vertex_handle> *l1 = getBoundary(seedingFace->neighbor(1), seedingFace->neighbor(1)->index(seedingFace), &interior);
			vertices->splice(vertices->end(), *l1);
			delete l1;
		} vertices->push_back(seedingFace->vertex(2));
		if (seedingFace->neighbor(0)->info().getTags() == &interior) {
			seedingFace->neighbor(0)->info().setTags(NULL);
			std::list<Triangulation::Vertex_handle> *l0 = getBoundary(seedingFace->neighbor(0), seedingFace->neighbor(0)->index(seedingFace), &interior);
			vertices->splice(vertices->end(), *l0);
			delete l0;
		} vertices->push_back(seedingFace->vertex(1));
    
		verticesList.push_back(vertices);
	}
	
	// From now on, process each list...
	for (std::vector<std::list<Triangulation::Vertex_handle> *>::iterator currentVerticesList = verticesList.begin(); currentVerticesList != verticesList.end(); ++currentVerticesList) {
		
		// STEP 5: Find cutting vertices
		std::set<Triangulation::Vertex_handle> visitedVertices;
    std::set<Triangulation::Vertex_handle> repeatedVertices;
		for (std::list<Triangulation::Vertex_handle>::iterator currentVertex = (*currentVerticesList)->begin(); currentVertex != (*currentVerticesList)->end(); ++currentVertex) {
			if (!visitedVertices.insert(*currentVertex).second) repeatedVertices.insert(*currentVertex);
		} visitedVertices.clear();
		
		// STEP 6: Cut and join rings in the correct order
    std::list<std::list<Triangulation::Vertex_handle> *> rings;
    std::stack<std::list<Triangulation::Vertex_handle> *> chainsStack;
    std::map<Triangulation::Vertex_handle, std::list<Triangulation::Vertex_handle> *> vertexChainMap;
    std::list<Triangulation::Vertex_handle> *newChain = new std::list<Triangulation::Vertex_handle>();
    for (std::list<Triangulation::Vertex_handle>::iterator currentVertex = (*currentVerticesList)->begin(); currentVertex != (*currentVerticesList)->end(); ++currentVertex) {
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
            if (repeatedVertices.count(newChain->front()) > 0) vertexChainMap[newChain->front()] = newChain;
            chainsStack.push(newChain);
          }
        } newChain = new std::list<Triangulation::Vertex_handle>();
      } newChain->push_back(*currentVertex);
    }
    
    // Final ring
    while (chainsStack.size() > 0) {
      newChain->splice(newChain->begin(), *chainsStack.top());
      chainsStack.pop();
    }
    
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
    
    // STEP 7: Make rings from these lists and save them
		for (std::list<std::list<Triangulation::Vertex_handle> *>::iterator currentRing = rings.begin(); currentRing != rings.end(); ++currentRing) {
			Ring *newRing = new Ring();
      for (std::list<Triangulation::Vertex_handle>::iterator currentVertex = (*currentRing)->begin(); currentVertex != (*currentRing)->end(); ++currentVertex) {
				newRing->push_back((*currentVertex)->point());
			} outputRings.push_back(newRing);
		}
		
		// Free memory from the chains
		for (std::list<std::list<Triangulation::Vertex_handle> *>::iterator currentRing = rings.begin(); currentRing != rings.end(); ++currentRing) delete *currentRing;
		
		// Free memory from the vertices lists
		delete *currentVerticesList;
	}
	
	// Free memory from the triangulation
	ringTriangulation.clear();
	
	std::cout << "\tCreated " << outputRings.size() << " rings." << std::endl;
	return outputRings;
}

void IOWorker::testRings(std::vector<Ring *> &outerRings, std::vector<Ring *> &innerRings, std::vector<std::vector<Ring> > &classification, long fid) {
	Triangulation ringsTriangulation;
  startingSearchFaceInRing = Triangulation::Face_handle();
	std::vector<std::vector<Triangulation::Vertex_handle> > ringsToTag;
	std::vector<PolygonHandle *> tags;
	tags.reserve(outerRings.size());
	std::map<PolygonHandle *, unsigned int> tagsMap;
	
	// STEP 1: Put the edges from the outer rings in the triangulation and save them for tagging
	for (std::vector<Ring *>::iterator currentRing = outerRings.begin(); currentRing != outerRings.end(); ++currentRing) {
		ringsToTag.push_back(std::vector<Triangulation::Vertex_handle>());
		tags.push_back(new PolygonHandle());
		tagsMap[tags.back()] = tags.size()-1;
		for (Ring::Edge_const_iterator currentEdge = (*currentRing)->edges_begin(); currentEdge != (*currentRing)->edges_end(); ++currentEdge) {
			Triangulation::Vertex_handle source = ringsTriangulation.insert(currentEdge->source(), startingSearchFaceInRing);
      startingSearchFaceInRing = ringsTriangulation.incident_faces(source);
			Triangulation::Vertex_handle target = ringsTriangulation.insert(currentEdge->target(), startingSearchFaceInRing);
			ringsTriangulation.insert_constraint(source, target);
      startingSearchFaceInRing = ringsTriangulation.incident_faces(target);
			ringsToTag.back().push_back(source);
		} ringsToTag.back().push_back(ringsToTag.back().front());
	}
	
	// STEP 2: Tag triangles
	Triangulation::Face_handle currentFace;
  int incident;
	std::stack<Triangulation::Face_handle> stack;
	for (unsigned int currentRing = 0; currentRing < ringsToTag.size(); ++currentRing) {
		for (unsigned int currentVertex = 1; currentVertex < ringsToTag[currentRing].size(); ++currentVertex) {
			if (!ringsTriangulation.is_edge(ringsToTag[currentRing].at(currentVertex-1), ringsToTag[currentRing].at(currentVertex), currentFace, incident)) {
				std::cout << "\tError: Cannot find adjoining face to an edge from the edge list!" << std::endl;
			} stack.push(currentFace);
		} tagStack(stack, tags[currentRing]);
	}
	
	// Tag the universe
	stack.push(ringsTriangulation.infinite_face());
	tagStack(stack, &universe);
	
	// Free memory
	ringsToTag.clear();
	
	// STEP 3: Check where inner rings belong
	for (std::vector<Ring *>::iterator currentRing = innerRings.begin(); currentRing != innerRings.end(); ++currentRing) {
		std::set<unsigned int> addedTo;
		for (Ring::Edge_const_iterator currentEdge = (*currentRing)->edges_begin(); currentEdge != (*currentRing)->edges_end(); ++currentEdge) {
			Triangulation::Locate_type locateType;
			int locateIndex;
			Triangulation::Face_handle location = ringsTriangulation.locate(currentEdge->source(), locateType, locateIndex, startingSearchFaceInRing);
      startingSearchFaceInRing = location;
			if (locateType == Triangulation::FACE) {
				PolygonHandle *tag = location->info().getTags();
				if (tagsMap.count(tag) > 0) {
					unsigned int tagIndex = tagsMap[tag];
					if (addedTo.count(tagIndex) == 0) {
						addedTo.insert(tagIndex);
						classification[tagIndex].push_back(**currentRing);
					}
				} else {
					std::cout << "\tFeature #" << fid << ": inner boundary vertex outside outer boundary. Using other vertices..." << std::endl;
				}
			}
		} if (addedTo.size() < 1) std::cout << "\tFeature #" << fid << ": inner boundary cannot fit in any outer boundary. Skipped." << std::endl;
		else if (addedTo.size() > 1) std::cout << "\tFeature #" << fid << ": inner boundary fits in more than one OB. Added to all." << std::endl;
	}
}

void IOWorker::copyFields(OGRFeature *ogrfeature, PolygonHandle *handle) {
	Field *newField;
	for (int i = 0; i < ogrfeature->GetFieldCount(); i++) {
		switch (ogrfeature->GetFieldDefnRef(i)->GetType()) {
			case OFTString:
				newField = new StringField(ogrfeature->GetFieldAsString(i));
				break;
			case OFTReal:
				newField = new DoubleField(ogrfeature->GetFieldAsDouble(i));
				break;
			case OFTInteger:
				newField = new IntField(ogrfeature->GetFieldAsInteger(i));
				break;
			default:
				std::cout << "\tError: Field type not supported. Skipped." << std::endl;
				continue;
				break;
		} handle->addField(newField);
	}
}

void IOWorker::tagStack(std::stack<Triangulation::Face_handle> &positiveStack, std::stack<Triangulation::Face_handle> &negativeStack, PolygonHandle *positiveHandle, PolygonHandle *negativeHandle) {
  //std::cout << "tagStack() Infinite vertex at: "  << std::endl;
	while (!positiveStack.empty() || !negativeStack.empty()) {
    //std::cout << "positiveStack: " << positiveStack.size() << " negativeStack: " << negativeStack.size() << std::endl;
		if (positiveStack.empty()) {
			Triangulation::Face_handle currentFace = negativeStack.top();
      //std::cout << "Triangle: <" << *(currentFace->vertex(0)) << ", " << *(currentFace->vertex(1)) << ", " << *(currentFace->vertex(2)) << ">" << std::endl;
			negativeStack.pop();
			currentFace->info().setTags(negativeHandle);
			if (currentFace->is_constrained(0)) {
				if (currentFace->neighbor(0)->info().getTags() != positiveHandle) {
					currentFace->neighbor(0)->info().setTags(positiveHandle);
					positiveStack.push(currentFace->neighbor(0));
				}
			} else {
				if (currentFace->neighbor(0)->info().getTags() != negativeHandle) {
					currentFace->neighbor(0)->info().setTags(negativeHandle);
					negativeStack.push(currentFace->neighbor(0));
				}
			} if (currentFace->is_constrained(1)) {
				if (currentFace->neighbor(1)->info().getTags() != positiveHandle) {
					currentFace->neighbor(1)->info().setTags(positiveHandle);
					positiveStack.push(currentFace->neighbor(1));
				}
			} else {
				if (currentFace->neighbor(1)->info().getTags() != negativeHandle) {
					currentFace->neighbor(1)->info().setTags(negativeHandle);
					negativeStack.push(currentFace->neighbor(1));
				}
			} if (currentFace->is_constrained(2)) {
				if (currentFace->neighbor(2)->info().getTags() != positiveHandle) {
					currentFace->neighbor(2)->info().setTags(positiveHandle);
					positiveStack.push(currentFace->neighbor(2));
				}
			} else {
				if (currentFace->neighbor(2)->info().getTags() != negativeHandle) {
					currentFace->neighbor(2)->info().setTags(negativeHandle);
					negativeStack.push(currentFace->neighbor(2));
				}
			}
		} else {
			Triangulation::Face_handle currentFace = positiveStack.top();
      //std::cout << "Triangle: <" << *(currentFace->vertex(0)) << ", " << *(currentFace->vertex(1)) << ", " << *(currentFace->vertex(2)) << ">" << std::endl;
			positiveStack.pop();
			currentFace->info().setTags(positiveHandle);
			if (currentFace->is_constrained(0)) {
				if (currentFace->neighbor(0)->info().getTags() != negativeHandle) {
					currentFace->neighbor(0)->info().setTags(negativeHandle);
					negativeStack.push(currentFace->neighbor(0));
				}
			} else {
				if (currentFace->neighbor(0)->info().getTags() != positiveHandle) {
					currentFace->neighbor(0)->info().setTags(positiveHandle);
					positiveStack.push(currentFace->neighbor(0));
				}
			} if (currentFace->is_constrained(1)) {
				if (currentFace->neighbor(1)->info().getTags() != negativeHandle) {
					currentFace->neighbor(1)->info().setTags(negativeHandle);
					negativeStack.push(currentFace->neighbor(1));
				}
			} else {
				if (currentFace->neighbor(1)->info().getTags() != positiveHandle) {
					currentFace->neighbor(1)->info().setTags(positiveHandle);
					positiveStack.push(currentFace->neighbor(1));
				}
			} if (currentFace->is_constrained(2)) {
				if (currentFace->neighbor(2)->info().getTags() != negativeHandle) {
					currentFace->neighbor(2)->info().setTags(negativeHandle);
					negativeStack.push(currentFace->neighbor(2));
				}
			} else {
				if (currentFace->neighbor(2)->info().getTags() != positiveHandle) {
					currentFace->neighbor(2)->info().setTags(positiveHandle);
					positiveStack.push(currentFace->neighbor(2));
				}
			}
		}
	}
}

void IOWorker::tagStack(std::stack<Triangulation::Face_handle> &stack, PolygonHandle *handle) {
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

std::list<Triangulation::Vertex_handle> * IOWorker::getBoundary(Triangulation::Face_handle face, int edge, PolygonHandle *polygon) {
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

void IOWorker::addtoCount(std::map<PolygonHandle *, unsigned int> &count, PolygonHandle *ph) {
	if (ph == NULL) return;
	if (ph->isMultiPolygonHandle()) {
		MultiPolygonHandle *mph = static_cast<MultiPolygonHandle *>(ph);
		for (std::list<PolygonHandle *>::const_iterator currentPolygonHandle = mph->getHandles()->begin(); currentPolygonHandle != mph->getHandles()->end(); ++currentPolygonHandle) {
			if (count.count(*currentPolygonHandle)) count[*currentPolygonHandle]++;
			else count[*currentPolygonHandle] = 1;
		}
	} else {
		if (count.count(ph)) count[ph]++;
		else count[ph] = 1;
	}
}

void IOWorker::addToLength(std::map<PolygonHandle *, double> &lengths, PolygonHandle *ph, double length) {
	if (ph == NULL) return;
	if (ph->isMultiPolygonHandle()) {
		MultiPolygonHandle *mph = static_cast<MultiPolygonHandle *>(ph);
		for (std::list<PolygonHandle *>::const_iterator currentPolygonHandle = mph->getHandles()->begin(); currentPolygonHandle != mph->getHandles()->end(); ++currentPolygonHandle) {
			if (lengths.count(*currentPolygonHandle)) lengths[*currentPolygonHandle] += length;
			else lengths[*currentPolygonHandle] = length;
		}
	} else {
		if (lengths.count(ph)) lengths[ph] += length;
		else lengths[ph] = length;
	}
}

void IOWorker::insertToStream(std::ostream &ostr, OGRFeatureDefn *layerDefinition, unsigned int indentation, int schemaIndex) {
  ostr << std::endl;
	for (int currentField = 0; currentField < layerDefinition->GetFieldCount(); currentField++) {
		OGRFieldDefn *fieldDefinition = layerDefinition->GetFieldDefn(currentField);
    if (currentField == schemaIndex) ostr << ">";
		for (unsigned int i = 0; i <= indentation; i++) ostr << "\t";
		insertToStream(ostr, fieldDefinition->GetType());
		ostr << "\t" << fieldDefinition->GetNameRef() << std::endl;
	}
}

void IOWorker::insertToStream(std::ostream &ostr, const OGRFieldType &ft) {
	switch (ft) {
		case OFTInteger:
			ostr << "int        ";
			break;
		case OFTIntegerList:
			ostr << "int[]      ";
			break;
		case OFTReal:
			ostr << "double     ";
			break;
		case OFTRealList:
			ostr << "double[]   ";
			break;
		case OFTString:
			ostr << "string     ";
			break;
		case OFTStringList:
			ostr << "string[]   ";
			break;
		case OFTWideString:
			ostr << "deprecated ";
			break;
		case OFTWideStringList:
			ostr << "deprecated ";
			break;
		case OFTBinary:
			ostr << "binary data";
			break;
		case OFTDate:
			ostr << "date       ";
			break;
		case OFTTime:
			ostr << "time       ";
			break;
		case OFTDateTime:
			ostr << "date & time";
			break;
		default:
			ostr << "unknown    ";
			break;
	}
}

void IOWorker::insertToStream(std::ostream &ostr, const OGRwkbGeometryType &gt) {
	switch (gt) {
		case wkbUnknown:
			ostr << "unknown";
			break;
		case wkbPoint:
			ostr << "point";
			break;
		case wkbLineString:
			ostr << "line string";
			break;
		case wkbPolygon:
			ostr << "polygon";
			break;
		case wkbMultiPoint:
			ostr << "multi point";
			break;
		case wkbMultiLineString:
			ostr << "multi line string";
			break;
		case wkbMultiPolygon:
			ostr << "multi polygon";
			break;
		case wkbGeometryCollection:
			ostr << "geometry collection";
			break;
		case wkbNone:
			ostr << "none";
			break;
		case wkbLinearRing:
			ostr << "linear ring";
			break;
		case wkbPoint25D:
			ostr << "2.5D point";
			break;
		case wkbLineString25D:
			ostr << "2.5D line string";
			break;
		case wkbPolygon25D:
			ostr << "2.5D polygon";
			break;
		case wkbMultiPoint25D:
			ostr << "2.5D multi point";
			break;
		case wkbMultiLineString25D:
			ostr << "2.5D multi line string";
			break;
		case wkbMultiPolygon25D:
			ostr << "2.5D multi polygon";
			break;
		case wkbGeometryCollection25D:
			ostr << "2.5D geometry collection";
			break;
		default:
			ostr << "other";
			break;
	}
}

void IOWorker::insertTriangulationInfo(std::ostream &ostr, const Triangulation &t) {
	// Number of tags
	unsigned int untagged = 0, onetag = 0, multipletags = 0, total;
	for (Triangulation::Finite_faces_iterator currentFace = t.finite_faces_begin(); currentFace != t.finite_faces_end(); ++currentFace) {
		if ((*currentFace).info().hasNoTags()) untagged++;
		else if ((*currentFace).info().hasOneTag()) onetag++;
		else multipletags++;
	} total = onetag + multipletags + untagged;
  ostr << "\tHoles:    " << untagged << " triangles (" << 100.0*untagged/total << " %)" << std::endl <<
  "\tOk:       " << onetag << " triangles (" << 100.0*onetag/total << " %)" << std::endl <<
  "\tOverlaps: " << multipletags << " triangles (" << 100.0*multipletags/total << " %)" << std::endl;
	
	// Other info?
}