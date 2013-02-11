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

#include "PlanarPartition.h"

PlanarPartition::PlanarPartition() {
	// Registers drivers for all supported formats in OGR
	OGRRegisterAll();
	
	// Set internal states
	state = CREATED;
	
	// std::cout precision (for debugging)
	std::cout.setf(std::ios::fixed,std::ios::floatfield);
	std::cout.precision(8);
}

PlanarPartition::~PlanarPartition() {
	triangulation.clear();
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

void PlanarPartition::printInfo() {
    io.insertTriangulationInfo(std::cout, triangulation);
}