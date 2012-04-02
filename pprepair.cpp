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

enum RepairMethod {
    VALIDATE_ONLY = -1,
    NUMBER_OF_NEIGHBOURS = 0,
    ABSOLUTE_MAJORITY = 1,
    LONGEST_BOUNDARY = 2,
    REGIONS_BY_LONGEST_BOUNDARY = 3,
    REGIONS_BY_RANDOM_NEIGHBOUR = 4,
    PRIORITY_LIST = 5,
    PRIORITY_LIST_EDGEMATCHING = 6
};

int main(int argc, const char *argv[]) {

    time_t startTime = time(NULL);
    PlanarPartition pp;
    bool processInOrder = false;
    
    std::list<std::pair<std::string, int> > inputFiles;
    std::string outputFile, outputFileWithProvenance, taggedTriangulationOutputFile, triangulationOutputFile, triangulationOutputFileWithProvenance;
    bool makeHolesValid = false, splitRegions = false, alsoUniverse = false;
    double splitRegionsRatio;
    std::list<std::pair<RepairMethod, std::string> > repairMethods;
    
    // Process help argument
    if (argc == 1 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
        std::cout << "=== pprepair Help ===" << std::endl;
        std::cout << "    Simple:   pprepair [options]" << std::endl;
        std::cout << "    Advanced: pprepair -p [processing steps in order]" << std::endl;
        std::cout << "    Example:  ./pprepair -i \"myInput.shp\" -o \"myOutput.shp\" -fix" << std::endl;
        std::cout << "== Basic options ==" << std::endl;
        std::cout << "    -i filename [schemaindex] Add this file to the triangulation using this schema index" << std::endl;
        std::cout << "    -o filename  Output the reconstructed polygons in this file" << std::endl;
        std::cout << "    -fix  Automagically repair (same as -rrlb -rrrn)" << std::endl;
        std::cout << "== Possible steps (in usual processing order) ==" << std::endl;
        std::cout << "    -i filename [schemaindex] Add this file to the triangulation using this schema index" << std::endl;
        std::cout << "    -t  Tag the triangulation" << std::endl;
        std::cout << "    -otnt filename  Output the tagged triangulation with the number of tags to this file" << std::endl;
        std::cout << "    -vh  Consider holes as valid" << std::endl;
        std::cout << "    -sr ratio  Split invalid regions at triangles with a higher aspect ratio than this" << std::endl;
        std::cout << "    -v  Validate" << std::endl;
        std::cout << "    -au  Allow removing invalid regions (where convenient)" << std::endl;
        std::cout << "    -rtnn  Repair triangles by assigning them to the neighbour present on most sides" << std::endl;
        std::cout << "    -rtam  Repair triangles by assigning them to a neighbour present on at least 2 sides" << std::endl;
        std::cout << "    -rtlb  Repair triangles by assigning them to the neighbour present along the longest part of their boundary" << std::endl;
        std::cout << "    -rrlb  Repair regions by assigning them to the neighbour present along the longest part of their boundary" << std::endl;
        std::cout << "    -rrrn  Repair regions by assigning them to a random neighbour" << std::endl;
        std::cout << "    -rpl filename  Repair by assigning according to the priority list in this file" << std::endl;
        std::cout << "    -rem filename  Repair for edge matching according to the priority list in this file" << std::endl;
        std::cout << "    -ot filename  Output the triangulation to this file" << std::endl;
        std::cout << "    -otwp filename  Output the triangulation to this file, including the input file where each triangle came from" << std::endl;
        std::cout << "    -rp  Reconstruct polygons" << std::endl;
        std::cout << "    -o filename  Output the reconstructed polygons in this file" << std::endl;
        std::cout << "    -owp filename  Output the reconstructed polygons in this file, including the input file where they came from" << std::endl;
        std::cout << "    -pi  Print triangulation information" << std::endl;
        return 0;
    }
    
    for (int argNum = 1; argNum < argc; ++argNum) {
        
        // Process in order
        if (strcmp(argv[argNum], "-p") == 0) {
            processInOrder = true;
        }
        
        // Input
        else if (strcmp(argv[argNum], "-i") == 0) {
            if (argNum + 2 <= argc - 1 && argv[argNum+2][0] != '-') {
                ++argNum;
                if (processInOrder) pp.addToTriangulation(argv[argNum], atoi(argv[argNum+1]));
                else inputFiles.push_back(std::pair<std::string, int>(argv[argNum], atoi(argv[argNum+1])));
                ++argNum;
            } else if (argNum + 1 <= argc - 1 && argv[argNum+1][0] != '-') {
                ++argNum;
                if (processInOrder) pp.addToTriangulation(argv[argNum]);
                else inputFiles.push_back(std::pair<std::string, int>(argv[argNum], 0));
            } else {
                std::cerr << "Error: Missing filename argument for -i";
                return 1;
            }
        }
        
        // Tag triangulation
        else if (strcmp(argv[argNum], "-fix") == 0) {
            if (!processInOrder) {
                repairMethods.push_back(std::pair<RepairMethod, std::string>(REGIONS_BY_LONGEST_BOUNDARY, std::string()));
                repairMethods.push_back(std::pair<RepairMethod, std::string>(REGIONS_BY_RANDOM_NEIGHBOUR, std::string()));
            }
        }
        
        // Tag triangulation
        else if (strcmp(argv[argNum], "-t") == 0) {
            if (processInOrder) pp.tagTriangulation();
        }
        
        // Output tagged triangulation
        else if (strcmp(argv[argNum], "-otnt") == 0) {
            if (argNum + 1 <= argc - 1 && argv[argNum+1][0] != '-') {
                ++argNum;
                if (processInOrder) pp.exportTriangulation(argv[argNum], true, false, false);
                else taggedTriangulationOutputFile = argv[argNum];
            } else {
                std::cerr << "Error: Missing filename argument for -otnt";
                return 1;
            }
        }
        
        // Make all holes valid
        else if (strcmp(argv[argNum], "-vh") == 0) {
            if (processInOrder) pp.makeAllHolesValid();
            else makeHolesValid = true;
        }
        
        // Split regions
        else if (strcmp(argv[argNum], "-sr") == 0) {
            if (argNum + 1 <= argc - 1 && argv[argNum+1][0] != '-') {
                ++argNum;
                if (processInOrder) pp.splitRegions(atof(argv[argNum]));
                else {
                    splitRegions = true;
                    splitRegionsRatio = atof(argv[argNum]);
                }
            } else {
                std::cerr << "Error: Missing ratio argument for -sr";
                return 1;
            }
        }
        
        // Validate
        else if (strcmp(argv[argNum], "-v") == 0) {
            if (processInOrder) pp.checkValidity();
            else repairMethods.push_back(std::pair<RepairMethod, std::string>(VALIDATE_ONLY, std::string()));
        }
        
        // Validate
        else if (strcmp(argv[argNum], "-au") == 0) {
            alsoUniverse = !alsoUniverse;
        }
        
        // Repair triangle by number of neighbours
        else if (strcmp(argv[argNum], "-rtnn") == 0) {
            if (processInOrder) pp.repairTrianglesByNumberOfNeighbours(alsoUniverse);
            else repairMethods.push_back(std::pair<RepairMethod, std::string>(NUMBER_OF_NEIGHBOURS, std::string()));
        }
        
        // Repair triangles by absolute majority
        else if (strcmp(argv[argNum], "-rtam") == 0) {
            if (processInOrder) pp.repairTrianglesByAbsoluteMajority(alsoUniverse);
            else repairMethods.push_back(std::pair<RepairMethod, std::string>(ABSOLUTE_MAJORITY, std::string()));
        }
        
        // Repair triangles by longest boundary
        else if (strcmp(argv[argNum], "-rtlb") == 0) {
            if (processInOrder) pp.repairTrianglesByLongestBoundary(alsoUniverse);
            else repairMethods.push_back(std::pair<RepairMethod, std::string>(LONGEST_BOUNDARY, std::string()));
        }
        
        // Repair regions by longest boundary
        else if (strcmp(argv[argNum], "-rrlb") == 0) {
            if (processInOrder) pp.repairRegionsByLongestBoundary(alsoUniverse);
            else repairMethods.push_back(std::pair<RepairMethod, std::string>(REGIONS_BY_LONGEST_BOUNDARY, std::string()));
        }
        
        // Repair regions by random neighbour
        else if (strcmp(argv[argNum], "-rrrn") == 0) {
            if (processInOrder) pp.repairRegionsByRandomNeighbour(alsoUniverse);
            else repairMethods.push_back(std::pair<RepairMethod, std::string>(REGIONS_BY_RANDOM_NEIGHBOUR, std::string()));
        }
        
        // Repair by priority list
        else if (strcmp(argv[argNum], "-rpl") == 0) {
            if (argNum + 1 <= argc - 1 && argv[argNum+1][0] != '-') {
                ++argNum;
                if (processInOrder) pp.repairByPriorityList(argv[argNum]);
                else repairMethods.push_back(std::pair<RepairMethod, std::string>(PRIORITY_LIST, argv[argNum]));
            } else {
                std::cerr << "Error: Missing priority list argument for -rpl";
                return 1;
            }
        }
        
        // Repair by priority list (edge matching)
        else if (strcmp(argv[argNum], "-rem") == 0) {
            if (argNum + 1 <= argc - 1 && argv[argNum+1][0] != '-') {
                ++argNum;
                if (processInOrder) pp.repairEdgeMatching(argv[argNum]);
                else repairMethods.push_back(std::pair<RepairMethod, std::string>(PRIORITY_LIST_EDGEMATCHING, argv[argNum]));
            } else {
                std::cerr << "Error: Missing priority list argument for -rem";
                return 1;
            }
        }
        
        // Output triangulation
        else if (strcmp(argv[argNum], "-ot") == 0) {
            if (argNum + 1 <= argc - 1 && argv[argNum+1][0] != '-') {
                ++argNum;
                if (processInOrder) pp.exportTriangulation(argv[argNum], false, true, false);
                else triangulationOutputFile = argv[argNum];
            } else {
                std::cerr << "Error: Missing filename argument for -ot";
                return 1;
            }
        }
        
        // Output triangulation with provenance
        else if (strcmp(argv[argNum], "-otwp") == 0) {
            if (argNum + 1 <= argc - 1 && argv[argNum+1][0] != '-') {
                ++argNum;
                if (processInOrder) pp.exportTriangulation(argv[argNum], false, true, true);
                else triangulationOutputFileWithProvenance = argv[argNum];
            } else {
                std::cerr << "Error: Missing filename argument for -otwp";
                return 1;
            }
        }
        
        // Reconstruct polygons
        else if (strcmp(argv[argNum], "-rp") == 0) {
            if (processInOrder) pp.reconstructPolygons();
        }
        
        // Output
        else if (strcmp(argv[argNum], "-o") == 0) {
            if (argNum + 1 <= argc - 1 && argv[argNum+1][0] != '-') {
                ++argNum;
                if (processInOrder) pp.exportPolygons(argv[argNum], false);
                else outputFile = argv[argNum];
            } else {
                std::cerr << "Error: Missing filename argument for -o";
                return 1;
            }
        }
        
        // Output with provenance
        else if (strcmp(argv[argNum], "-owp") == 0) {
            if (argNum + 1 <= argc - 1 && argv[argNum+1][0] != '-') {
                ++argNum;
                if (processInOrder) pp.exportPolygons(argv[argNum], true);
                else outputFileWithProvenance = argv[argNum];
            } else {
                std::cerr << "Error: Missing filename argument for -o";
                return 1;
            }
        }
        
        // Print triangulation information
        else if (strcmp(argv[argNum], "-pi") == 0) {
            if (processInOrder) pp.printInfo();
        }
        
        // Unrecognised option
        else {
            std::cerr << "Error: unrecognised option " << argv[argNum] << std::endl;
        }
    }
    
    // For the simple mode
    if (!processInOrder) {
        
        // Process input
        if (inputFiles.size() == 0) {
            std::cerr << "Error: No input files given.";
            return 1;
        } for (std::list<std::pair<std::string, int> >::iterator currentFile = inputFiles.begin(); currentFile != inputFiles.end(); ++currentFile)
            pp.addToTriangulation(currentFile->first.c_str(), currentFile->second);
        
        // Tag
        pp.tagTriangulation();
        
        // Print info
        std::cout << "Input triangulation:" << std::endl;
        pp.printInfo();
        
        // Output triangulation with number of tags
        if (taggedTriangulationOutputFile.size() > 0) pp.exportTriangulation(taggedTriangulationOutputFile.c_str(), true, false, false);
        
        // Consider holes as valid
        if (makeHolesValid) pp.makeAllHolesValid();
        
        // Split regions
        if (splitRegions) pp.splitRegions(splitRegionsRatio);
        
        // Repair
        bool outputResults = false;
        for (std::list<std::pair<RepairMethod, std::string> >::iterator currentFile = repairMethods.begin(); currentFile != repairMethods.end(); ++currentFile) {
            switch (currentFile->first) {
                case VALIDATE_ONLY:
                    pp.checkValidity();
                    break;
                    
                case NUMBER_OF_NEIGHBOURS:
                    pp.repairTrianglesByNumberOfNeighbours(alsoUniverse);
                    outputResults = true;
                    break;
                    
                case ABSOLUTE_MAJORITY:
                    pp.repairTrianglesByAbsoluteMajority(alsoUniverse);
                    outputResults = true;
                    break;
                    
                case LONGEST_BOUNDARY:
                    pp.repairTrianglesByLongestBoundary(alsoUniverse);
                    outputResults = true;
                    break;
                    
                case REGIONS_BY_LONGEST_BOUNDARY:
                    pp.repairRegionsByLongestBoundary(alsoUniverse);
                    outputResults = true;
                    break;
                    
                case REGIONS_BY_RANDOM_NEIGHBOUR:
                    pp.repairRegionsByRandomNeighbour(alsoUniverse);
                    outputResults = true;
                    break;
                    
                case PRIORITY_LIST:
                    pp.repairByPriorityList(currentFile->second.c_str());
                    outputResults = true;
                    break;
                    
                case PRIORITY_LIST_EDGEMATCHING:
                    pp.repairEdgeMatching(currentFile->second.c_str());
                    outputResults = true;
                    break;
                    
                default:
                    break;
            }
        }
        
        // Print info
        if (outputResults) {
            std::cout << "Repaired triangulation:" << std::endl;
            pp.printInfo();
        }
        
        // Output the triangulation
        if (triangulationOutputFile.size() > 0) pp.exportTriangulation(triangulationOutputFile.c_str(), false, true, false);
        
        // Output the triangulation with provenance
        if (triangulationOutputFileWithProvenance.size() > 0) pp.exportTriangulation(triangulationOutputFileWithProvenance.c_str(), false, true, true);
        
        // Reconstruct polygons
        if (outputFile.size() > 0 || outputFileWithProvenance.size() > 0) pp.reconstructPolygons();
        
        // Output
        if (outputFile.size() > 0) pp.exportPolygons(outputFile.c_str(), false);
        
        // Output with provenance
        if (outputFileWithProvenance.size() > 0) pp.exportPolygons(outputFileWithProvenance.c_str(), true);
    }
    
    time_t totalTime = time(NULL)-startTime;
	std::cout << "Done! Process finished in " << totalTime/60 << " minutes " << totalTime%60 << " seconds." << std::endl;
    
    return 0;
}

