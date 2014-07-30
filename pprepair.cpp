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
#include <tclap/CmdLine.h>


int main (int argc, char* const argv[]) {
  
  std::vector<std::string> repairMethods;
  repairMethods.push_back("RN"); //-- random neighbour
  repairMethods.push_back("LB"); //-- longest boundary
  repairMethods.push_back("PL"); //-- priority list
  repairMethods.push_back("EM"); //-- edge-matching
  TCLAP::ValuesConstraint<std::string> rmVals(repairMethods);

  TCLAP::CmdLine cmd("Allowed options", ' ', "");
  try {
    TCLAP::MultiArg<std::string> inputDSs       ("i", "input", "input OGR dataset (this can be used more than once)", true, "string");
    TCLAP::ValueArg<std::string> extent         ("e", "extent", "spatial extent", false, "", "string");
    TCLAP::ValueArg<std::string> outfiles       ("o", "output", "folder for repaired file(s) (SHP only)", false, "","string");
    TCLAP::ValueArg<std::string> repair         ("r", "repair", "repair method used: RN/LB/PL/EM", false, "", &rmVals);
    TCLAP::SwitchArg             validation     ("v", "validation", "validation only (gaps and overlaps reported)", false);
    TCLAP::ValueArg<std::string> priority       ("",  "priority", "priority list for repairing", false, "", "string");
    TCLAP::ValueArg<std::string> outerrors      ("",  "outerrors", "errors to a shapefile", false, "","string");
    TCLAP::ValueArg<std::string> outtr          ("",  "outtr", "output triangulation to a shapefile", false, "","string");

    cmd.add(inputDSs);
    cmd.add(extent);
    cmd.add(outfiles);
    cmd.add(outerrors);
    cmd.add(outtr);
    cmd.add(priority);
    cmd.xorAdd(repair, validation);
    cmd.parse( argc, argv );
    
    //-- add input datasets to PP
    PlanarPartition pp;      
    std::vector<std::string> inputs = inputDSs.getValue();
    for (std::vector<std::string>::iterator it = inputs.begin() ; it != inputs.end(); ++it) {
      if (pp.addOGRdataset(*it) == false)
        throw false;
    }
    std::cout << "Total input polygons: " << pp.noPolygons() << std::endl;
    //-- add spatial extent
    if (extent.getValue() != "") {
      if (pp.addOGRdatasetExtent(extent.getValue()) == false)
        throw false;
    }

    //-- tag the triangulation
    pp.buildPP();
    
    //-- validation only
    if (validation.getValue() == true) {
      if (pp.isValid() == false) {
        std::cout << "\nValidation:\n\t planar partition is NOT valid." << std::endl;
        pp.printInfo();
      }
      else {
        std::cout << "\nValidation:\n\t planar partition is valid." << std::endl;
      }
    }
    else { //-- repairing
      pp.printInfo();
      
      if (repair.getValue() == "PL") {
        pp.repair("PL", true, priority.getValue());
      }
      else {
        pp.repair(repair.getValue());
      }
      //-- if there was a 'tie' then fix with RN
      if (pp.isValid() == false) {
        std::cout << "Reparing 'ties'..." << std::endl;
        pp.repair("RN");
      }
      pp.printInfo();
      
      //-- output repaired SHP files
      if (outfiles.getValue() != "") {
        pp.reconstructPolygons();
        if (pp.exportPolygonsSHP(outfiles.getValue()) == false) {
          return(0);
        }
      }
    }
    //-- output triangulation in SHP
    if (outtr.getValue() != "") {
      pp.exportTriangulation(outtr.getValue());
    }
	}
  catch (TCLAP::ArgException &e) {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    return(0);
  }
  catch (bool problems) {
    std::cerr << "Abort." << std::endl;
    return(0);
  }
  
  return(1);
}



