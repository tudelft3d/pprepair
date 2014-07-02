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

#ifndef POLYGONHANDLE_H
#define POLYGONHANDLE_H

#include "definitions/definitions.h"

// Store a tag
class PolygonHandle {
public:
	PolygonHandle(OGRFeature* f = NULL);
	virtual ~PolygonHandle();
  
  OGRFeature* feature;
  std::string getValueAttributeAsString(std::string attr);
  
	virtual const bool isMultiPolygonHandle();
};

// Trick to save memory when multiple tags are not present
class MultiPolygonHandle : public PolygonHandle {
public:
	// Constructor. No need to start with less than two
  MultiPolygonHandle(PolygonHandle *ph);
	~MultiPolygonHandle();
	
  // Checking whether it's a MultiPolygonHandle
	virtual const bool isMultiPolygonHandle();
  
  // Access functions to the individual PolygonHandles
  bool hasHandle(PolygonHandle *handle);
  void addHandle(PolygonHandle *handle);
  const std::list<PolygonHandle *> *getHandles();
  unsigned int numberOfHandles();
private:
	std::list<PolygonHandle *> handles;
};

#endif