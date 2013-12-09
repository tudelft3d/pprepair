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

#ifndef FACEINFO_H
#define FACEINFO_H

#include "PolygonHandle.h"

class FaceInfo {
public:
  // Constructors and destructors
	FaceInfo();
	~FaceInfo();
  
  // Clean (and expensive) access operations
	bool hasTag(PolygonHandle *handle);
  bool hasNoTags() const;
  bool hasOneTag() const;
  unsigned int numberOfTags() const;
  void addTag(PolygonHandle *handle);
  void removeAllTags();
  void substituteTagsWith(PolygonHandle *handle);
  PolygonHandle * getOneTag() const;
  
  // Dirty (and cheap) access operations
	PolygonHandle * getTags() const;
	void setTags(PolygonHandle *handle);
  
private:
	// Tags to the polygons it belongs to.
	// If more than one, it points to a MultiPolygonHandle with a set of pointers to PolygonHandles.
	PolygonHandle *tag;
};

#endif