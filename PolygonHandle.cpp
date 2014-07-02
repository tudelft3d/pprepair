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

#include "PolygonHandle.h"


PolygonHandle::PolygonHandle(OGRFeature* f) {
  feature = f;
}

PolygonHandle::~PolygonHandle() {
  OGRFeature::DestroyFeature(feature);
//	for (unsigned int i = 0; i < fields.size(); ++i) {
//		delete fields[i];
//	}
  // TODO: clear OGR memory?
}

std::string PolygonHandle::getValueAttributeAsString(std::string attr) {
  int i = feature->GetFieldIndex(attr.c_str());
  if (i == -1)
    return "";
  else
    return feature->GetFieldAsString(i);
}


const bool PolygonHandle::isMultiPolygonHandle() {
	return false;
}

MultiPolygonHandle::MultiPolygonHandle(PolygonHandle *ph) {
	if (ph->isMultiPolygonHandle()) {
		for (std::list<PolygonHandle *>::iterator currentHandle = ((MultiPolygonHandle *)ph)->handles.begin();
         currentHandle != ((MultiPolygonHandle *)ph)->handles.end();
         ++currentHandle) {
			handles.push_back(*currentHandle);
		}
	}
  else if (ph != NULL) {
		handles.push_back(ph);
	}
}

MultiPolygonHandle::~MultiPolygonHandle() {
//	for (unsigned int i = 0; i < fields.size(); ++i) {
//		delete fields[i];
//	}
  // TODO: clear OGR memory here perhaps?
}

const bool MultiPolygonHandle::isMultiPolygonHandle() {
	return true;
}

bool MultiPolygonHandle::hasHandle(PolygonHandle *handle) {
	if (find(handles.begin(), handles.end(), handle) != handles.end()) return true;
	return false;
}

void MultiPolygonHandle::addHandle(PolygonHandle *handle) {
	if (handle->isMultiPolygonHandle()) {
		for (std::list<PolygonHandle *>::iterator currentHandle = static_cast<MultiPolygonHandle *>(handle)->handles.begin();
         currentHandle != static_cast<MultiPolygonHandle *>(handle)->handles.end();
         ++currentHandle) {
			handles.push_back(*currentHandle);
		}
	}
  else handles.push_back(handle);
}

const std::list<PolygonHandle *> *MultiPolygonHandle::getHandles() {
	return &handles;
}

unsigned int MultiPolygonHandle::numberOfHandles() {
	return handles.size();
}