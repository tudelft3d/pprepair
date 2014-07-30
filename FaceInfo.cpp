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

#include "FaceInfo.h"

FaceInfo::FaceInfo() {
	tag = NULL;
}

FaceInfo::~FaceInfo() {
	if (tag != NULL) {
		if (tag->isMultiPolygonHandle()) {
			delete tag;
		}
	}
}

bool FaceInfo::hasTag(PolygonHandle *handle) {
	if (tag == NULL) return false;
	if (tag->isMultiPolygonHandle()) return static_cast<MultiPolygonHandle *>(tag)->hasHandle(handle);
	if (tag == handle) return true;
	return false;
}

bool FaceInfo::hasNoTags() const {
	if (tag == NULL) return true;
	return false;
}

bool FaceInfo::hasOneTag() const {
	if (tag == NULL) return false;
	if (tag->isMultiPolygonHandle()) return false;
	return true;
}

unsigned int FaceInfo::numberOfTags() const {
	if (tag == NULL)
    return 0;
	if (tag->isMultiPolygonHandle())
    return static_cast<MultiPolygonHandle *>(tag)->numberOfHandles();
	return 1;
}

void FaceInfo::addTag(PolygonHandle *handle) {
	if (tag == NULL) tag = handle;
	else if (tag->isMultiPolygonHandle()) static_cast<MultiPolygonHandle *>(tag)->addHandle(handle);
	else if (tag != handle) {
		MultiPolygonHandle *multiTag = new MultiPolygonHandle(tag);
		multiTag->addHandle(handle);
		tag = multiTag;
	}
}

void FaceInfo::removeAllTags() {
	if (tag == NULL) return;
	if (tag->isMultiPolygonHandle()) delete tag;
	tag = NULL;
}

void FaceInfo::substituteTagsWith(PolygonHandle *handle) {
	removeAllTags();
	addTag(handle);
}

PolygonHandle * FaceInfo::getOneTag() const {
  if (tag == NULL) return NULL;
	if (tag->isMultiPolygonHandle()) {
		return *static_cast<MultiPolygonHandle *>(tag)->getHandles()->begin();
	} return tag;
}

PolygonHandle * FaceInfo::getTags() const {
	return tag;
}

void FaceInfo::setTags(PolygonHandle *handle) {
	tag = handle;
}