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
	if (tag == NULL) return 0;
	if (tag->isMultiPolygonHandle()) return static_cast<MultiPolygonHandle *>(tag)->numberOfHandles();
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