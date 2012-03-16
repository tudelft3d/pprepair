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

#include "PolygonHandle.h"

Field::~Field() {
    
}

bool Field::operator<(Field &f) {
	const Field *fp = &f;
	Field *tp = this;
	if (f.getType() != getType()) return fp < tp;
	switch (getType()) {
		case OFTString:
			return (*(StringField *)fp) < (*(StringField *)tp);
			break;
		case OFTReal:
			return (*(DoubleField *)fp) < (*(DoubleField *)tp);
			break;
		case OFTInteger:
			return (*(IntField *)fp) < (*(IntField *)tp);
			break;
		default:
			return false;
			break;
	}
}

bool Field::operator==(Field &f) {
	const Field *fp = &f;
	Field *tp = this;
	if (f.getType() != getType()) return fp == tp;
	switch (getType()) {
		case OFTString:
			return (*(StringField *)fp) == (*(StringField *)tp);
			break;
		case OFTReal:
			return (*(DoubleField *)fp) == (*(DoubleField *)tp);
			break;
		case OFTInteger:
			return (*(IntField *)fp) == (*(IntField *)tp);
			break;
		default:
			return false;
			break;
	}
}

const char * Field::getValueAsString() {
	std::cout << "Error: Getting value from abstract base class." << std::endl;
	return "";
}

double Field::getValueAsDouble() {
	std::cout << "Error: Getting value from abstract base class." << std::endl;
	return 0.0;
}

int Field::getValueAsInt() {
	std::cout << "Error: Getting value from abstract base class." << std::endl;
	return 0;
}

void Field::setValueFromString(const char *v) {
	std::cout << "Error: Setting value to abstract base class." << std::endl;
}

void Field::setValueFromDouble(double v) {
	std::cout << "Error: Setting value to abstract base class." << std::endl;
}

void Field::setValueFromInt(int v) {
	std::cout << "Error: Setting value to abstract base class." << std::endl;
}

StringField::StringField(const char *v) {
	setValueFromString(v);
}

StringField::~StringField() {
	free(contents);
}

bool StringField::operator<(const StringField &f) {
	return strcmp(f.contents, contents) < 0;
}

bool StringField::operator==(const StringField &f) {
	return strcmp(f.contents, contents) == 0;
}

const OGRFieldType StringField::getType() {
	return OFTString;
}

const char * StringField::getValueAsString() {
	return contents;
}

void StringField::setValueFromString(const char *v) {
	contents = new char[(strlen(v)+1)];
	strcpy(contents, v);
}

DoubleField::DoubleField(double v) {
	setValueFromDouble(v);
}

bool DoubleField::operator<(const DoubleField &f) {
	return f.contents < contents;
}

bool DoubleField::operator==(const DoubleField &f) {
	return f.contents == contents;
}

const OGRFieldType DoubleField::getType() {
	return OFTReal;
}

double DoubleField::getValueAsDouble() {
	return contents;
}

void DoubleField::setValueFromDouble(double v) {
	contents = v;
}

IntField::IntField(int v) {
	setValueFromInt(v);
}

bool IntField::operator<(const IntField &f) {
	return f.contents < contents;
}

bool IntField::operator==(const IntField &f) {
	return f.contents == contents;
}

const OGRFieldType IntField::getType() {
	return OFTInteger;
}

int IntField::getValueAsInt() {
	return contents;
}

void IntField::setValueFromInt(int v) {
	contents = v;
}

PolygonHandle::PolygonHandle(unsigned int si, char *of, unsigned int l, long fid) {
	schemaIndex = si;
	originalFile = of;
	layer = l;
}

PolygonHandle::~PolygonHandle() {
	for (unsigned int i = 0; i < fields.size(); ++i) {
		delete fields[i];
	}
}

void PolygonHandle::addField(Field *field) {
	fields.push_back(field);
}

Field * PolygonHandle::getSchemaField() {
	return fields[schemaIndex];
}

Field * PolygonHandle::getField(unsigned int i) {
	return fields[i];
}

unsigned int PolygonHandle::getNumberOfFields() {
	return fields.size();
}

char * PolygonHandle::getOriginalFile() {
	return originalFile;
}

unsigned int PolygonHandle::getLayer() {
	return layer;
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
	} else if (ph != NULL) {
		handles.push_back(ph);
	}
}

MultiPolygonHandle::~MultiPolygonHandle() {
	for (unsigned int i = 0; i < fields.size(); ++i) {
		delete fields[i];
	}
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
	} else handles.push_back(handle);
}

const std::list<PolygonHandle *> *MultiPolygonHandle::getHandles() {
	return &handles;
}

unsigned int MultiPolygonHandle::numberOfHandles() {
	return handles.size();
}