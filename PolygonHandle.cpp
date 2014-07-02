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

PolygonHandle::PolygonHandle(OGRFeature* f) {
  feature = f;
}

PolygonHandle::~PolygonHandle() {
  OGRFeature::DestroyFeature(feature);
//	for (unsigned int i = 0; i < fields.size(); ++i) {
//		delete fields[i];
//	}
}

std::string PolygonHandle::getValueAttributeAsString(std::string attr) {
  int i = feature->GetFieldIndex(attr.c_str());
  if (i == -1)
    return "";
  else
    return feature->GetFieldAsString(i);
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
	}
  else if (ph != NULL) {
		handles.push_back(ph);
	}
}

MultiPolygonHandle::~MultiPolygonHandle() {
//	for (unsigned int i = 0; i < fields.size(); ++i) {
//		delete fields[i];
//	}
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