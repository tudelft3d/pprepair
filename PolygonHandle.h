/*
 Copyright (c) 2009-2022,
 Ken Arroyo Ohori    k.ohori@tudelft.nl
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

// Abstract class to store different types of attributes
class Field {
public:
  // Destructor
  virtual ~Field() = 0;
  
	// Find the type of field
	const virtual OGRFieldType getType() = 0;
	
	// Comparison
	bool operator<(Field &f);
	bool operator==(Field &f);
	
	// Getters
	virtual const char * getValueAsString();
	virtual double getValueAsDouble();
	virtual int getValueAsInt();
	
	// Setters
	virtual void setValueFromString(const char *v);
	virtual void setValueFromDouble(double v);
	virtual void setValueFromInt(int v);
};

// Specific classes
class StringField : public Field {
public:
	StringField(const char *v);
	~StringField();
	
	bool operator<(const StringField &f);
	bool operator==(const StringField &f);
	
	const OGRFieldType getType();
	const char * getValueAsString();
	void setValueFromString(const char *v);
private:
	char * contents;
};
class DoubleField : public Field {
public:
	DoubleField(double v);
	
	bool operator<(const DoubleField &f);
	bool operator==(const DoubleField &f);
	
	const OGRFieldType getType();
	double getValueAsDouble();
	void setValueFromDouble(double v);
private:
	double contents;
};
class IntField : public Field {
public:
	IntField(int v);
	
	bool operator<(const IntField &f);
	bool operator==(const IntField &f);
	
	const OGRFieldType getType();
	int getValueAsInt();
	void setValueFromInt(int v);
private:
	int contents;
};

// Store a tag
class PolygonHandle {
public:
	// Constructors and destructors
	PolygonHandle(unsigned int si = 0, char *of = NULL, unsigned int l = 0, long fid = 0);
	virtual ~PolygonHandle();
  
  // References
	char * getOriginalFile();
	unsigned int getLayer();
  
  // Field information
	void addField(Field *field);
  Field * getSchemaField();
  Field * getField(unsigned int i);
	unsigned long getNumberOfFields();
  
  // Checking whether it's a MultiPolygonHandle
	virtual const bool isMultiPolygonHandle();
protected:
	char *originalFile;
	unsigned int layer;
	//long featureID;
	
	// The field to use as schema
	// (could be changed to a set or regex to represent complex criteria)
	unsigned int schemaIndex;
	
	// Fields it contains
	std::vector<Field *> fields;
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
  unsigned long numberOfHandles();
private:
	std::list<PolygonHandle *> handles;
};

#endif
