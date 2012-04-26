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

#ifndef POLYGONHANDLE_H
#define POLYGONHANDLE_H

#include "definitions/Definitions.h"

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
	unsigned int getNumberOfFields();
    
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
    unsigned int numberOfHandles();
private:
	std::list<PolygonHandle *> handles;
};

#endif