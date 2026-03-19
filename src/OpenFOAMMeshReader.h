// File       : OpenFOAMMeshReader.h
// Created    : Thu Mar 19 2026
// Author     : Mhamad Mahdi Alloush
// Description: 
// Copyright (c) 2026 CCFNUM, Lucerne University of Applied Sciences and
// Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef OPENFOAM_MESH_READER_H
#define OPENFOAM_MESH_READER_H

#include <string>
#include <vector>
#include <map>
#include <array>

struct Point {
    double x, y, z;
};

struct Face {
    std::vector<int> pointIndices;
};

struct Cell {
    std::vector<int> faceIndices;
    std::string type;
};

struct BoundaryPatch {
    std::string name;
    std::string type;
    int startFace;
    int nFaces;
};

struct CellZone {
    std::string name;
    std::vector<int> cellIndices;
};

class OpenFOAMMeshReader {
public:
    OpenFOAMMeshReader(const std::string& casePath);
    
    void readMesh();
    
    int getNumPoints() const { return points.size(); }
    int getNumCells() const { return cells.size(); }
    int getNumBoundaryPatches() const { return boundaryPatches.size(); }
    int getNumCellZones() const { return cellZones.size(); }
    
    const std::vector<Point>& getPoints() const { return points; }
    const std::vector<Face>& getFaces() const { return faces; }
    const std::vector<Cell>& getCells() const { return cells; }
    const std::vector<BoundaryPatch>& getBoundaryPatches() const { return boundaryPatches; }
    const std::vector<CellZone>& getCellZones() const { return cellZones; }
    const std::vector<int>& getOwner() const { return owner; }
    const std::vector<int>& getNeighbour() const { return neighbour; }

private:
    std::string meshPath;
    std::vector<Point> points;
    std::vector<Face> faces;
    std::vector<Cell> cells;
    std::vector<BoundaryPatch> boundaryPatches;
    std::vector<CellZone> cellZones;
    std::vector<int> owner;
    std::vector<int> neighbour;
    
    void readPoints();
    void readFaces();
    void readOwner();
    void readNeighbour();
    void readBoundary();
    void readCellZones();
    void constructCells();
    
    std::string stripComments(const std::string& line);
    bool skipHeaderAndDetectFormat(std::ifstream& file);
    int parseLeadingInt(const std::string& line, size_t* endPos = nullptr);
    std::vector<int> parseIntListFromString(const std::string& line);
    
    // Binary reading helpers
    template<typename T>
    T readBinary(std::ifstream& file);
    void readPointsBinary(std::ifstream& file, int nPoints);
    void readFacesBinary(std::ifstream& file, int nFaces);
    void readOwnerBinary(std::ifstream& file, int nOwner);
    void readNeighbourBinary(std::ifstream& file, int nNeighbour);
};

#endif
