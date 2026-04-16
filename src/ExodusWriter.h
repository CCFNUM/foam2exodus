// File       : ExodusWriter.h
// Created    : Thu Mar 19 2026
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2026 CCFNUM, Lucerne University of Applied Sciences and
// Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef EXODUS_WRITER_H
#define EXODUS_WRITER_H

#include "OpenFOAMMeshReader.h"
#include <map>
#include <netcdf.h>
#include <string>
#include <vector>

// Forward declaration for MergedMeshReader
class MergedMeshReader;

class ExodusWriter
{
public:
    ExodusWriter(const std::string& filename);
    ~ExodusWriter();

    void writeMesh(const OpenFOAMMeshReader& reader);
    void writeMesh(const OpenFOAMMeshReader& reader,
                   const std::map<std::string, std::string>& elementBlockNames,
                   const std::map<std::string, std::string>& sidesetNames);

    // Overloads for MergedMeshReader
    void writeMesh(const MergedMeshReader& reader);
    void writeMesh(const MergedMeshReader& reader,
                   const std::map<std::string, std::string>& elementBlockNames,
                   const std::map<std::string, std::string>& sidesetNames);

private:
    std::string filename;
    int ncid;
    std::map<std::string, std::string> customElementBlockNames;
    std::map<std::string, std::string> customSidesetNames;

    // Maps OpenFOAM cell index to 1-based Exodus element ID. Populated by
    // writeElements as cells are assigned to blocks (which reorders them by
    // zone/type), and consumed by writeSideSets so that sideset entries refer
    // to the correct Exodus element.
    std::vector<int> cellToExodusElem;

    void initializeExodusFile(int numNodes,
                              int numElems,
                              int numElemBlocks,
                              int numNodeSets,
                              int numSideSets);
    void writeNodes(const std::vector<Point>& points);

    // Template method to work with any reader type (OpenFOAMMeshReader or
    // MergedMeshReader)
    template <typename ReaderType>
    void writeElementsImpl(const ReaderType& reader);

    template <typename ReaderType>
    void writeSideSetsImpl(const ReaderType& reader);

    void writeElements(const OpenFOAMMeshReader& reader);
    void writeSideSets(const OpenFOAMMeshReader& reader);
    void writeElements(const MergedMeshReader& reader);
    void writeSideSets(const MergedMeshReader& reader);

    void checkError(int status, const std::string& message);
    std::string getBlockName(const std::string& originalName);
    std::string getSidesetName(const std::string& originalName);
    std::vector<int> orderHexNodes(const Cell& cell,
                                   const std::vector<Face>& faces,
                                   const std::vector<Point>& points);
    std::vector<int> orderTetNodes(const Cell& cell,
                                   const std::vector<Face>& faces,
                                   const std::vector<Point>& points);
    std::vector<int> orderPyramidNodes(const Cell& cell,
                                       const std::vector<Face>& faces,
                                       const std::vector<Point>& points);
    std::vector<int> orderWedgeNodes(const Cell& cell,
                                     const std::vector<Face>& faces,
                                     const std::vector<Point>& points);
    int getHexFaceId(const std::vector<int>& faceNodes,
                     const std::vector<int>& hexNodes);
    int getTetFaceId(const std::vector<int>& faceNodes,
                     const std::vector<int>& tetNodes);
    int getPyramidFaceId(const std::vector<int>& faceNodes,
                         const std::vector<int>& pyrNodes);
    int getWedgeFaceId(const std::vector<int>& faceNodes,
                       const std::vector<int>& wedgeNodes);
};

#endif
