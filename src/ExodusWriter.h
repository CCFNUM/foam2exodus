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

    // Polyhedral ("unknown") cells cannot be written as a standard Exodus
    // element, so they are split into conformal tets/pyramids. A tri face
    // becomes a tet, a quad face a pyramid, and an n-gon face a fan of tets
    // around a shared per-face centroid; the apex is a per-cell centroid.
    struct SubElem
    {
        char type;               // 'T' tet, 'P' pyramid
        std::vector<int> nodes;  // 0-based indices into base + extra points
    };
    std::vector<Point> polyExtraPoints;  // appended after the reader points
    std::vector<SubElem> polySubElems;
    // Originating cell index of each sub-element, so decomposed elements can be
    // grouped into per-region (zone/mesh) blocks rather than one global block.
    std::vector<int> polySubElemCell;
    // Per cell: 1 if the cell is split into sub-elements instead of being
    // written as a standard element. True for "unknown" cells and for any
    // standard cell whose ordering would produce a non-positive volume.
    std::vector<char> cellDecomposed;

    // cellIdx -> block group name (zone name, "unzoned", or "fluid"), matching
    // how standard cells are grouped so poly blocks align with them.
    std::vector<std::string> buildCellGroups(
        const std::vector<Cell>& cells,
        const std::vector<CellZone>& cellZones) const;
    // 1-based Exodus element ID of each sub-element, filled while writing
    // connectivity and consumed by writeSideSets.
    std::vector<int> polySubElemExoId;
    // Boundary face index -> (sub-element index, local side id) entries so a
    // boundary face of a polyhedral cell resolves to the sub-element side(s)
    // that actually carry it (an n-gon face yields several entries).
    std::map<int, std::vector<std::pair<int, int>>> polyFaceToSubs;

    void buildPolyDecomposition(const std::vector<Point>& points,
                                const std::vector<Face>& faces,
                                const std::vector<Cell>& cells,
                                const std::vector<int>& owner,
                                int boundaryStart);

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
