// File       : MergedMeshReader.cpp
// Created    : Thu Mar 19 2026
// Author     : Mhamad Mahdi Alloush
// Description: 
// Copyright (c) 2026 CCFNUM, Lucerne University of Applied Sciences and
// Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "MergedMeshReader.h"
#include <iostream>
#include <stdexcept>
#include <filesystem>

MergedMeshReader::MergedMeshReader(const std::vector<std::string>& casePaths)
    : casePaths(casePaths) {
    if (casePaths.empty()) {
        throw std::runtime_error("No case paths provided for merged mesh reader");
    }
}

void MergedMeshReader::readMeshes() {
    std::cout << "Reading " << casePaths.size() << " OpenFOAM meshes for merging..." << std::endl;

    // Read all meshes
    for (size_t i = 0; i < casePaths.size(); ++i) {
        std::cout << "  [" << (i + 1) << "/" << casePaths.size() << "] Reading: "
                  << casePaths[i] << std::endl;

        OpenFOAMMeshReader* reader = new OpenFOAMMeshReader(casePaths[i]);
        reader->readMesh();

        std::cout << "    Points: " << reader->getNumPoints()
                  << ", Cells: " << reader->getNumCells()
                  << ", Patches: " << reader->getNumBoundaryPatches() << std::endl;

        readers.push_back(reader);
    }

    // Calculate offsets for each mesh
    offsets.resize(readers.size());
    int nodeOffset = 0;
    int faceOffset = 0;
    int cellOffset = 0;

    for (size_t i = 0; i < readers.size(); ++i) {
        offsets[i].nodeOffset = nodeOffset;
        offsets[i].faceOffset = faceOffset;
        offsets[i].cellOffset = cellOffset;

        nodeOffset += readers[i]->getNumPoints();
        faceOffset += readers[i]->getFaces().size();
        cellOffset += readers[i]->getNumCells();
    }

    // Merge all meshes
    mergeMeshes();

    std::cout << "Merged mesh statistics:" << std::endl;
    std::cout << "  Total Points: " << mergedPoints.size() << std::endl;
    std::cout << "  Total Cells: " << mergedCells.size() << std::endl;
    std::cout << "  Total Boundary Patches: " << mergedBoundaryPatches.size() << std::endl;
    std::cout << "  Total Cell Zones: " << mergedCellZones.size() << std::endl;
}

void MergedMeshReader::mergeMeshes() {
    mergePoints();
    mergeFaces();
    mergeCells();
    mergeBoundaryPatches();
    mergeCellZones();
    mergeOwnerNeighbour();
}

void MergedMeshReader::mergePoints() {
    // Simply concatenate all points
    for (size_t meshIdx = 0; meshIdx < readers.size(); ++meshIdx) {
        const auto& points = readers[meshIdx]->getPoints();
        mergedPoints.insert(mergedPoints.end(), points.begin(), points.end());
    }
}

void MergedMeshReader::mergeFaces() {
    // Merge faces with point index offsetting
    for (size_t meshIdx = 0; meshIdx < readers.size(); ++meshIdx) {
        const auto& faces = readers[meshIdx]->getFaces();
        int nodeOffset = offsets[meshIdx].nodeOffset;

        for (const auto& face : faces) {
            Face offsetFace;
            offsetFace.pointIndices.reserve(face.pointIndices.size());

            // Apply node offset to all point indices
            for (int pointIdx : face.pointIndices) {
                offsetFace.pointIndices.push_back(pointIdx + nodeOffset);
            }

            mergedFaces.push_back(offsetFace);
        }
    }
}

void MergedMeshReader::mergeCells() {
    // Merge cells with face index offsetting
    for (size_t meshIdx = 0; meshIdx < readers.size(); ++meshIdx) {
        const auto& cells = readers[meshIdx]->getCells();
        int faceOffset = offsets[meshIdx].faceOffset;

        for (const auto& cell : cells) {
            Cell offsetCell;
            offsetCell.type = cell.type;
            offsetCell.faceIndices.reserve(cell.faceIndices.size());

            // Apply face offset to all face indices
            for (int faceIdx : cell.faceIndices) {
                offsetCell.faceIndices.push_back(faceIdx + faceOffset);
            }

            mergedCells.push_back(offsetCell);
        }
    }
}

void MergedMeshReader::mergeBoundaryPatches() {
    // Track original boundary names to detect duplicates across meshes
    std::map<std::string, std::vector<std::string>> originalNameToMeshes;

    // First pass: collect all boundary names and their source meshes
    for (size_t meshIdx = 0; meshIdx < readers.size(); ++meshIdx) {
        const auto& patches = readers[meshIdx]->getBoundaryPatches();
        std::string meshPrefix = getMeshPrefix(meshIdx);

        for (const auto& patch : patches) {
            originalNameToMeshes[patch.name].push_back(meshPrefix);
        }
    }

    // Report any duplicate boundary names
    bool hasDuplicates = false;
    for (const auto& [name, meshes] : originalNameToMeshes) {
        if (meshes.size() > 1) {
            if (!hasDuplicates) {
                std::cout << "\nDetected duplicate boundary names across meshes:" << std::endl;
                hasDuplicates = true;
            }
            std::cout << "  '" << name << "' found in: ";
            for (size_t i = 0; i < meshes.size(); ++i) {
                std::cout << meshes[i];
                if (i < meshes.size() - 1) std::cout << ", ";
            }
            std::cout << std::endl;
        }
    }
    if (hasDuplicates) {
        std::cout << "  → Only duplicate names will be prefixed with mesh name" << std::endl;
        std::cout << std::endl;
    }

    // Merge boundary patches with face offset
    // Only prefix names that have duplicates across meshes
    for (size_t meshIdx = 0; meshIdx < readers.size(); ++meshIdx) {
        const auto& patches = readers[meshIdx]->getBoundaryPatches();
        int faceOffset = offsets[meshIdx].faceOffset;
        std::string meshPrefix = getMeshPrefix(meshIdx);

        for (const auto& patch : patches) {
            BoundaryPatch offsetPatch;

            // Only prefix if this name exists in multiple meshes
            if (originalNameToMeshes[patch.name].size() > 1) {
                offsetPatch.name = meshPrefix + "_" + patch.name;
            } else {
                offsetPatch.name = patch.name;
            }
            offsetPatch.type = patch.type;
            offsetPatch.startFace = patch.startFace + faceOffset;
            offsetPatch.nFaces = patch.nFaces;

            mergedBoundaryPatches.push_back(offsetPatch);

            // Add to sideset naming map for custom names
            sidesetNames[offsetPatch.name] = offsetPatch.name;
        }
    }
}

void MergedMeshReader::mergeCellZones() {
    // Track original zone names to detect duplicates across meshes
    std::map<std::string, std::vector<std::string>> originalZoneNameToMeshes;

    // First pass: collect all zone names and their source meshes
    for (size_t meshIdx = 0; meshIdx < readers.size(); ++meshIdx) {
        const auto& zones = readers[meshIdx]->getCellZones();
        std::string meshPrefix = getMeshPrefix(meshIdx);

        for (const auto& zone : zones) {
            originalZoneNameToMeshes[zone.name].push_back(meshPrefix);
        }
    }

    // Report any duplicate zone names
    bool hasDuplicates = false;
    for (const auto& [name, meshes] : originalZoneNameToMeshes) {
        if (meshes.size() > 1) {
            if (!hasDuplicates) {
                std::cout << "\nDetected duplicate cell zone names across meshes:" << std::endl;
                hasDuplicates = true;
            }
            std::cout << "  '" << name << "' found in: ";
            for (size_t i = 0; i < meshes.size(); ++i) {
                std::cout << meshes[i];
                if (i < meshes.size() - 1) std::cout << ", ";
            }
            std::cout << std::endl;
        }
    }
    if (hasDuplicates) {
        std::cout << "  → All zones prefixed with mesh name to ensure unique element blocks" << std::endl;
        std::cout << std::endl;
    }

    // Merge cell zones with unique naming and cell index offsetting
    for (size_t meshIdx = 0; meshIdx < readers.size(); ++meshIdx) {
        const auto& zones = readers[meshIdx]->getCellZones();
        int cellOffset = offsets[meshIdx].cellOffset;
        std::string meshPrefix = getMeshPrefix(meshIdx);

        for (const auto& zone : zones) {
            CellZone offsetZone;

            // Create unique name for this zone
            offsetZone.name = meshPrefix + "_" + zone.name;
            offsetZone.cellIndices.reserve(zone.cellIndices.size());

            // Apply cell offset to all cell indices
            for (int cellIdx : zone.cellIndices) {
                offsetZone.cellIndices.push_back(cellIdx + cellOffset);
            }

            mergedCellZones.push_back(offsetZone);

            // Add to element block naming map
            // The ExodusWriter creates element blocks based on zone names and element types
            // We need to ensure each mesh's zones produce distinct element blocks
            elementBlockNames[offsetZone.name] = offsetZone.name;
        }
    }

    // If any mesh has no cell zones, create a zone named after the case
    for (size_t meshIdx = 0; meshIdx < readers.size(); ++meshIdx) {
        if (readers[meshIdx]->getNumCellZones() == 0) {
            std::string meshPrefix = getMeshPrefix(meshIdx);
            CellZone defaultZone;
            defaultZone.name = meshPrefix;

            // Add all cells from this mesh to the default zone
            int cellOffset = offsets[meshIdx].cellOffset;
            int numCells = readers[meshIdx]->getNumCells();

            defaultZone.cellIndices.reserve(numCells);
            for (int i = 0; i < numCells; ++i) {
                defaultZone.cellIndices.push_back(cellOffset + i);
            }

            mergedCellZones.push_back(defaultZone);
            elementBlockNames[defaultZone.name] = defaultZone.name;
        }
    }
}

void MergedMeshReader::mergeOwnerNeighbour() {
    // Merge owner and neighbour with cell index offsetting
    for (size_t meshIdx = 0; meshIdx < readers.size(); ++meshIdx) {
        const auto& owner = readers[meshIdx]->getOwner();
        const auto& neighbour = readers[meshIdx]->getNeighbour();
        int cellOffset = offsets[meshIdx].cellOffset;

        // Apply cell offset to owner
        for (int ownerIdx : owner) {
            mergedOwner.push_back(ownerIdx + cellOffset);
        }

        // Apply cell offset to neighbour (handle -1 for boundary faces)
        for (int neighbourIdx : neighbour) {
            if (neighbourIdx >= 0) {
                mergedNeighbour.push_back(neighbourIdx + cellOffset);
            } else {
                mergedNeighbour.push_back(neighbourIdx);  // Keep -1 as is
            }
        }
    }
}

std::string MergedMeshReader::getMeshPrefix(int meshIndex) const {
    // Extract a meaningful name from the case path
    std::string meshName = extractMeshName(casePaths[meshIndex]);

    // If we couldn't extract a name or it's generic, use mesh index
    if (meshName.empty() || meshName == "polyMesh" || meshName == "constant") {
        return "mesh" + std::to_string(meshIndex);
    }

    return meshName;
}

std::string MergedMeshReader::extractMeshName(const std::string& path) const {
    // Try to extract a meaningful name from the path
    // Remove trailing slashes
    std::string cleanPath = path;
    while (!cleanPath.empty() && (cleanPath.back() == '/' || cleanPath.back() == '\\')) {
        cleanPath.pop_back();
    }

    // Try using filesystem to get the directory name
    try {
        std::filesystem::path fsPath(cleanPath);
        std::string dirName = fsPath.filename().string();

        // If it's empty or generic, try parent directory
        if (dirName.empty() || dirName == "." || dirName == "..") {
            dirName = fsPath.parent_path().filename().string();
        }

        // Replace any non-alphanumeric characters with underscores
        for (char& c : dirName) {
            if (!std::isalnum(c)) {
                c = '_';
            }
        }

        return dirName;
    } catch (...) {
        // If filesystem operations fail, fall back to simple parsing
        size_t lastSlash = cleanPath.find_last_of("/\\");
        if (lastSlash != std::string::npos && lastSlash < cleanPath.length() - 1) {
            return cleanPath.substr(lastSlash + 1);
        }
        return cleanPath;
    }
}
