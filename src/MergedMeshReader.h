#ifndef MERGED_MESH_READER_H
#define MERGED_MESH_READER_H

#include "OpenFOAMMeshReader.h"
#include <string>
#include <vector>
#include <map>

/**
 * MergedMeshReader - Reads multiple OpenFOAM meshes and merges them into a single mesh
 *
 * This class maintains separate element blocks for each source mesh by:
 * - Offsetting all indices (nodes, cells, faces) appropriately
 * - Creating unique names for zones and boundary patches from each mesh
 * - Providing a unified interface compatible with OpenFOAMMeshReader
 */
class MergedMeshReader {
public:
    MergedMeshReader(const std::vector<std::string>& casePaths);

    void readMeshes();

    // Interface compatible with OpenFOAMMeshReader
    int getNumPoints() const { return mergedPoints.size(); }
    int getNumCells() const { return mergedCells.size(); }
    int getNumBoundaryPatches() const { return mergedBoundaryPatches.size(); }
    int getNumCellZones() const { return mergedCellZones.size(); }

    const std::vector<Point>& getPoints() const { return mergedPoints; }
    const std::vector<Face>& getFaces() const { return mergedFaces; }
    const std::vector<Cell>& getCells() const { return mergedCells; }
    const std::vector<BoundaryPatch>& getBoundaryPatches() const { return mergedBoundaryPatches; }
    const std::vector<CellZone>& getCellZones() const { return mergedCellZones; }
    const std::vector<int>& getOwner() const { return mergedOwner; }
    const std::vector<int>& getNeighbour() const { return mergedNeighbour; }

    // Custom naming maps for distinct element blocks and sidesets
    const std::map<std::string, std::string>& getElementBlockNames() const { return elementBlockNames; }
    const std::map<std::string, std::string>& getSidesetNames() const { return sidesetNames; }

private:
    std::vector<std::string> casePaths;
    std::vector<OpenFOAMMeshReader*> readers;

    // Merged data
    std::vector<Point> mergedPoints;
    std::vector<Face> mergedFaces;
    std::vector<Cell> mergedCells;
    std::vector<BoundaryPatch> mergedBoundaryPatches;
    std::vector<CellZone> mergedCellZones;
    std::vector<int> mergedOwner;
    std::vector<int> mergedNeighbour;

    // Custom naming maps
    std::map<std::string, std::string> elementBlockNames;
    std::map<std::string, std::string> sidesetNames;

    // Offset tracking
    struct MeshOffsets {
        int nodeOffset;
        int faceOffset;
        int cellOffset;
    };
    std::vector<MeshOffsets> offsets;

    void mergeMeshes();
    void mergePoints();
    void mergeFaces();
    void mergeCells();
    void mergeBoundaryPatches();
    void mergeCellZones();
    void mergeOwnerNeighbour();

    std::string getMeshPrefix(int meshIndex) const;
    std::string extractMeshName(const std::string& path) const;
};

#endif
