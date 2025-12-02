#ifndef EXODUS_WRITER_H
#define EXODUS_WRITER_H

#include <string>
#include <map>
#include <netcdf.h>
#include "OpenFOAMMeshReader.h"

class ExodusWriter {
public:
    ExodusWriter(const std::string& filename);
    ~ExodusWriter();

    void writeMesh(const OpenFOAMMeshReader& reader);
    void writeMesh(const OpenFOAMMeshReader& reader,
                   const std::map<std::string, std::string>& elementBlockNames,
                   const std::map<std::string, std::string>& sidesetNames);

private:
    std::string filename;
    int ncid;
    std::map<std::string, std::string> customElementBlockNames;
    std::map<std::string, std::string> customSidesetNames;

    void initializeExodusFile(int numNodes, int numElems, int numElemBlocks,
                             int numNodeSets, int numSideSets);
    void writeNodes(const std::vector<Point>& points);
    void writeElements(const OpenFOAMMeshReader& reader);
    void writeSideSets(const OpenFOAMMeshReader& reader);
    void checkError(int status, const std::string& message);
    std::string getBlockName(const std::string& originalName);
    std::string getSidesetName(const std::string& originalName);
    std::vector<int> orderHexNodes(const Cell& cell, const std::vector<Face>& faces, 
                                    const std::vector<Point>& points);
    std::vector<int> orderTetNodes(const Cell& cell, const std::vector<Face>& faces,
                                    const std::vector<Point>& points);
    std::vector<int> orderPyramidNodes(const Cell& cell, const std::vector<Face>& faces,
                                        const std::vector<Point>& points);
    int getHexFaceId(const std::vector<int>& faceNodes, const std::vector<int>& hexNodes);
};

#endif
