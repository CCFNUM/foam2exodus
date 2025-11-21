#include "ExodusWriter.h"
#include <iostream>
#include <stdexcept>
#include <cstring>
#include <algorithm>
#include <set>
#include <cmath>

ExodusWriter::ExodusWriter(const std::string& filename) : filename(filename), ncid(-1) {}

ExodusWriter::~ExodusWriter() {
    if (ncid >= 0) {
        nc_close(ncid);
    }
}

void ExodusWriter::checkError(int status, const std::string& message) {
    if (status != NC_NOERR) {
        throw std::runtime_error(message + ": " + nc_strerror(status));
    }
}

void ExodusWriter::initializeExodusFile(int numNodes, int numElems, int numElemBlocks,
                                       int numNodeSets, int numSideSets) {
    int status = nc_create(filename.c_str(), NC_CLOBBER | NC_64BIT_OFFSET, &ncid);
    checkError(status, "Failed to create Exodus file");
    
    int dim_len_string, dim_len_line, dim_len_name, dim_four;
    int dim_num_nodes, dim_num_elem, dim_num_el_blk;
    int dim_num_side_sets, dim_time_step;
    
    status = nc_def_dim(ncid, "len_string", 33, &dim_len_string);
    checkError(status, "Failed to define len_string dimension");
    status = nc_def_dim(ncid, "len_line", 81, &dim_len_line);
    checkError(status, "Failed to define len_line dimension");
    status = nc_def_dim(ncid, "len_name", 33, &dim_len_name);
    checkError(status, "Failed to define len_name dimension");
    status = nc_def_dim(ncid, "four", 4, &dim_four);
    checkError(status, "Failed to define four dimension");
    status = nc_def_dim(ncid, "time_step", NC_UNLIMITED, &dim_time_step);
    checkError(status, "Failed to define time_step dimension");
    status = nc_def_dim(ncid, "num_dim", 3, NULL);
    checkError(status, "Failed to define num_dim dimension");
    status = nc_def_dim(ncid, "num_nodes", numNodes, &dim_num_nodes);
    checkError(status, "Failed to define num_nodes dimension");
    status = nc_def_dim(ncid, "num_elem", numElems, &dim_num_elem);
    checkError(status, "Failed to define num_elem dimension");
    status = nc_def_dim(ncid, "num_el_blk", numElemBlocks, &dim_num_el_blk);
    checkError(status, "Failed to define num_el_blk dimension");
    status = nc_def_dim(ncid, "num_side_sets", numSideSets, &dim_num_side_sets);
    checkError(status, "Failed to define num_side_sets dimension");
    
    int var_title, var_coord_x, var_coord_y, var_coord_z, var_coor_names;
    int var_eb_status, var_eb_prop1, var_eb_names, var_ss_status, var_ss_prop1;
    int dim3;
    
    status = nc_def_dim(ncid, "3", 3, &dim3);
    checkError(status, "Failed to define 3 dimension");
    
    int dims_coor[2] = {dim3, dim_len_string};
    status = nc_def_var(ncid, "coor_names", NC_CHAR, 2, dims_coor, &var_coor_names);
    checkError(status, "Failed to define coor_names variable");
    status = nc_def_var(ncid, "title", NC_CHAR, 1, &dim_len_line, &var_title);
    checkError(status, "Failed to define title variable");
    status = nc_def_var(ncid, "coordx", NC_DOUBLE, 1, &dim_num_nodes, &var_coord_x);
    checkError(status, "Failed to define coordx variable");
    status = nc_def_var(ncid, "coordy", NC_DOUBLE, 1, &dim_num_nodes, &var_coord_y);
    checkError(status, "Failed to define coordy variable");
    status = nc_def_var(ncid, "coordz", NC_DOUBLE, 1, &dim_num_nodes, &var_coord_z);
    checkError(status, "Failed to define coordz variable");
    status = nc_def_var(ncid, "eb_status", NC_INT, 1, &dim_num_el_blk, &var_eb_status);
    checkError(status, "Failed to define eb_status variable");
    status = nc_def_var(ncid, "eb_prop1", NC_INT, 1, &dim_num_el_blk, &var_eb_prop1);
    checkError(status, "Failed to define eb_prop1 variable");
    int dims_eb_names[2] = {dim_num_el_blk, dim_len_name};
    status = nc_def_var(ncid, "eb_names", NC_CHAR, 2, dims_eb_names, &var_eb_names);
    checkError(status, "Failed to define eb_names variable");
    status = nc_def_var(ncid, "ss_status", NC_INT, 1, &dim_num_side_sets, &var_ss_status);
    checkError(status, "Failed to define ss_status variable");
    status = nc_def_var(ncid, "ss_prop1", NC_INT, 1, &dim_num_side_sets, &var_ss_prop1);
    checkError(status, "Failed to define ss_prop1 variable");
    
    // Write required Exodus II global attributes
    float version = 8.11f;
    status = nc_put_att_float(ncid, NC_GLOBAL, "api_version", NC_FLOAT, 1, &version);
    checkError(status, "Failed to write api_version attribute");
    status = nc_put_att_float(ncid, NC_GLOBAL, "version", NC_FLOAT, 1, &version);
    checkError(status, "Failed to write version attribute");
    
    int word_size = 8;
    status = nc_put_att_int(ncid, NC_GLOBAL, "floating_point_word_size", NC_INT, 1, &word_size);
    checkError(status, "Failed to write floating_point_word_size attribute");
    
    int file_size = 1;
    status = nc_put_att_int(ncid, NC_GLOBAL, "file_size", NC_INT, 1, &file_size);
    checkError(status, "Failed to write file_size attribute");
    
    status = nc_put_att_text(ncid, NC_GLOBAL, "title", 28, "OpenFOAM to Exodus II mesh");
    checkError(status, "Failed to write title attribute");
    
    status = nc_put_att_text(ncid, var_eb_prop1, "name", 2, "ID");
    checkError(status, "Failed to write eb_prop1 name attribute");
    status = nc_put_att_text(ncid, var_ss_prop1, "name", 2, "ID");
    checkError(status, "Failed to write ss_prop1 name attribute");
    
    status = nc_enddef(ncid);
    checkError(status, "Failed to end define mode");
    
    const char* coord_names = "x                                y                                z                                ";
    status = nc_put_var_text(ncid, var_coor_names, coord_names);
    checkError(status, "Failed to write coordinate names");
    
    const char* title_str = "OpenFOAM mesh converted to Exodus II";
    status = nc_put_var_text(ncid, var_title, title_str);
    checkError(status, "Failed to write title");
}

void ExodusWriter::writeNodes(const std::vector<Point>& points) {
    std::vector<double> x_coords, y_coords, z_coords;
    x_coords.reserve(points.size());
    y_coords.reserve(points.size());
    z_coords.reserve(points.size());
    
    for (const auto& p : points) {
        x_coords.push_back(p.x);
        y_coords.push_back(p.y);
        z_coords.push_back(p.z);
    }
    
    int var_id;
    int status;
    
    status = nc_inq_varid(ncid, "coordx", &var_id);
    checkError(status, "Failed to get coordx variable");
    status = nc_put_var_double(ncid, var_id, x_coords.data());
    checkError(status, "Failed to write coordx data");
    
    status = nc_inq_varid(ncid, "coordy", &var_id);
    checkError(status, "Failed to get coordy variable");
    status = nc_put_var_double(ncid, var_id, y_coords.data());
    checkError(status, "Failed to write coordy data");
    
    status = nc_inq_varid(ncid, "coordz", &var_id);
    checkError(status, "Failed to get coordz variable");
    status = nc_put_var_double(ncid, var_id, z_coords.data());
    checkError(status, "Failed to write coordz data");
    
    std::cout << "Wrote " << points.size() << " nodes to Exodus file" << std::endl;
}

std::vector<int> ExodusWriter::orderTetNodes(const Cell& cell, const std::vector<Face>& faces,
                                              const std::vector<Point>& points) {
    if (cell.faceIndices.size() < 4) {
        return std::vector<int>(4, 0);
    }
    
    std::vector<std::vector<int>> triFaces;
    for (int faceIdx : cell.faceIndices) {
        if (faceIdx >= 0 && faceIdx < faces.size()) {
            const auto& face = faces[faceIdx];
            if (face.pointIndices.size() == 3) {
                triFaces.push_back(face.pointIndices);
            }
        }
    }
    
    if (triFaces.size() < 4) {
        std::set<int> nodeSet;
        for (const auto& face : triFaces) {
            for (int n : face) nodeSet.insert(n);
        }
        return std::vector<int>(nodeSet.begin(), nodeSet.end());
    }
    
    std::vector<int> face0 = triFaces[0];
    std::set<int> nodes0(face0.begin(), face0.end());
    
    int apexNode = -1;
    for (const auto& face : triFaces) {
        for (int node : face) {
            if (nodes0.find(node) == nodes0.end()) {
                apexNode = node;
                break;
            }
        }
        if (apexNode != -1) break;
    }
    
    if (apexNode == -1) {
        apexNode = face0[0];
    }
    
    std::vector<int> orderedNodes(4);
    orderedNodes[0] = face0[0];
    orderedNodes[1] = face0[1];
    orderedNodes[2] = face0[2];
    orderedNodes[3] = apexNode;
    
    return orderedNodes;
}

std::vector<int> ExodusWriter::orderPyramidNodes(const Cell& cell, const std::vector<Face>& faces,
                                                  const std::vector<Point>& points) {
    if (cell.faceIndices.size() < 5) {
        return std::vector<int>(5, 0);
    }
    
    std::vector<int> quadFace;
    std::vector<std::vector<int>> triFaces;
    
    for (int faceIdx : cell.faceIndices) {
        if (faceIdx >= 0 && faceIdx < faces.size()) {
            const auto& face = faces[faceIdx];
            if (face.pointIndices.size() == 4) {
                quadFace = face.pointIndices;
            } else if (face.pointIndices.size() == 3) {
                triFaces.push_back(face.pointIndices);
            }
        }
    }
    
    if (quadFace.empty() || triFaces.empty()) {
        std::set<int> nodeSet;
        for (int faceIdx : cell.faceIndices) {
            if (faceIdx >= 0 && faceIdx < faces.size()) {
                for (int n : faces[faceIdx].pointIndices) {
                    nodeSet.insert(n);
                }
            }
        }
        return std::vector<int>(nodeSet.begin(), nodeSet.end());
    }
    
    std::set<int> baseSet(quadFace.begin(), quadFace.end());
    int apexNode = -1;
    for (const auto& triFace : triFaces) {
        for (int node : triFace) {
            if (baseSet.find(node) == baseSet.end()) {
                apexNode = node;
                break;
            }
        }
        if (apexNode != -1) break;
    }
    
    if (apexNode == -1) {
        apexNode = quadFace[0];
    }
    
    std::vector<int> orderedNodes(5);
    orderedNodes[0] = quadFace[0];
    orderedNodes[1] = quadFace[1];
    orderedNodes[2] = quadFace[2];
    orderedNodes[3] = quadFace[3];
    orderedNodes[4] = apexNode;
    
    return orderedNodes;
}

std::vector<int> ExodusWriter::orderHexNodes(const Cell& cell, const std::vector<Face>& faces, 
                                              const std::vector<Point>& points) {
    if (cell.faceIndices.size() < 6) {
        return std::vector<int>(8, 0);
    }
    
    std::vector<std::vector<int>> quadFaces;
    for (int faceIdx : cell.faceIndices) {
        if (faceIdx >= 0 && faceIdx < faces.size()) {
            const auto& face = faces[faceIdx];
            if (face.pointIndices.size() == 4) {
                quadFaces.push_back(face.pointIndices);
            }
        }
    }
    
    if (quadFaces.size() != 6) {
        std::set<int> nodeSet;
        for (const auto& face : quadFaces) {
            for (int n : face) nodeSet.insert(n);
        }
        return std::vector<int>(nodeSet.begin(), nodeSet.end());
    }
    
    std::vector<int> bottomFace = quadFaces[0];
    std::set<int> bottomSet(bottomFace.begin(), bottomFace.end());
    
    std::vector<int> topFace;
    int topFaceIdx = -1;
    for (size_t i = 1; i < quadFaces.size(); ++i) {
        std::set<int> faceSet(quadFaces[i].begin(), quadFaces[i].end());
        std::vector<int> intersection;
        std::set_intersection(bottomSet.begin(), bottomSet.end(),
                            faceSet.begin(), faceSet.end(),
                            std::back_inserter(intersection));
        
        if (intersection.empty()) {
            topFace = quadFaces[i];
            topFaceIdx = i;
            break;
        }
    }
    
    if (topFace.empty()) {
        std::set<int> nodeSet;
        for (const auto& face : quadFaces) {
            for (int n : face) nodeSet.insert(n);
        }
        return std::vector<int>(nodeSet.begin(), nodeSet.end());
    }
    
    std::vector<std::vector<int>> sideFaces;
    for (size_t i = 1; i < quadFaces.size(); ++i) {
        if (i != topFaceIdx) {
            sideFaces.push_back(quadFaces[i]);
        }
    }
    
    std::vector<int> orderedNodes(8);
    for (int i = 0; i < 4; i++) {
        orderedNodes[i] = bottomFace[i];
    }
    
    for (int i = 0; i < 4; i++) {
        int n0 = bottomFace[i];
        int n1 = bottomFace[(i + 1) % 4];
        
        int topNode = -1;
        for (const auto& sideFace : sideFaces) {
            int n0Idx = -1, n1Idx = -1;
            
            for (size_t j = 0; j < sideFace.size(); ++j) {
                if (sideFace[j] == n0) n0Idx = j;
                if (sideFace[j] == n1) n1Idx = j;
            }
            
            if (n0Idx >= 0 && n1Idx >= 0) {
                int nextIdx = (n0Idx + 1) % 4;
                int prevIdx = (n0Idx + 3) % 4;
                
                int candidate = -1;
                if (nextIdx == n1Idx) {
                    candidate = sideFace[prevIdx];
                } else if (prevIdx == n1Idx) {
                    candidate = sideFace[nextIdx];
                }
                
                if (candidate != -1) {
                    auto it = std::find(topFace.begin(), topFace.end(), candidate);
                    if (it != topFace.end()) {
                        topNode = candidate;
                        break;
                    }
                }
            }
        }
        
        if (topNode == -1) {
            double minDist = 1e10;
            for (int node : topFace) {
                bool alreadyUsed = false;
                for (int j = 0; j < i; j++) {
                    if (orderedNodes[4 + j] == node) {
                        alreadyUsed = true;
                        break;
                    }
                }
                if (!alreadyUsed) {
                    const Point& pb = points[n0];
                    const Point& pt = points[node];
                    double dx = pt.x - pb.x;
                    double dy = pt.y - pb.y;
                    double dz = pt.z - pb.z;
                    double dist = dx*dx + dy*dy + dz*dz;
                    
                    if (dist < minDist) {
                        minDist = dist;
                        topNode = node;
                    }
                }
            }
        }
        
        if (topNode != -1) {
            orderedNodes[4 + i] = topNode;
        } else {
            orderedNodes[4 + i] = topFace[i];
        }
    }
    
    auto computeVolume = [&](const std::vector<int>& nodes) -> double {
        const Point& p0 = points[nodes[0]];
        const Point& p1 = points[nodes[1]];
        const Point& p3 = points[nodes[3]];
        const Point& p4 = points[nodes[4]];
        
        double v1x = p1.x - p0.x, v1y = p1.y - p0.y, v1z = p1.z - p0.z;
        double v2x = p3.x - p0.x, v2y = p3.y - p0.y, v2z = p3.z - p0.z;
        double v3x = p4.x - p0.x, v3y = p4.y - p0.y, v3z = p4.z - p0.z;
        
        return v1x * (v2y * v3z - v2z * v3y) +
               v1y * (v2z * v3x - v2x * v3z) +
               v1z * (v2x * v3y - v2y * v3x);
    };
    
    double volume = computeVolume(orderedNodes);
    if (volume < 0) {
        std::swap(orderedNodes[1], orderedNodes[3]);
        std::swap(orderedNodes[5], orderedNodes[7]);
    }
    
    return orderedNodes;
}

void ExodusWriter::writeElements(const OpenFOAMMeshReader& reader) {
    const auto& cells = reader.getCells();
    const auto& faces = reader.getFaces();
    const auto& points = reader.getPoints();
    const auto& cellZones = reader.getCellZones();
    
    nc_redef(ncid);
    
    struct BlockInfo {
        std::string name;
        std::vector<int> cellIndices;
        std::string cellType;
    };
    
    std::vector<BlockInfo> blocks;
    
    if (!cellZones.empty() && cellZones.size() > 1) {
        std::set<int> zonedCells;
        
        for (const auto& zone : cellZones) {
            std::map<std::string, std::vector<int>> zoneElemsByType;
            
            for (int cellIdx : zone.cellIndices) {
                if (cellIdx < cells.size()) {
                    zoneElemsByType[cells[cellIdx].type].push_back(cellIdx);
                    zonedCells.insert(cellIdx);
                }
            }
            
            for (const auto& [cellType, cellIndices] : zoneElemsByType) {
                BlockInfo block;
                block.name = zone.name + "-" + cellType;
                block.cellIndices = cellIndices;
                block.cellType = cellType;
                blocks.push_back(block);
            }
        }
        
        std::map<std::string, std::vector<int>> unzonedElemsByType;
        for (size_t i = 0; i < cells.size(); ++i) {
            if (zonedCells.find(i) == zonedCells.end()) {
                unzonedElemsByType[cells[i].type].push_back(i);
            }
        }
        
        for (const auto& [cellType, cellIndices] : unzonedElemsByType) {
            BlockInfo block;
            block.name = "unzoned-" + cellType;
            block.cellIndices = cellIndices;
            block.cellType = cellType;
            blocks.push_back(block);
        }
    } else {
        std::map<std::string, std::vector<int>> elemsByType;
        for (size_t i = 0; i < cells.size(); ++i) {
            elemsByType[cells[i].type].push_back(i);
        }
        
        for (const auto& [cellType, cellIndices] : elemsByType) {
            BlockInfo block;
            block.name = "fluid-" + cellType;
            block.cellIndices = cellIndices;
            block.cellType = cellType;
            blocks.push_back(block);
        }
    }
    
    int blockId = 1;
    for (const auto& block : blocks) {
        int numElemsInBlock = block.cellIndices.size();
        
        int numNodesPerElem = 8;
        std::string exoType = "HEX8";
        
        if (block.cellType == "tet") {
            numNodesPerElem = 4;
            exoType = "TETRA4";
        } else if (block.cellType == "pyr") {
            numNodesPerElem = 5;
            exoType = "PYRAMID5";
        } else if (block.cellType == "hex") {
            numNodesPerElem = 8;
            exoType = "HEX8";
        }
        
        int dim_num_el_in_blk, dim_num_nod_per_el;
        std::string blk_num = std::to_string(blockId);
        
        int status;
        status = nc_def_dim(ncid, ("num_el_in_blk" + blk_num).c_str(), numElemsInBlock, &dim_num_el_in_blk);
        checkError(status, "Failed to define element block dimension");
        
        status = nc_def_dim(ncid, ("num_nod_per_el" + blk_num).c_str(), numNodesPerElem, &dim_num_nod_per_el);
        checkError(status, "Failed to define nodes per element dimension");
        
        int var_connect;
        int dims_connect[2] = {dim_num_el_in_blk, dim_num_nod_per_el};
        status = nc_def_var(ncid, ("connect" + blk_num).c_str(), NC_INT, 2, dims_connect, &var_connect);
        checkError(status, "Failed to define connectivity variable");
        
        status = nc_put_att_text(ncid, var_connect, "elem_type", exoType.length(), exoType.c_str());
        checkError(status, "Failed to set element type attribute");
        
        blockId++;
    }
    
    nc_enddef(ncid);
    
    blockId = 1;
    for (const auto& block : blocks) {
        std::vector<int> connectivity;
        
        for (int cellIdx : block.cellIndices) {
            const auto& cell = cells[cellIdx];
            
            std::vector<int> nodes;
            if (block.cellType == "hex") {
                nodes = orderHexNodes(cell, faces, points);
            } else if (block.cellType == "tet") {
                nodes = orderTetNodes(cell, faces, points);
            } else if (block.cellType == "pyr") {
                nodes = orderPyramidNodes(cell, faces, points);
            } else {
                std::set<int> nodeSet;
                for (int faceIdx : cell.faceIndices) {
                    if (faceIdx >= 0 && faceIdx < faces.size()) {
                        for (int nodeIdx : faces[faceIdx].pointIndices) {
                            nodeSet.insert(nodeIdx);
                        }
                    }
                }
                nodes.assign(nodeSet.begin(), nodeSet.end());
            }
            
            int numNodesPerElem = (block.cellType == "tet") ? 4 : (block.cellType == "pyr") ? 5 : 8;
            for (int i = 0; i < numNodesPerElem; ++i) {
                connectivity.push_back(i < nodes.size() ? nodes[i] + 1 : 1);
            }
        }
        
        int var_id;
        nc_inq_varid(ncid, ("connect" + std::to_string(blockId)).c_str(), &var_id);
        nc_put_var_int(ncid, var_id, connectivity.data());
        
        int eb_status = 1;
        int var_status;
        nc_inq_varid(ncid, "eb_status", &var_status);
        size_t index = blockId - 1;
        nc_put_var1_int(ncid, var_status, &index, &eb_status);
        
        int var_prop;
        nc_inq_varid(ncid, "eb_prop1", &var_prop);
        nc_put_var1_int(ncid, var_prop, &index, &blockId);
        
        int var_eb_names;
        nc_inq_varid(ncid, "eb_names", &var_eb_names);
        char name_buffer[33];
        memset(name_buffer, ' ', 33);
        size_t copy_len = std::min(block.name.length(), (size_t)32);
        memcpy(name_buffer, block.name.c_str(), copy_len);
        name_buffer[32] = '\0';
        size_t start[2] = {index, 0};
        size_t count[2] = {1, 33};
        nc_put_vara_text(ncid, var_eb_names, start, count, name_buffer);
        
        std::cout << "Wrote element block " << blockId << " (" << block.name << ") with " 
                  << block.cellIndices.size() << " elements" << std::endl;
        
        blockId++;
    }
}

int ExodusWriter::getHexFaceId(const std::vector<int>& faceNodes, const std::vector<int>& hexNodes) {
    std::set<int> faceSet(faceNodes.begin(), faceNodes.end());
    
    std::vector<std::vector<int>> hexFaces = {
        {hexNodes[0], hexNodes[1], hexNodes[5], hexNodes[4]},
        {hexNodes[1], hexNodes[2], hexNodes[6], hexNodes[5]},
        {hexNodes[2], hexNodes[3], hexNodes[7], hexNodes[6]},
        {hexNodes[3], hexNodes[0], hexNodes[4], hexNodes[7]},
        {hexNodes[0], hexNodes[3], hexNodes[2], hexNodes[1]},
        {hexNodes[4], hexNodes[5], hexNodes[6], hexNodes[7]}
    };
    
    for (int i = 0; i < 6; ++i) {
        std::set<int> hexFaceSet(hexFaces[i].begin(), hexFaces[i].end());
        if (faceSet == hexFaceSet) {
            return i + 1;
        }
    }
    
    return 1;
}

void ExodusWriter::writeSideSets(const OpenFOAMMeshReader& reader) {
    const auto& patches = reader.getBoundaryPatches();
    const auto& faces = reader.getFaces();
    const auto& owner = reader.getOwner();
    const auto& cells = reader.getCells();
    const auto& points = reader.getPoints();
    
    if (patches.empty()) {
        std::cout << "No boundary patches to write as sidesets" << std::endl;
        return;
    }
    
    std::map<int, std::vector<int>> cellToOrderedNodes;
    for (size_t i = 0; i < cells.size(); ++i) {
        const auto& cell = cells[i];
        std::vector<int> orderedNodes;
        
        if (cell.type == "hex") {
            orderedNodes = orderHexNodes(cell, faces, points);
        } else if (cell.type == "tet") {
            orderedNodes = orderTetNodes(cell, faces, points);
        } else if (cell.type == "pyr") {
            orderedNodes = orderPyramidNodes(cell, faces, points);
        } else {
            std::set<int> nodeSet;
            for (int faceIdx : cell.faceIndices) {
                if (faceIdx >= 0 && faceIdx < faces.size()) {
                    for (int nodeIdx : faces[faceIdx].pointIndices) {
                        nodeSet.insert(nodeIdx);
                    }
                }
            }
            orderedNodes.assign(nodeSet.begin(), nodeSet.end());
        }
        cellToOrderedNodes[i] = orderedNodes;
    }
    
    nc_redef(ncid);
    
    for (size_t i = 0; i < patches.size(); ++i) {
        const auto& patch = patches[i];
        std::string ss_num = std::to_string(i + 1);
        
        int dim_num_side_ss;
        nc_def_dim(ncid, ("num_side_ss" + ss_num).c_str(), patch.nFaces, &dim_num_side_ss);
        
        int var_elem_ss, var_side_ss;
        nc_def_var(ncid, ("elem_ss" + ss_num).c_str(), NC_INT, 1, &dim_num_side_ss, &var_elem_ss);
        nc_def_var(ncid, ("side_ss" + ss_num).c_str(), NC_INT, 1, &dim_num_side_ss, &var_side_ss);
        
        nc_put_att_text(ncid, var_elem_ss, "name", patch.name.length(), patch.name.c_str());
    }
    
    int var_ss_names;
    int dim_num_side_sets, dim_len_name;
    int status;
    status = nc_inq_dimid(ncid, "num_side_sets", &dim_num_side_sets);
    checkError(status, "Failed to get num_side_sets dimension");
    
    status = nc_inq_dimid(ncid, "len_name", &dim_len_name);
    checkError(status, "Failed to get len_name dimension");
    
    int dims_ss_names[2] = {dim_num_side_sets, dim_len_name};
    status = nc_def_var(ncid, "ss_names", NC_CHAR, 2, dims_ss_names, &var_ss_names);
    checkError(status, "Failed to define ss_names variable");
    
    nc_enddef(ncid);
    
    for (size_t i = 0; i < patches.size(); ++i) {
        const auto& patch = patches[i];
        std::string ss_num = std::to_string(i + 1);
        
        std::vector<int> elem_list, side_list;
        elem_list.reserve(patch.nFaces);
        side_list.reserve(patch.nFaces);
        
        for (int j = 0; j < patch.nFaces; ++j) {
            int faceIdx = patch.startFace + j;
            if (faceIdx < owner.size()) {
                int cellIdx = owner[faceIdx];
                elem_list.push_back(cellIdx + 1);
                
                const auto& face = faces[faceIdx];
                const auto& cell = cells[cellIdx];
                int sideId = 1;
                
                if (cell.type == "hex" && cellToOrderedNodes.count(cellIdx)) {
                    sideId = getHexFaceId(face.pointIndices, cellToOrderedNodes[cellIdx]);
                }
                
                side_list.push_back(sideId);
            }
        }
        
        int var_id;
        nc_inq_varid(ncid, ("elem_ss" + ss_num).c_str(), &var_id);
        nc_put_var_int(ncid, var_id, elem_list.data());
        
        nc_inq_varid(ncid, ("side_ss" + ss_num).c_str(), &var_id);
        nc_put_var_int(ncid, var_id, side_list.data());
        
        int ss_status = 1;
        int var_status;
        nc_inq_varid(ncid, "ss_status", &var_status);
        size_t index = i;
        nc_put_var1_int(ncid, var_status, &index, &ss_status);
        
        int ss_id = i + 1;
        int var_prop;
        nc_inq_varid(ncid, "ss_prop1", &var_prop);
        nc_put_var1_int(ncid, var_prop, &index, &ss_id);
        
        char name_padded[33];
        std::memset(name_padded, 0, 33);
        std::strncpy(name_padded, patch.name.c_str(), 32);
        
        nc_inq_varid(ncid, "ss_names", &var_id);
        size_t start[2] = {i, 0};
        size_t count[2] = {1, 33};
        nc_put_vara_text(ncid, var_id, start, count, name_padded);
        
        std::cout << "Wrote sideset " << (i + 1) << ": " << patch.name 
                  << " with " << patch.nFaces << " faces" << std::endl;
    }
}

void ExodusWriter::writeMesh(const OpenFOAMMeshReader& reader) {
    int numNodes = reader.getNumPoints();
    int numElems = reader.getNumCells();
    int numNodeSets = 0;
    int numSideSets = reader.getNumBoundaryPatches();
    
    const auto& cells = reader.getCells();
    const auto& cellZones = reader.getCellZones();
    
    int numElemBlocks = 0;
    if (!cellZones.empty()) {
        std::set<int> zonedCells;
        
        for (const auto& zone : cellZones) {
            std::set<std::string> typesInZone;
            for (int cellIdx : zone.cellIndices) {
                if (cellIdx < cells.size()) {
                    typesInZone.insert(cells[cellIdx].type);
                    zonedCells.insert(cellIdx);
                }
            }
            numElemBlocks += typesInZone.size();
        }
        
        std::set<std::string> unzonedTypes;
        for (size_t i = 0; i < cells.size(); ++i) {
            if (zonedCells.find(i) == zonedCells.end()) {
                unzonedTypes.insert(cells[i].type);
            }
        }
        numElemBlocks += unzonedTypes.size();
    } else {
        std::map<std::string, std::vector<int>> elemsByType;
        for (size_t i = 0; i < cells.size(); ++i) {
            elemsByType[cells[i].type].push_back(i);
        }
        numElemBlocks = elemsByType.size();
    }
    
    initializeExodusFile(numNodes, numElems, numElemBlocks, numNodeSets, numSideSets);
    writeNodes(reader.getPoints());
    writeElements(reader);
    writeSideSets(reader);
    
    std::cout << "Successfully wrote Exodus II file: " << filename << std::endl;
}
