// File       : ExodusWriter.cpp
// Created    : Thu Mar 19 2026
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2026 CCFNUM, Lucerne University of Applied Sciences and
// Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "ExodusWriter.h"
#include "MergedMeshReader.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <set>
#include <stdexcept>
#include <unordered_map>

ExodusWriter::ExodusWriter(const std::string& filename)
    : filename(filename), ncid(-1)
{
}

ExodusWriter::~ExodusWriter()
{
    if (ncid >= 0)
    {
        nc_close(ncid);
    }
}

void ExodusWriter::checkError(int status, const std::string& message)
{
    if (status != NC_NOERR)
    {
        throw std::runtime_error(message + ": " + nc_strerror(status));
    }
}

void ExodusWriter::initializeExodusFile(int numNodes,
                                        int numElems,
                                        int numElemBlocks,
                                        int numNodeSets,
                                        int numSideSets)
{
    int status =
        nc_create(filename.c_str(), NC_CLOBBER | NC_64BIT_OFFSET, &ncid);
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
    status =
        nc_def_var(ncid, "coor_names", NC_CHAR, 2, dims_coor, &var_coor_names);
    checkError(status, "Failed to define coor_names variable");
    status = nc_def_var(ncid, "title", NC_CHAR, 1, &dim_len_line, &var_title);
    checkError(status, "Failed to define title variable");
    status =
        nc_def_var(ncid, "coordx", NC_DOUBLE, 1, &dim_num_nodes, &var_coord_x);
    checkError(status, "Failed to define coordx variable");
    status =
        nc_def_var(ncid, "coordy", NC_DOUBLE, 1, &dim_num_nodes, &var_coord_y);
    checkError(status, "Failed to define coordy variable");
    status =
        nc_def_var(ncid, "coordz", NC_DOUBLE, 1, &dim_num_nodes, &var_coord_z);
    checkError(status, "Failed to define coordz variable");
    status = nc_def_var(
        ncid, "eb_status", NC_INT, 1, &dim_num_el_blk, &var_eb_status);
    checkError(status, "Failed to define eb_status variable");
    status =
        nc_def_var(ncid, "eb_prop1", NC_INT, 1, &dim_num_el_blk, &var_eb_prop1);
    checkError(status, "Failed to define eb_prop1 variable");
    int dims_eb_names[2] = {dim_num_el_blk, dim_len_name};
    status =
        nc_def_var(ncid, "eb_names", NC_CHAR, 2, dims_eb_names, &var_eb_names);
    checkError(status, "Failed to define eb_names variable");
    status = nc_def_var(
        ncid, "ss_status", NC_INT, 1, &dim_num_side_sets, &var_ss_status);
    checkError(status, "Failed to define ss_status variable");
    status = nc_def_var(
        ncid, "ss_prop1", NC_INT, 1, &dim_num_side_sets, &var_ss_prop1);
    checkError(status, "Failed to define ss_prop1 variable");

    // Write required Exodus II global attributes
    float version = 8.11f;
    status =
        nc_put_att_float(ncid, NC_GLOBAL, "api_version", NC_FLOAT, 1, &version);
    checkError(status, "Failed to write api_version attribute");
    status =
        nc_put_att_float(ncid, NC_GLOBAL, "version", NC_FLOAT, 1, &version);
    checkError(status, "Failed to write version attribute");

    int word_size = 8;
    status = nc_put_att_int(
        ncid, NC_GLOBAL, "floating_point_word_size", NC_INT, 1, &word_size);
    checkError(status, "Failed to write floating_point_word_size attribute");

    int file_size = 1;
    status =
        nc_put_att_int(ncid, NC_GLOBAL, "file_size", NC_INT, 1, &file_size);
    checkError(status, "Failed to write file_size attribute");

    status = nc_put_att_text(
        ncid, NC_GLOBAL, "title", 28, "OpenFOAM to Exodus II mesh");
    checkError(status, "Failed to write title attribute");

    status = nc_put_att_text(ncid, var_eb_prop1, "name", 2, "ID");
    checkError(status, "Failed to write eb_prop1 name attribute");
    status = nc_put_att_text(ncid, var_ss_prop1, "name", 2, "ID");
    checkError(status, "Failed to write ss_prop1 name attribute");

    status = nc_enddef(ncid);
    checkError(status, "Failed to end define mode");

    const char* coord_names =
        "x                                y                                z   "
        "                             ";
    status = nc_put_var_text(ncid, var_coor_names, coord_names);
    checkError(status, "Failed to write coordinate names");

    const char* title_str = "OpenFOAM mesh converted to Exodus II";
    status = nc_put_var_text(ncid, var_title, title_str);
    checkError(status, "Failed to write title");
}

void ExodusWriter::writeNodes(const std::vector<Point>& points)
{
    size_t total = points.size() + polyExtraPoints.size();
    std::vector<double> x_coords, y_coords, z_coords;
    x_coords.reserve(total);
    y_coords.reserve(total);
    z_coords.reserve(total);

    for (const auto& p : points)
    {
        x_coords.push_back(p.x);
        y_coords.push_back(p.y);
        z_coords.push_back(p.z);
    }
    // Centroid nodes added by the polyhedral decomposition (empty otherwise).
    for (const auto& p : polyExtraPoints)
    {
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

    std::cout << "Wrote " << total << " nodes to Exodus file" << std::endl;
}

std::vector<int> ExodusWriter::orderTetNodes(const Cell& cell,
                                             const std::vector<Face>& faces,
                                             const std::vector<Point>& points)
{
    if (cell.faceIndices.size() < 4)
    {
        return std::vector<int>(4, 0);
    }

    std::vector<std::vector<int>> triFaces;
    for (int faceIdx : cell.faceIndices)
    {
        if (faceIdx >= 0 && faceIdx < faces.size())
        {
            const auto& face = faces[faceIdx];
            if (face.pointIndices.size() == 3)
            {
                triFaces.push_back(face.pointIndices);
            }
        }
    }

    if (triFaces.size() < 4)
    {
        std::set<int> nodeSet;
        for (const auto& face : triFaces)
        {
            for (int n : face)
                nodeSet.insert(n);
        }
        return std::vector<int>(nodeSet.begin(), nodeSet.end());
    }

    std::vector<int> face0 = triFaces[0];
    std::set<int> nodes0(face0.begin(), face0.end());

    int apexNode = -1;
    for (const auto& face : triFaces)
    {
        for (int node : face)
        {
            if (nodes0.find(node) == nodes0.end())
            {
                apexNode = node;
                break;
            }
        }
        if (apexNode != -1)
            break;
    }

    if (apexNode == -1)
    {
        apexNode = face0[0];
    }

    std::vector<int> orderedNodes(4);
    orderedNodes[0] = face0[0];
    orderedNodes[1] = face0[1];
    orderedNodes[2] = face0[2];
    orderedNodes[3] = apexNode;

    // Ensure positive orientation: scalar triple product
    // (p1-p0)·((p2-p0)×(p3-p0)) > 0
    const Point& p0 = points[orderedNodes[0]];
    const Point& p1 = points[orderedNodes[1]];
    const Point& p2 = points[orderedNodes[2]];
    const Point& p3 = points[orderedNodes[3]];
    double v1x = p1.x - p0.x, v1y = p1.y - p0.y, v1z = p1.z - p0.z;
    double v2x = p2.x - p0.x, v2y = p2.y - p0.y, v2z = p2.z - p0.z;
    double v3x = p3.x - p0.x, v3y = p3.y - p0.y, v3z = p3.z - p0.z;
    double vol = v1x * (v2y * v3z - v2z * v3y) + v1y * (v2z * v3x - v2x * v3z) +
                 v1z * (v2x * v3y - v2y * v3x);
    if (vol < 0)
    {
        std::swap(orderedNodes[1], orderedNodes[2]);
    }

    return orderedNodes;
}

std::vector<int>
ExodusWriter::orderPyramidNodes(const Cell& cell,
                                const std::vector<Face>& faces,
                                const std::vector<Point>& points)
{
    if (cell.faceIndices.size() < 5)
    {
        return std::vector<int>(5, 0);
    }

    std::vector<int> quadFace;
    std::vector<std::vector<int>> triFaces;

    for (int faceIdx : cell.faceIndices)
    {
        if (faceIdx >= 0 && faceIdx < faces.size())
        {
            const auto& face = faces[faceIdx];
            if (face.pointIndices.size() == 4)
            {
                quadFace = face.pointIndices;
            }
            else if (face.pointIndices.size() == 3)
            {
                triFaces.push_back(face.pointIndices);
            }
        }
    }

    if (quadFace.empty() || triFaces.empty())
    {
        std::set<int> nodeSet;
        for (int faceIdx : cell.faceIndices)
        {
            if (faceIdx >= 0 && faceIdx < faces.size())
            {
                for (int n : faces[faceIdx].pointIndices)
                {
                    nodeSet.insert(n);
                }
            }
        }
        return std::vector<int>(nodeSet.begin(), nodeSet.end());
    }

    std::set<int> baseSet(quadFace.begin(), quadFace.end());
    int apexNode = -1;
    for (const auto& triFace : triFaces)
    {
        for (int node : triFace)
        {
            if (baseSet.find(node) == baseSet.end())
            {
                apexNode = node;
                break;
            }
        }
        if (apexNode != -1)
            break;
    }

    if (apexNode == -1)
    {
        apexNode = quadFace[0];
    }

    std::vector<int> orderedNodes(5);
    orderedNodes[0] = quadFace[0];
    orderedNodes[1] = quadFace[1];
    orderedNodes[2] = quadFace[2];
    orderedNodes[3] = quadFace[3];
    orderedNodes[4] = apexNode;

    // Ensure positive orientation: the cross product of two base-quad edges
    // crossed with the vector to apex should be positive (apex is on the
    // positive-normal side of the base).
    const Point& p0 = points[orderedNodes[0]];
    const Point& p1 = points[orderedNodes[1]];
    const Point& p3 = points[orderedNodes[3]];
    const Point& p4 = points[orderedNodes[4]];
    // Base-quad normal: (p1-p0) x (p3-p0)
    double v1x = p1.x - p0.x, v1y = p1.y - p0.y, v1z = p1.z - p0.z;
    double v2x = p3.x - p0.x, v2y = p3.y - p0.y, v2z = p3.z - p0.z;
    double nx = v1y * v2z - v1z * v2y;
    double ny = v1z * v2x - v1x * v2z;
    double nz = v1x * v2y - v1y * v2x;
    // Vector from p0 to apex
    double ax = p4.x - p0.x, ay = p4.y - p0.y, az = p4.z - p0.z;
    // dot(normal, apex_vec) should be positive
    if (nx * ax + ny * ay + nz * az < 0)
    {
        std::swap(orderedNodes[1], orderedNodes[3]);
    }

    return orderedNodes;
}

std::vector<int> ExodusWriter::orderHexNodes(const Cell& cell,
                                             const std::vector<Face>& faces,
                                             const std::vector<Point>& points)
{
    if (cell.faceIndices.size() < 6)
    {
        return std::vector<int>(8, 0);
    }

    std::vector<std::vector<int>> quadFaces;
    for (int faceIdx : cell.faceIndices)
    {
        if (faceIdx >= 0 && faceIdx < faces.size())
        {
            const auto& face = faces[faceIdx];
            if (face.pointIndices.size() == 4)
            {
                quadFaces.push_back(face.pointIndices);
            }
        }
    }

    if (quadFaces.size() != 6)
    {
        std::set<int> nodeSet;
        for (const auto& face : quadFaces)
        {
            for (int n : face)
                nodeSet.insert(n);
        }
        return std::vector<int>(nodeSet.begin(), nodeSet.end());
    }

    // Pick bottom face as the one with minimum mean Z centroid so that
    // all elements use the same geometric face as "bottom" regardless of
    // the arbitrary face-index order in the OpenFOAM cell's face list.
    int bottomIdx = 0;
    double minZ = 1e30;
    for (int i = 0; i < (int)quadFaces.size(); ++i)
    {
        double z = 0.0;
        for (int n : quadFaces[i])
            z += points[n].z;
        z /= 4.0;
        if (z < minZ)
        {
            minZ = z;
            bottomIdx = i;
        }
    }
    std::vector<int> bottomFace = quadFaces[bottomIdx];
    std::set<int> bottomSet(bottomFace.begin(), bottomFace.end());

    std::vector<int> topFace;
    int topFaceIdx = -1;
    for (size_t i = 0; i < quadFaces.size(); ++i)
    {
        if ((int)i == bottomIdx)
            continue;
        std::set<int> faceSet(quadFaces[i].begin(), quadFaces[i].end());
        std::vector<int> intersection;
        std::set_intersection(bottomSet.begin(),
                              bottomSet.end(),
                              faceSet.begin(),
                              faceSet.end(),
                              std::back_inserter(intersection));

        if (intersection.empty())
        {
            topFace = quadFaces[i];
            topFaceIdx = i;
            break;
        }
    }

    if (topFace.empty())
    {
        std::set<int> nodeSet;
        for (const auto& face : quadFaces)
        {
            for (int n : face)
                nodeSet.insert(n);
        }
        return std::vector<int>(nodeSet.begin(), nodeSet.end());
    }

    std::vector<std::vector<int>> sideFaces;
    for (size_t i = 0; i < quadFaces.size(); ++i)
    {
        if ((int)i != bottomIdx && (int)i != topFaceIdx)
        {
            sideFaces.push_back(quadFaces[i]);
        }
    }

    std::vector<int> orderedNodes(8);
    for (int i = 0; i < 4; i++)
    {
        orderedNodes[i] = bottomFace[i];
    }

    for (int i = 0; i < 4; i++)
    {
        int n0 = bottomFace[i];
        int n1 = bottomFace[(i + 1) % 4];

        int topNode = -1;
        for (const auto& sideFace : sideFaces)
        {
            int n0Idx = -1, n1Idx = -1;

            for (size_t j = 0; j < sideFace.size(); ++j)
            {
                if (sideFace[j] == n0)
                    n0Idx = j;
                if (sideFace[j] == n1)
                    n1Idx = j;
            }

            if (n0Idx >= 0 && n1Idx >= 0)
            {
                int nextIdx = (n0Idx + 1) % 4;
                int prevIdx = (n0Idx + 3) % 4;

                int candidate = -1;
                if (nextIdx == n1Idx)
                {
                    candidate = sideFace[prevIdx];
                }
                else if (prevIdx == n1Idx)
                {
                    candidate = sideFace[nextIdx];
                }

                if (candidate != -1)
                {
                    auto it =
                        std::find(topFace.begin(), topFace.end(), candidate);
                    if (it != topFace.end())
                    {
                        topNode = candidate;
                        break;
                    }
                }
            }
        }

        if (topNode == -1)
        {
            double minDist = 1e10;
            for (int node : topFace)
            {
                bool alreadyUsed = false;
                for (int j = 0; j < i; j++)
                {
                    if (orderedNodes[4 + j] == node)
                    {
                        alreadyUsed = true;
                        break;
                    }
                }
                if (!alreadyUsed)
                {
                    const Point& pb = points[n0];
                    const Point& pt = points[node];
                    double dx = pt.x - pb.x;
                    double dy = pt.y - pb.y;
                    double dz = pt.z - pb.z;
                    double dist = dx * dx + dy * dy + dz * dz;

                    if (dist < minDist)
                    {
                        minDist = dist;
                        topNode = node;
                    }
                }
            }
        }

        if (topNode != -1)
        {
            orderedNodes[4 + i] = topNode;
        }
        else
        {
            orderedNodes[4 + i] = topFace[i];
        }
    }

    auto computeVolume = [&](const std::vector<int>& nodes) -> double
    {
        const Point& p0 = points[nodes[0]];
        const Point& p1 = points[nodes[1]];
        const Point& p3 = points[nodes[3]];
        const Point& p4 = points[nodes[4]];

        double v1x = p1.x - p0.x, v1y = p1.y - p0.y, v1z = p1.z - p0.z;
        double v2x = p3.x - p0.x, v2y = p3.y - p0.y, v2z = p3.z - p0.z;
        double v3x = p4.x - p0.x, v3y = p4.y - p0.y, v3z = p4.z - p0.z;

        return v1x * (v2y * v3z - v2z * v3y) + v1y * (v2z * v3x - v2x * v3z) +
               v1z * (v2x * v3y - v2y * v3x);
    };

    double volume = computeVolume(orderedNodes);
    if (volume < 0)
    {
        std::swap(orderedNodes[1], orderedNodes[3]);
        std::swap(orderedNodes[5], orderedNodes[7]);
    }

    return orderedNodes;
}

std::vector<std::string>
ExodusWriter::buildCellGroups(const std::vector<Cell>& cells,
                             const std::vector<CellZone>& cellZones) const
{
    // Mirror the standard-block grouping in writeElements: with >1 zone, cells
    // are grouped by zone name (or "unzoned"); otherwise a single "fluid" group.
    std::vector<std::string> grp(cells.size(), "fluid");
    if (cellZones.size() > 1)
    {
        std::fill(grp.begin(), grp.end(), "unzoned");
        for (const auto& z : cellZones)
            for (int ci : z.cellIndices)
                if (ci >= 0 && ci < (int)cells.size())
                    grp[ci] = z.name;
    }
    return grp;
}

void ExodusWriter::buildPolyDecomposition(const std::vector<Point>& points,
                                          const std::vector<Face>& faces,
                                          const std::vector<Cell>& cells,
                                          const std::vector<int>& owner,
                                          int boundaryStart)
{
    polyExtraPoints.clear();
    polySubElems.clear();
    polySubElemCell.clear();
    polySubElemExoId.clear();
    polyFaceToSubs.clear();
    cellDecomposed.assign(cells.size(), 0);

    const int nBase = (int)points.size();

    // Coordinate of a node index in the combined [base ; extra] space. Returns
    // by value so it stays valid across push_backs into polyExtraPoints.
    auto pointAt = [&](int idx) -> Point
    {
        return (idx < nBase) ? points[idx] : polyExtraPoints[idx - nBase];
    };

    // Signed volume (x6) of tet (a,b,c,d); >0 for positive Exodus orientation.
    auto tetVol = [&](int a, int b, int c, int d) -> double
    {
        Point p0 = pointAt(a), p1 = pointAt(b), p2 = pointAt(c), p3 = pointAt(d);
        double v1x = p1.x - p0.x, v1y = p1.y - p0.y, v1z = p1.z - p0.z;
        double v2x = p2.x - p0.x, v2y = p2.y - p0.y, v2z = p2.z - p0.z;
        double v3x = p3.x - p0.x, v3y = p3.y - p0.y, v3z = p3.z - p0.z;
        return v1x * (v2y * v3z - v2z * v3y) + v1y * (v2z * v3x - v2x * v3z) +
               v1z * (v2x * v3y - v2y * v3x);
    };

    // Divergence-theorem contribution of an outward triangle, used to get the
    // signed volume of a standard element from its Exodus (outward) faces.
    auto triDiv = [&](int a, int b, int c) -> double
    {
        Point A = pointAt(a), B = pointAt(b), C = pointAt(c);
        return A.x * (B.y * C.z - B.z * C.y) + A.y * (B.z * C.x - B.x * C.z) +
               A.z * (B.x * C.y - B.y * C.x);
    };
    // Signed volume of a standard element given its ordered nodes; <=0 means
    // the ordering is invalid/inverted and the cell should be decomposed.
    auto stdVol = [&](const std::vector<int>& n, const std::string& t) -> double
    {
        double v = 0.0;
        auto Q = [&](int a, int b, int c, int d)
        { v += triDiv(n[a], n[b], n[c]) + triDiv(n[a], n[c], n[d]); };
        auto T = [&](int a, int b, int c) { v += triDiv(n[a], n[b], n[c]); };
        if (t == "hex")
        {
            Q(0, 1, 5, 4); Q(1, 2, 6, 5); Q(2, 3, 7, 6);
            Q(3, 0, 4, 7); Q(0, 3, 2, 1); Q(4, 5, 6, 7);
        }
        else if (t == "tet")
        {
            T(0, 1, 3); T(1, 2, 3); T(2, 0, 3); T(0, 2, 1);
        }
        else if (t == "pyr")
        {
            T(0, 1, 4); T(1, 2, 4); T(2, 3, 4); T(3, 0, 4); Q(0, 3, 2, 1);
        }
        else if (t == "wedge")
        {
            Q(0, 1, 4, 3); Q(1, 2, 5, 4); Q(2, 0, 3, 5); T(0, 2, 1); T(3, 4, 5);
        }
        return v / 6.0;
    };

    // Face centroids are shared per face so both adjacent poly cells fan an
    // n-gon identically, keeping the split conformal.
    std::unordered_map<int, int> faceCentroidNode;

    for (int c = 0; c < (int)cells.size(); ++c)
    {
        const auto& cell = cells[c];

        // Decompose "unknown" cells, and also any standard cell whose ordering
        // would yield a non-positive volume (mirrors what writeElements emits,
        // so those otherwise-inverted elements are rescued as valid sub-cells).
        if (cell.type != "unknown")
        {
            std::vector<int> on;
            int need = 8;
            if (cell.type == "hex")
                on = orderHexNodes(cell, faces, points);
            else if (cell.type == "tet")
                on = orderTetNodes(cell, faces, points), need = 4;
            else if (cell.type == "pyr")
                on = orderPyramidNodes(cell, faces, points), need = 5;
            else if (cell.type == "wedge")
                on = orderWedgeNodes(cell, faces, points), need = 6;
            if ((int)on.size() == need && stdVol(on, cell.type) > 0.0)
                continue;  // valid standard element; leave it as-is
        }
        cellDecomposed[c] = 1;

        std::set<int> verts;
        for (int fi : cell.faceIndices)
            if (fi >= 0 && fi < (int)faces.size())
                for (int n : faces[fi].pointIndices)
                    verts.insert(n);
        if (verts.empty())
            continue;

        // Cell centroid (apex of every sub-element) = mean of cell vertices.
        double cx = 0, cy = 0, cz = 0;
        for (int n : verts)
        {
            Point p = pointAt(n);
            cx += p.x;
            cy += p.y;
            cz += p.z;
        }
        cx /= verts.size();
        cy /= verts.size();
        cz /= verts.size();
        int cNode = nBase + (int)polyExtraPoints.size();
        polyExtraPoints.push_back({cx, cy, cz});

        for (int fi : cell.faceIndices)
        {
            if (fi < 0 || fi >= (int)faces.size())
                continue;
            const auto& fv = faces[fi].pointIndices;
            int nv = (int)fv.size();
            bool onBoundary = (fi >= boundaryStart) && (owner[fi] == c);

            if (nv == 3)
            {
                int a = fv[0], b = fv[1], cc = fv[2];
                if (tetVol(a, b, cc, cNode) < 0)
                    std::swap(b, cc);
                int gid = (int)polySubElems.size();
                polySubElems.push_back({'T', {a, b, cc, cNode}});
                if (onBoundary)
                    polyFaceToSubs[fi].push_back(
                        {gid, getTetFaceId(fv, polySubElems[gid].nodes)});
            }
            else if (nv == 4)
            {
                int q0 = fv[0], q1 = fv[1], q2 = fv[2], q3 = fv[3];
                // Orient the base from the true split volume (same q0-q2
                // diagonal the pyramid is later evaluated on), which is robust
                // for warped quads unlike a single-triangle normal.
                if (tetVol(q0, q1, q2, cNode) + tetVol(q0, q2, q3, cNode) < 0)
                    std::swap(q1, q3);
                int gid = (int)polySubElems.size();
                polySubElems.push_back({'P', {q0, q1, q2, q3, cNode}});
                if (onBoundary)
                    polyFaceToSubs[fi].push_back(
                        {gid, getPyramidFaceId(fv, polySubElems[gid].nodes)});
            }
            else if (nv >= 5)
            {
                int fNode;
                auto it = faceCentroidNode.find(fi);
                if (it == faceCentroidNode.end())
                {
                    double fx = 0, fy = 0, fz = 0;
                    for (int n : fv)
                    {
                        Point p = pointAt(n);
                        fx += p.x;
                        fy += p.y;
                        fz += p.z;
                    }
                    fNode = nBase + (int)polyExtraPoints.size();
                    polyExtraPoints.push_back({fx / nv, fy / nv, fz / nv});
                    faceCentroidNode[fi] = fNode;
                }
                else
                {
                    fNode = it->second;
                }
                for (int k = 0; k < nv; ++k)
                {
                    int a = fv[k], b = fv[(k + 1) % nv];
                    int y = b, z = fNode;
                    if (tetVol(a, y, z, cNode) < 0)
                        std::swap(y, z);
                    int gid = (int)polySubElems.size();
                    polySubElems.push_back({'T', {a, y, z, cNode}});
                    if (onBoundary)
                    {
                        std::vector<int> tri = {a, b, fNode};
                        polyFaceToSubs[fi].push_back(
                            {gid, getTetFaceId(tri, polySubElems[gid].nodes)});
                    }
                }
            }
        }
        // All sub-elements just appended belong to this cell.
        polySubElemCell.resize(polySubElems.size(), c);
    }

    polySubElemExoId.assign(polySubElems.size(), 0);
}

void ExodusWriter::writeElements(const OpenFOAMMeshReader& reader)
{
    const auto& cells = reader.getCells();
    const auto& faces = reader.getFaces();
    const auto& points = reader.getPoints();
    const auto& cellZones = reader.getCellZones();

    nc_redef(ncid);

    struct BlockInfo
    {
        std::string name;
        std::vector<int> cellIndices;
        std::string cellType;
        // Decomposed polyhedral blocks carry sub-element indices instead of
        // cell indices; connectivity is taken from polySubElems.
        bool decomposed = false;
        std::vector<int> subIndices;
    };

    std::vector<BlockInfo> blocks;

    // "unknown" cells never form a block of their own; they are decomposed
    // into the tet/pyramid blocks appended below.
    if (!cellZones.empty() && cellZones.size() > 1)
    {
        std::set<int> zonedCells;

        for (const auto& zone : cellZones)
        {
            std::map<std::string, std::vector<int>> zoneElemsByType;

            for (int cellIdx : zone.cellIndices)
            {
                if (cellIdx < cells.size())
                {
                    if (!cellDecomposed[cellIdx])
                        zoneElemsByType[cells[cellIdx].type].push_back(cellIdx);
                    zonedCells.insert(cellIdx);
                }
            }

            for (const auto& [cellType, cellIndices] : zoneElemsByType)
            {
                BlockInfo block;
                std::string originalName = zone.name + "-" + cellType;
                block.name = getBlockName(originalName);
                block.cellIndices = cellIndices;
                block.cellType = cellType;
                blocks.push_back(block);
            }
        }

        std::map<std::string, std::vector<int>> unzonedElemsByType;
        for (size_t i = 0; i < cells.size(); ++i)
        {
            if (zonedCells.find(i) == zonedCells.end() && !cellDecomposed[i])
            {
                unzonedElemsByType[cells[i].type].push_back(i);
            }
        }

        for (const auto& [cellType, cellIndices] : unzonedElemsByType)
        {
            BlockInfo block;
            std::string originalName = "unzoned-" + cellType;
            block.name = getBlockName(originalName);
            block.cellIndices = cellIndices;
            block.cellType = cellType;
            blocks.push_back(block);
        }
    }
    else
    {
        std::map<std::string, std::vector<int>> elemsByType;
        for (size_t i = 0; i < cells.size(); ++i)
        {
            if (!cellDecomposed[i])
                elemsByType[cells[i].type].push_back(i);
        }

        for (const auto& [cellType, cellIndices] : elemsByType)
        {
            BlockInfo block;
            std::string originalName = "fluid-" + cellType;
            block.name = getBlockName(originalName);
            block.cellIndices = cellIndices;
            block.cellType = cellType;
            blocks.push_back(block);
        }
    }

    // Append decomposed tet/pyramid blocks grouped by region (zone/mesh), so a
    // rotor poly block never mixes with a stator one.
    {
        std::vector<std::string> grp = buildCellGroups(cells, cellZones);
        std::map<std::string, std::vector<int>> tetByGrp, pyrByGrp;
        for (int g = 0; g < (int)polySubElems.size(); ++g)
        {
            const std::string& k = grp[polySubElemCell[g]];
            (polySubElems[g].type == 'T' ? tetByGrp : pyrByGrp)[k].push_back(g);
        }
        for (auto& [k, gids] : tetByGrp)
        {
            BlockInfo block;
            block.name = getBlockName(k + "-poly-tet");
            block.cellType = "tet";
            block.decomposed = true;
            block.subIndices = std::move(gids);
            blocks.push_back(std::move(block));
        }
        for (auto& [k, gids] : pyrByGrp)
        {
            BlockInfo block;
            block.name = getBlockName(k + "-poly-pyr");
            block.cellType = "pyr";
            block.decomposed = true;
            block.subIndices = std::move(gids);
            blocks.push_back(std::move(block));
        }
    }

    int blockId = 1;
    for (const auto& block : blocks)
    {
        int numElemsInBlock =
            block.decomposed ? (int)block.subIndices.size()
                             : (int)block.cellIndices.size();

        int numNodesPerElem = 8;
        std::string exoType = "HEX8";

        if (block.cellType == "tet")
        {
            numNodesPerElem = 4;
            exoType = "TETRA4";
        }
        else if (block.cellType == "pyr")
        {
            numNodesPerElem = 5;
            exoType = "PYRAMID5";
        }
        else if (block.cellType == "wedge")
        {
            numNodesPerElem = 6;
            exoType = "WEDGE6";
        }
        else if (block.cellType == "hex")
        {
            numNodesPerElem = 8;
            exoType = "HEX8";
        }

        int dim_num_el_in_blk, dim_num_nod_per_el;
        std::string blk_num = std::to_string(blockId);

        int status;
        status = nc_def_dim(ncid,
                            ("num_el_in_blk" + blk_num).c_str(),
                            numElemsInBlock,
                            &dim_num_el_in_blk);
        checkError(status, "Failed to define element block dimension");

        status = nc_def_dim(ncid,
                            ("num_nod_per_el" + blk_num).c_str(),
                            numNodesPerElem,
                            &dim_num_nod_per_el);
        checkError(status, "Failed to define nodes per element dimension");

        int var_connect;
        int dims_connect[2] = {dim_num_el_in_blk, dim_num_nod_per_el};
        status = nc_def_var(ncid,
                            ("connect" + blk_num).c_str(),
                            NC_INT,
                            2,
                            dims_connect,
                            &var_connect);
        checkError(status, "Failed to define connectivity variable");

        status = nc_put_att_text(
            ncid, var_connect, "elem_type", exoType.length(), exoType.c_str());
        checkError(status, "Failed to set element type attribute");

        blockId++;
    }

    nc_enddef(ncid);

    cellToExodusElem.assign(cells.size(), 0);
    int nextExodusElem = 1;

    blockId = 1;
    for (const auto& block : blocks)
    {
        std::vector<int> connectivity;

        if (block.decomposed)
        {
            for (int gid : block.subIndices)
            {
                polySubElemExoId[gid] = nextExodusElem++;
                for (int n : polySubElems[gid].nodes)
                    connectivity.push_back(n + 1);
            }
        }
        else
            for (int cellIdx : block.cellIndices)
            {
                cellToExodusElem[cellIdx] = nextExodusElem++;

                const auto& cell = cells[cellIdx];

                std::vector<int> nodes;
                if (block.cellType == "hex")
                {
                    nodes = orderHexNodes(cell, faces, points);
                }
                else if (block.cellType == "tet")
                {
                    nodes = orderTetNodes(cell, faces, points);
                }
                else if (block.cellType == "pyr")
                {
                    nodes = orderPyramidNodes(cell, faces, points);
                }
                else if (block.cellType == "wedge")
                {
                    nodes = orderWedgeNodes(cell, faces, points);
                }

                int numNodesPerElem = (block.cellType == "tet")     ? 4
                                      : (block.cellType == "pyr")   ? 5
                                      : (block.cellType == "wedge") ? 6
                                                                    : 8;
                if ((int)nodes.size() != numNodesPerElem)
                {
                    std::cerr << "Warning: cell " << cellIdx
                              << " (type=" << block.cellType << ") returned "
                              << nodes.size() << " nodes, expected "
                              << numNodesPerElem << std::endl;
                }
                for (int i = 0; i < numNodesPerElem; ++i)
                {
                    connectivity.push_back(i < (int)nodes.size() ? nodes[i] + 1
                                                                 : 1);
                }
            }

        int var_id;
        nc_inq_varid(
            ncid, ("connect" + std::to_string(blockId)).c_str(), &var_id);
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

        std::cout << "Wrote element block " << blockId << " (" << block.name
                  << ") with "
                  << (block.decomposed ? block.subIndices.size()
                                       : block.cellIndices.size())
                  << " elements"
                  << std::endl;

        blockId++;
    }
}

std::vector<int> ExodusWriter::orderWedgeNodes(const Cell& cell,
                                               const std::vector<Face>& faces,
                                               const std::vector<Point>& points)
{
    if (cell.faceIndices.size() < 5)
    {
        return std::vector<int>(6, 0);
    }

    std::vector<int> triFaceIndices;
    std::vector<int> quadFaceIndices;
    for (int faceIdx : cell.faceIndices)
    {
        if (faceIdx >= 0 && faceIdx < (int)faces.size())
        {
            int nPts = faces[faceIdx].pointIndices.size();
            if (nPts == 3)
                triFaceIndices.push_back(faceIdx);
            else if (nPts == 4)
                quadFaceIndices.push_back(faceIdx);
        }
    }

    if (triFaceIndices.size() != 2 || quadFaceIndices.size() != 3)
    {
        std::set<int> nodeSet;
        for (int fi : cell.faceIndices)
        {
            if (fi >= 0 && fi < (int)faces.size())
                for (int n : faces[fi].pointIndices)
                    nodeSet.insert(n);
        }
        std::vector<int> nodes(nodeSet.begin(), nodeSet.end());
        nodes.resize(6, 0);
        return nodes;
    }

    // Pick the tri face with lower mean-Z as "bottom" for consistent
    // orientation
    auto triCentroidZ = [&](int fi)
    {
        double z = 0.0;
        for (int n : faces[fi].pointIndices)
            z += points[n].z;
        return z / 3.0;
    };
    int bottomTriIdx =
        (triCentroidZ(triFaceIndices[0]) <= triCentroidZ(triFaceIndices[1]))
            ? triFaceIndices[0]
            : triFaceIndices[1];
    int topTriIdx = (bottomTriIdx == triFaceIndices[0]) ? triFaceIndices[1]
                                                        : triFaceIndices[0];

    const std::vector<int>& bottomFace = faces[bottomTriIdx].pointIndices;
    const std::vector<int>& topFace = faces[topTriIdx].pointIndices;

    std::vector<int> orderedNodes(6);
    orderedNodes[0] = bottomFace[0];
    orderedNodes[1] = bottomFace[1];
    orderedNodes[2] = bottomFace[2];

    // Map each bottom node to its corresponding top node via shared quad side
    // faces
    for (int i = 0; i < 3; ++i)
    {
        int n0 = bottomFace[i];
        int n1 = bottomFace[(i + 1) % 3];

        int topNode = -1;
        for (int qfi : quadFaceIndices)
        {
            const auto& qf = faces[qfi].pointIndices;
            int n0Idx = -1, n1Idx = -1;
            for (int j = 0; j < 4; ++j)
            {
                if (qf[j] == n0)
                    n0Idx = j;
                if (qf[j] == n1)
                    n1Idx = j;
            }
            if (n0Idx >= 0 && n1Idx >= 0)
            {
                int nextIdx = (n0Idx + 1) % 4;
                int prevIdx = (n0Idx + 3) % 4;
                int candidate = -1;
                if (nextIdx == n1Idx)
                    candidate = qf[prevIdx];
                else if (prevIdx == n1Idx)
                    candidate = qf[nextIdx];
                if (candidate != -1)
                {
                    auto it =
                        std::find(topFace.begin(), topFace.end(), candidate);
                    if (it != topFace.end())
                    {
                        topNode = candidate;
                        break;
                    }
                }
            }
        }

        if (topNode == -1)
        {
            // Fallback: nearest unused top node to n0
            double minDist = 1e30;
            for (int node : topFace)
            {
                bool used = false;
                for (int j = 0; j < i; ++j)
                    if (orderedNodes[3 + j] == node)
                    {
                        used = true;
                        break;
                    }
                if (!used)
                {
                    const Point& pb = points[n0];
                    const Point& pt = points[node];
                    double dist = (pt.x - pb.x) * (pt.x - pb.x) +
                                  (pt.y - pb.y) * (pt.y - pb.y) +
                                  (pt.z - pb.z) * (pt.z - pb.z);
                    if (dist < minDist)
                    {
                        minDist = dist;
                        topNode = node;
                    }
                }
            }
        }
        orderedNodes[3 + i] = (topNode != -1) ? topNode : topFace[i];
    }

    // Orientation check: normal of bottom triangle should point toward top
    // nodes. n = (p1-p0) x (p2-p0);  h = centroid(top) - centroid(bottom)
    const Point& p0 = points[orderedNodes[0]];
    const Point& p1 = points[orderedNodes[1]];
    const Point& p2 = points[orderedNodes[2]];
    const Point& p3 = points[orderedNodes[3]];
    const Point& p4 = points[orderedNodes[4]];
    const Point& p5 = points[orderedNodes[5]];
    double v1x = p1.x - p0.x, v1y = p1.y - p0.y, v1z = p1.z - p0.z;
    double v2x = p2.x - p0.x, v2y = p2.y - p0.y, v2z = p2.z - p0.z;
    double nx = v1y * v2z - v1z * v2y;
    double ny = v1z * v2x - v1x * v2z;
    double nz = v1x * v2y - v1y * v2x;
    double hx = (p3.x + p4.x + p5.x) / 3.0 - (p0.x + p1.x + p2.x) / 3.0;
    double hy = (p3.y + p4.y + p5.y) / 3.0 - (p0.y + p1.y + p2.y) / 3.0;
    double hz = (p3.z + p4.z + p5.z) / 3.0 - (p0.z + p1.z + p2.z) / 3.0;
    if (nx * hx + ny * hy + nz * hz < 0)
    {
        std::swap(orderedNodes[1], orderedNodes[2]);
        std::swap(orderedNodes[4], orderedNodes[5]);
    }

    return orderedNodes;
}

int ExodusWriter::getHexFaceId(const std::vector<int>& faceNodes,
                               const std::vector<int>& hexNodes)
{
    std::set<int> faceSet(faceNodes.begin(), faceNodes.end());

    std::vector<std::vector<int>> hexFaces = {
        {hexNodes[0], hexNodes[1], hexNodes[5], hexNodes[4]},
        {hexNodes[1], hexNodes[2], hexNodes[6], hexNodes[5]},
        {hexNodes[2], hexNodes[3], hexNodes[7], hexNodes[6]},
        {hexNodes[3], hexNodes[0], hexNodes[4], hexNodes[7]},
        {hexNodes[0], hexNodes[3], hexNodes[2], hexNodes[1]},
        {hexNodes[4], hexNodes[5], hexNodes[6], hexNodes[7]}};

    for (int i = 0; i < 6; ++i)
    {
        std::set<int> hexFaceSet(hexFaces[i].begin(), hexFaces[i].end());
        if (faceSet == hexFaceSet)
        {
            return i + 1;
        }
    }

    return 1;
}

int ExodusWriter::getTetFaceId(const std::vector<int>& faceNodes,
                               const std::vector<int>& tetNodes)
{
    std::set<int> faceSet(faceNodes.begin(), faceNodes.end());
    // Exodus II TETRA4 face ordering
    std::vector<std::vector<int>> tetFaces = {
        {tetNodes[0], tetNodes[1], tetNodes[3]},
        {tetNodes[1], tetNodes[2], tetNodes[3]},
        {tetNodes[2], tetNodes[0], tetNodes[3]},
        {tetNodes[0], tetNodes[2], tetNodes[1]}};
    for (int i = 0; i < 4; ++i)
    {
        std::set<int> fs(tetFaces[i].begin(), tetFaces[i].end());
        if (faceSet == fs)
            return i + 1;
    }
    return 1;
}

int ExodusWriter::getPyramidFaceId(const std::vector<int>& faceNodes,
                                   const std::vector<int>& pyrNodes)
{
    std::set<int> faceSet(faceNodes.begin(), faceNodes.end());
    // Exodus II PYRAMID5 face ordering
    std::vector<std::vector<int>> pyrFaces = {
        {pyrNodes[0], pyrNodes[1], pyrNodes[4]},
        {pyrNodes[1], pyrNodes[2], pyrNodes[4]},
        {pyrNodes[2], pyrNodes[3], pyrNodes[4]},
        {pyrNodes[3], pyrNodes[0], pyrNodes[4]},
        {pyrNodes[0], pyrNodes[3], pyrNodes[2], pyrNodes[1]}};
    for (int i = 0; i < 5; ++i)
    {
        std::set<int> fs(pyrFaces[i].begin(), pyrFaces[i].end());
        if (faceSet == fs)
            return i + 1;
    }
    return 1;
}

int ExodusWriter::getWedgeFaceId(const std::vector<int>& faceNodes,
                                 const std::vector<int>& wedgeNodes)
{
    std::set<int> faceSet(faceNodes.begin(), faceNodes.end());
    // Exodus II WEDGE6 face ordering
    std::vector<std::vector<int>> wedgeFaces = {
        {wedgeNodes[0], wedgeNodes[1], wedgeNodes[4], wedgeNodes[3]},
        {wedgeNodes[1], wedgeNodes[2], wedgeNodes[5], wedgeNodes[4]},
        {wedgeNodes[2], wedgeNodes[0], wedgeNodes[3], wedgeNodes[5]},
        {wedgeNodes[0], wedgeNodes[2], wedgeNodes[1]},
        {wedgeNodes[3], wedgeNodes[4], wedgeNodes[5]}};
    for (int i = 0; i < 5; ++i)
    {
        std::set<int> fs(wedgeFaces[i].begin(), wedgeFaces[i].end());
        if (faceSet == fs)
            return i + 1;
    }
    return 1;
}

void ExodusWriter::writeSideSets(const OpenFOAMMeshReader& reader)
{
    const auto& patches = reader.getBoundaryPatches();
    const auto& faces = reader.getFaces();
    const auto& owner = reader.getOwner();
    const auto& cells = reader.getCells();
    const auto& points = reader.getPoints();

    if (patches.empty())
    {
        std::cout << "No boundary patches to write as sidesets" << std::endl;
        return;
    }

    // Ordered nodes for standard cells only; polyhedral cells resolve their
    // boundary faces through polyFaceToSubs instead.
    std::map<int, std::vector<int>> cellToOrderedNodes;
    for (size_t i = 0; i < cells.size(); ++i)
    {
        const auto& cell = cells[i];
        if (cell.type == "hex")
            cellToOrderedNodes[i] = orderHexNodes(cell, faces, points);
        else if (cell.type == "tet")
            cellToOrderedNodes[i] = orderTetNodes(cell, faces, points);
        else if (cell.type == "pyr")
            cellToOrderedNodes[i] = orderPyramidNodes(cell, faces, points);
        else if (cell.type == "wedge")
            cellToOrderedNodes[i] = orderWedgeNodes(cell, faces, points);
    }

    // Build each patch's (elem, side) entries first. A boundary face of a
    // polyhedral cell maps to the base side(s) of its sub-element(s); for an
    // n-gon face that is several sides, so the per-sideset count is not known
    // until the faces are walked.
    std::vector<std::vector<int>> patchElem(patches.size());
    std::vector<std::vector<int>> patchSide(patches.size());
    for (size_t i = 0; i < patches.size(); ++i)
    {
        const auto& patch = patches[i];
        auto& elem_list = patchElem[i];
        auto& side_list = patchSide[i];

        for (int j = 0; j < patch.nFaces; ++j)
        {
            int faceIdx = patch.startFace + j;
            if (faceIdx >= (int)owner.size())
                continue;

            auto pit = polyFaceToSubs.find(faceIdx);
            if (pit != polyFaceToSubs.end())
            {
                for (const auto& entry : pit->second)
                {
                    elem_list.push_back(polySubElemExoId[entry.first]);
                    side_list.push_back(entry.second);
                }
                continue;
            }

            int cellIdx = owner[faceIdx];
            int exoId =
                (cellIdx >= 0 && cellIdx < (int)cellToExodusElem.size())
                    ? cellToExodusElem[cellIdx]
                    : (cellIdx + 1);

            const auto& face = faces[faceIdx];
            const auto& cell = cells[cellIdx];
            int sideId = 1;
            auto cit = cellToOrderedNodes.find(cellIdx);
            if (cit != cellToOrderedNodes.end())
            {
                const auto& ordNodes = cit->second;
                if (cell.type == "hex")
                    sideId = getHexFaceId(face.pointIndices, ordNodes);
                else if (cell.type == "tet")
                    sideId = getTetFaceId(face.pointIndices, ordNodes);
                else if (cell.type == "pyr")
                    sideId = getPyramidFaceId(face.pointIndices, ordNodes);
                else if (cell.type == "wedge")
                    sideId = getWedgeFaceId(face.pointIndices, ordNodes);
            }
            elem_list.push_back(exoId);
            side_list.push_back(sideId);
        }
    }

    nc_redef(ncid);

    for (size_t i = 0; i < patches.size(); ++i)
    {
        const auto& patch = patches[i];
        std::string ss_num = std::to_string(i + 1);

        int dim_num_side_ss;
        nc_def_dim(ncid,
                   ("num_side_ss" + ss_num).c_str(),
                   patchElem[i].size(),
                   &dim_num_side_ss);

        int var_elem_ss, var_side_ss;
        nc_def_var(ncid,
                   ("elem_ss" + ss_num).c_str(),
                   NC_INT,
                   1,
                   &dim_num_side_ss,
                   &var_elem_ss);
        nc_def_var(ncid,
                   ("side_ss" + ss_num).c_str(),
                   NC_INT,
                   1,
                   &dim_num_side_ss,
                   &var_side_ss);

        std::string sidesetName = getSidesetName(patch.name);
        nc_put_att_text(ncid,
                        var_elem_ss,
                        "name",
                        sidesetName.length(),
                        sidesetName.c_str());
    }

    int var_ss_names;
    int dim_num_side_sets, dim_len_name;
    int status;
    status = nc_inq_dimid(ncid, "num_side_sets", &dim_num_side_sets);
    checkError(status, "Failed to get num_side_sets dimension");

    status = nc_inq_dimid(ncid, "len_name", &dim_len_name);
    checkError(status, "Failed to get len_name dimension");

    int dims_ss_names[2] = {dim_num_side_sets, dim_len_name};
    status =
        nc_def_var(ncid, "ss_names", NC_CHAR, 2, dims_ss_names, &var_ss_names);
    checkError(status, "Failed to define ss_names variable");

    nc_enddef(ncid);

    for (size_t i = 0; i < patches.size(); ++i)
    {
        const auto& patch = patches[i];
        std::string ss_num = std::to_string(i + 1);

        const std::vector<int>& elem_list = patchElem[i];
        const std::vector<int>& side_list = patchSide[i];

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

        std::string sidesetName = getSidesetName(patch.name);
        char name_padded[33];
        std::memset(name_padded, 0, 33);
        std::strncpy(name_padded, sidesetName.c_str(), 32);

        nc_inq_varid(ncid, "ss_names", &var_id);
        size_t start[2] = {i, 0};
        size_t count[2] = {1, 33};
        nc_put_vara_text(ncid, var_id, start, count, name_padded);

        std::cout << "Wrote sideset " << (i + 1) << ": " << sidesetName
                  << " with " << elem_list.size() << " sides ("
                  << patch.nFaces << " boundary faces)" << std::endl;
    }
}

void ExodusWriter::writeMesh(const OpenFOAMMeshReader& reader)
{
    customElementBlockNames.clear();
    customSidesetNames.clear();

    const auto& cells = reader.getCells();
    const auto& cellZones = reader.getCellZones();
    const auto& patches = reader.getBoundaryPatches();

    // Split polyhedral cells into conformal tets/pyramids before counting, as
    // it adds centroid nodes and replaces each poly cell with sub-elements.
    int boundaryStart = (int)reader.getFaces().size();
    for (const auto& p : patches)
        boundaryStart = std::min(boundaryStart, p.startFace);
    buildPolyDecomposition(reader.getPoints(),
                           reader.getFaces(),
                           cells,
                           reader.getOwner(),
                           boundaryStart);

    int nStandardCells = 0;
    for (size_t i = 0; i < cells.size(); ++i)
        if (!cellDecomposed[i])
            ++nStandardCells;

    int numNodes = reader.getNumPoints() + (int)polyExtraPoints.size();
    int numElems = nStandardCells + (int)polySubElems.size();
    int numNodeSets = 0;
    int numSideSets = reader.getNumBoundaryPatches();

    int numElemBlocks = 0;
    if (!cellZones.empty())
    {
        std::set<int> zonedCells;

        for (const auto& zone : cellZones)
        {
            std::set<std::string> typesInZone;
            for (int cellIdx : zone.cellIndices)
            {
                if (cellIdx < cells.size())
                {
                    if (!cellDecomposed[cellIdx])
                        typesInZone.insert(cells[cellIdx].type);
                    zonedCells.insert(cellIdx);
                }
            }
            numElemBlocks += typesInZone.size();
        }

        std::set<std::string> unzonedTypes;
        for (size_t i = 0; i < cells.size(); ++i)
        {
            if (zonedCells.find(i) == zonedCells.end() && !cellDecomposed[i])
            {
                unzonedTypes.insert(cells[i].type);
            }
        }
        numElemBlocks += unzonedTypes.size();
    }
    else
    {
        std::set<std::string> types;
        for (size_t i = 0; i < cells.size(); ++i)
            if (!cellDecomposed[i])
                types.insert(cells[i].type);
        numElemBlocks = (int)types.size();
    }

    // One poly block per (region, tet/pyr) pair present; must match the block
    // grouping in writeElements.
    {
        std::vector<std::string> grp = buildCellGroups(cells, cellZones);
        std::set<std::string> polyKeys;
        for (int g = 0; g < (int)polySubElems.size(); ++g)
            polyKeys.insert(grp[polySubElemCell[g]] +
                            (polySubElems[g].type == 'T' ? "\tT" : "\tP"));
        numElemBlocks += (int)polyKeys.size();
    }

    initializeExodusFile(
        numNodes, numElems, numElemBlocks, numNodeSets, numSideSets);
    writeNodes(reader.getPoints());
    writeElements(reader);
    writeSideSets(reader);

    std::cout << "Successfully wrote Exodus II file: " << filename << std::endl;
}

void ExodusWriter::writeMesh(
    const OpenFOAMMeshReader& reader,
    const std::map<std::string, std::string>& elementBlockNames,
    const std::map<std::string, std::string>& sidesetNames)
{
    customElementBlockNames = elementBlockNames;
    customSidesetNames = sidesetNames;
    writeMesh(reader);
}

std::string ExodusWriter::getBlockName(const std::string& originalName)
{
    auto it = customElementBlockNames.find(originalName);
    if (it != customElementBlockNames.end() && !it->second.empty())
    {
        return it->second;
    }
    return originalName;
}

std::string ExodusWriter::getSidesetName(const std::string& originalName)
{
    auto it = customSidesetNames.find(originalName);
    if (it != customSidesetNames.end() && !it->second.empty())
    {
        return it->second;
    }
    return originalName;
}

// MergedMeshReader overloads - delegate to template implementations
void ExodusWriter::writeMesh(const MergedMeshReader& reader)
{
    customElementBlockNames.clear();
    customSidesetNames.clear();

    const auto& cells = reader.getCells();
    const auto& cellZones = reader.getCellZones();
    const auto& patches = reader.getBoundaryPatches();

    // Split polyhedral cells into conformal tets/pyramids before counting.
    int boundaryStart = (int)reader.getFaces().size();
    for (const auto& p : patches)
        boundaryStart = std::min(boundaryStart, p.startFace);
    buildPolyDecomposition(reader.getPoints(),
                           reader.getFaces(),
                           cells,
                           reader.getOwner(),
                           boundaryStart);

    int nStandardCells = 0;
    for (size_t i = 0; i < cells.size(); ++i)
        if (!cellDecomposed[i])
            ++nStandardCells;

    int numNodes = reader.getNumPoints() + (int)polyExtraPoints.size();
    int numElems = nStandardCells + (int)polySubElems.size();
    int numNodeSets = 0;
    int numSideSets = reader.getNumBoundaryPatches();

    int numElemBlocks = 0;
    if (!cellZones.empty())
    {
        std::set<int> zonedCells;

        for (const auto& zone : cellZones)
        {
            std::set<std::string> typesInZone;
            for (int cellIdx : zone.cellIndices)
            {
                if (cellIdx < cells.size())
                {
                    if (!cellDecomposed[cellIdx])
                        typesInZone.insert(cells[cellIdx].type);
                }
            }
            numElemBlocks += typesInZone.size();
            for (int cellIdx : zone.cellIndices)
            {
                zonedCells.insert(cellIdx);
            }
        }

        std::set<std::string> unzonedTypes;
        for (size_t i = 0; i < cells.size(); ++i)
        {
            if (zonedCells.find(i) == zonedCells.end() && !cellDecomposed[i])
            {
                unzonedTypes.insert(cells[i].type);
            }
        }
        numElemBlocks += unzonedTypes.size();
    }
    else
    {
        std::set<std::string> types;
        for (size_t i = 0; i < cells.size(); ++i)
        {
            if (!cellDecomposed[i])
                types.insert(cells[i].type);
        }
        numElemBlocks = types.size();
    }

    // One poly block per (region, tet/pyr) pair present; must match the block
    // grouping in writeElements.
    {
        std::vector<std::string> grp = buildCellGroups(cells, cellZones);
        std::set<std::string> polyKeys;
        for (int g = 0; g < (int)polySubElems.size(); ++g)
            polyKeys.insert(grp[polySubElemCell[g]] +
                            (polySubElems[g].type == 'T' ? "\tT" : "\tP"));
        numElemBlocks += (int)polyKeys.size();
    }

    initializeExodusFile(
        numNodes, numElems, numElemBlocks, numNodeSets, numSideSets);

    const auto& points = reader.getPoints();
    writeNodes(points);

    writeElements(reader);
    writeSideSets(reader);

    int status = nc_close(ncid);
    checkError(status, "Failed to close Exodus file");
    ncid = -1;
}

void ExodusWriter::writeMesh(
    const MergedMeshReader& reader,
    const std::map<std::string, std::string>& elementBlockNames,
    const std::map<std::string, std::string>& sidesetNames)
{
    customElementBlockNames = elementBlockNames;
    customSidesetNames = sidesetNames;
    writeMesh(reader);
}

void ExodusWriter::writeElements(const MergedMeshReader& reader)
{
    const auto& cells = reader.getCells();
    const auto& faces = reader.getFaces();
    const auto& points = reader.getPoints();
    const auto& cellZones = reader.getCellZones();

    nc_redef(ncid);

    struct BlockInfo
    {
        std::string name;
        std::vector<int> cellIndices;
        std::string cellType;
        // Decomposed polyhedral blocks carry sub-element indices instead of
        // cell indices; connectivity is taken from polySubElems.
        bool decomposed = false;
        std::vector<int> subIndices;
    };

    std::vector<BlockInfo> blocks;

    // "unknown" cells never form a block of their own; they are decomposed
    // into the tet/pyramid blocks appended below.
    if (!cellZones.empty() && cellZones.size() > 1)
    {
        std::set<int> zonedCells;

        for (const auto& zone : cellZones)
        {
            std::map<std::string, std::vector<int>> zoneElemsByType;

            for (int cellIdx : zone.cellIndices)
            {
                if (cellIdx < cells.size())
                {
                    if (!cellDecomposed[cellIdx])
                        zoneElemsByType[cells[cellIdx].type].push_back(cellIdx);
                    zonedCells.insert(cellIdx);
                }
            }

            for (const auto& [cellType, cellIndices] : zoneElemsByType)
            {
                BlockInfo block;
                std::string originalName = zone.name + "-" + cellType;
                block.name = getBlockName(originalName);
                block.cellIndices = cellIndices;
                block.cellType = cellType;
                blocks.push_back(block);
            }
        }

        std::map<std::string, std::vector<int>> unzonedElemsByType;
        for (size_t i = 0; i < cells.size(); ++i)
        {
            if (zonedCells.find(i) == zonedCells.end() && !cellDecomposed[i])
            {
                unzonedElemsByType[cells[i].type].push_back(i);
            }
        }

        for (const auto& [cellType, cellIndices] : unzonedElemsByType)
        {
            BlockInfo block;
            std::string originalName = "unzoned-" + cellType;
            block.name = getBlockName(originalName);
            block.cellIndices = cellIndices;
            block.cellType = cellType;
            blocks.push_back(block);
        }
    }
    else
    {
        std::map<std::string, std::vector<int>> elemsByType;
        for (size_t i = 0; i < cells.size(); ++i)
        {
            if (!cellDecomposed[i])
                elemsByType[cells[i].type].push_back(i);
        }

        for (const auto& [cellType, cellIndices] : elemsByType)
        {
            BlockInfo block;
            std::string originalName = "fluid-" + cellType;
            block.name = getBlockName(originalName);
            block.cellIndices = cellIndices;
            block.cellType = cellType;
            blocks.push_back(block);
        }
    }

    // Append decomposed tet/pyramid blocks grouped by region (zone/mesh), so a
    // rotor poly block never mixes with a stator one.
    {
        std::vector<std::string> grp = buildCellGroups(cells, cellZones);
        std::map<std::string, std::vector<int>> tetByGrp, pyrByGrp;
        for (int g = 0; g < (int)polySubElems.size(); ++g)
        {
            const std::string& k = grp[polySubElemCell[g]];
            (polySubElems[g].type == 'T' ? tetByGrp : pyrByGrp)[k].push_back(g);
        }
        for (auto& [k, gids] : tetByGrp)
        {
            BlockInfo block;
            block.name = getBlockName(k + "-poly-tet");
            block.cellType = "tet";
            block.decomposed = true;
            block.subIndices = std::move(gids);
            blocks.push_back(std::move(block));
        }
        for (auto& [k, gids] : pyrByGrp)
        {
            BlockInfo block;
            block.name = getBlockName(k + "-poly-pyr");
            block.cellType = "pyr";
            block.decomposed = true;
            block.subIndices = std::move(gids);
            blocks.push_back(std::move(block));
        }
    }

    int blockId = 1;
    for (const auto& block : blocks)
    {
        int numElemsInBlock =
            block.decomposed ? (int)block.subIndices.size()
                             : (int)block.cellIndices.size();

        int numNodesPerElem = 8;
        std::string exoType = "HEX8";

        if (block.cellType == "tet")
        {
            numNodesPerElem = 4;
            exoType = "TETRA4";
        }
        else if (block.cellType == "pyr")
        {
            numNodesPerElem = 5;
            exoType = "PYRAMID5";
        }
        else if (block.cellType == "wedge")
        {
            numNodesPerElem = 6;
            exoType = "WEDGE6";
        }
        else if (block.cellType == "hex")
        {
            numNodesPerElem = 8;
            exoType = "HEX8";
        }

        int dim_num_el_in_blk, dim_num_nod_per_el;
        std::string blk_num = std::to_string(blockId);

        int status;
        status = nc_def_dim(ncid,
                            ("num_el_in_blk" + blk_num).c_str(),
                            numElemsInBlock,
                            &dim_num_el_in_blk);
        checkError(status, "Failed to define element block dimension");

        status = nc_def_dim(ncid,
                            ("num_nod_per_el" + blk_num).c_str(),
                            numNodesPerElem,
                            &dim_num_nod_per_el);
        checkError(status, "Failed to define nodes per element dimension");

        int var_connect;
        int dims_connect[2] = {dim_num_el_in_blk, dim_num_nod_per_el};
        status = nc_def_var(ncid,
                            ("connect" + blk_num).c_str(),
                            NC_INT,
                            2,
                            dims_connect,
                            &var_connect);
        checkError(status, "Failed to define connectivity variable");

        status = nc_put_att_text(
            ncid, var_connect, "elem_type", exoType.length(), exoType.c_str());
        checkError(status, "Failed to set element type attribute");

        blockId++;
    }

    nc_enddef(ncid);

    cellToExodusElem.assign(cells.size(), 0);
    int nextExodusElem = 1;

    blockId = 1;
    for (const auto& block : blocks)
    {
        std::vector<int> connectivity;

        if (block.decomposed)
        {
            for (int gid : block.subIndices)
            {
                polySubElemExoId[gid] = nextExodusElem++;
                for (int n : polySubElems[gid].nodes)
                    connectivity.push_back(n + 1);
            }
        }
        else
            for (int cellIdx : block.cellIndices)
            {
                cellToExodusElem[cellIdx] = nextExodusElem++;

                const auto& cell = cells[cellIdx];

                std::vector<int> nodes;
                if (block.cellType == "hex")
                {
                    nodes = orderHexNodes(cell, faces, points);
                }
                else if (block.cellType == "tet")
                {
                    nodes = orderTetNodes(cell, faces, points);
                }
                else if (block.cellType == "pyr")
                {
                    nodes = orderPyramidNodes(cell, faces, points);
                }
                else if (block.cellType == "wedge")
                {
                    nodes = orderWedgeNodes(cell, faces, points);
                }

                int numNodesPerElem = (block.cellType == "tet")     ? 4
                                      : (block.cellType == "pyr")   ? 5
                                      : (block.cellType == "wedge") ? 6
                                                                    : 8;
                if ((int)nodes.size() != numNodesPerElem)
                {
                    std::cerr << "Warning: cell " << cellIdx
                              << " (type=" << block.cellType << ") returned "
                              << nodes.size() << " nodes, expected "
                              << numNodesPerElem << std::endl;
                }
                for (int i = 0; i < numNodesPerElem; ++i)
                {
                    connectivity.push_back(i < (int)nodes.size() ? nodes[i] + 1
                                                                 : 1);
                }
            }

        int var_id;
        nc_inq_varid(
            ncid, ("connect" + std::to_string(blockId)).c_str(), &var_id);
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

        std::cout << "Wrote element block " << blockId << " (" << block.name
                  << ") with "
                  << (block.decomposed ? block.subIndices.size()
                                       : block.cellIndices.size())
                  << " elements"
                  << std::endl;

        blockId++;
    }
}

void ExodusWriter::writeSideSets(const MergedMeshReader& reader)
{
    const auto& patches = reader.getBoundaryPatches();
    const auto& faces = reader.getFaces();
    const auto& owner = reader.getOwner();
    const auto& cells = reader.getCells();
    const auto& points = reader.getPoints();

    if (patches.empty())
    {
        std::cout << "No boundary patches to write as sidesets" << std::endl;
        return;
    }

    // Ordered nodes for standard cells only; polyhedral cells resolve their
    // boundary faces through polyFaceToSubs instead.
    std::map<int, std::vector<int>> cellToOrderedNodes;
    for (size_t i = 0; i < cells.size(); ++i)
    {
        const auto& cell = cells[i];
        if (cell.type == "hex")
            cellToOrderedNodes[i] = orderHexNodes(cell, faces, points);
        else if (cell.type == "tet")
            cellToOrderedNodes[i] = orderTetNodes(cell, faces, points);
        else if (cell.type == "pyr")
            cellToOrderedNodes[i] = orderPyramidNodes(cell, faces, points);
        else if (cell.type == "wedge")
            cellToOrderedNodes[i] = orderWedgeNodes(cell, faces, points);
    }

    // Build each patch's (elem, side) entries first. A boundary face of a
    // polyhedral cell maps to the base side(s) of its sub-element(s); for an
    // n-gon face that is several sides, so the per-sideset count is not known
    // until the faces are walked.
    std::vector<std::vector<int>> patchElem(patches.size());
    std::vector<std::vector<int>> patchSide(patches.size());
    for (size_t i = 0; i < patches.size(); ++i)
    {
        const auto& patch = patches[i];
        auto& elem_list = patchElem[i];
        auto& side_list = patchSide[i];

        for (int j = 0; j < patch.nFaces; ++j)
        {
            int faceIdx = patch.startFace + j;
            if (faceIdx >= (int)owner.size())
                continue;

            auto pit = polyFaceToSubs.find(faceIdx);
            if (pit != polyFaceToSubs.end())
            {
                for (const auto& entry : pit->second)
                {
                    elem_list.push_back(polySubElemExoId[entry.first]);
                    side_list.push_back(entry.second);
                }
                continue;
            }

            int cellIdx = owner[faceIdx];
            int exoId =
                (cellIdx >= 0 && cellIdx < (int)cellToExodusElem.size())
                    ? cellToExodusElem[cellIdx]
                    : (cellIdx + 1);

            const auto& face = faces[faceIdx];
            const auto& cell = cells[cellIdx];
            int sideId = 1;
            auto cit = cellToOrderedNodes.find(cellIdx);
            if (cit != cellToOrderedNodes.end())
            {
                const auto& ordNodes = cit->second;
                if (cell.type == "hex")
                    sideId = getHexFaceId(face.pointIndices, ordNodes);
                else if (cell.type == "tet")
                    sideId = getTetFaceId(face.pointIndices, ordNodes);
                else if (cell.type == "pyr")
                    sideId = getPyramidFaceId(face.pointIndices, ordNodes);
                else if (cell.type == "wedge")
                    sideId = getWedgeFaceId(face.pointIndices, ordNodes);
            }
            elem_list.push_back(exoId);
            side_list.push_back(sideId);
        }
    }

    nc_redef(ncid);

    for (size_t i = 0; i < patches.size(); ++i)
    {
        const auto& patch = patches[i];
        std::string ss_num = std::to_string(i + 1);

        int dim_num_side_ss;
        nc_def_dim(ncid,
                   ("num_side_ss" + ss_num).c_str(),
                   patchElem[i].size(),
                   &dim_num_side_ss);

        int var_elem_ss, var_side_ss;
        nc_def_var(ncid,
                   ("elem_ss" + ss_num).c_str(),
                   NC_INT,
                   1,
                   &dim_num_side_ss,
                   &var_elem_ss);
        nc_def_var(ncid,
                   ("side_ss" + ss_num).c_str(),
                   NC_INT,
                   1,
                   &dim_num_side_ss,
                   &var_side_ss);

        std::string sidesetName = getSidesetName(patch.name);
        nc_put_att_text(ncid,
                        var_elem_ss,
                        "name",
                        sidesetName.length(),
                        sidesetName.c_str());
    }

    int var_ss_names;
    int dim_num_side_sets, dim_len_name;
    int status;
    status = nc_inq_dimid(ncid, "num_side_sets", &dim_num_side_sets);
    checkError(status, "Failed to get num_side_sets dimension");

    status = nc_inq_dimid(ncid, "len_name", &dim_len_name);
    checkError(status, "Failed to get len_name dimension");

    int dims_ss_names[2] = {dim_num_side_sets, dim_len_name};
    status =
        nc_def_var(ncid, "ss_names", NC_CHAR, 2, dims_ss_names, &var_ss_names);
    checkError(status, "Failed to define ss_names variable");

    nc_enddef(ncid);

    for (size_t i = 0; i < patches.size(); ++i)
    {
        const auto& patch = patches[i];
        std::string ss_num = std::to_string(i + 1);

        const std::vector<int>& elem_list = patchElem[i];
        const std::vector<int>& side_list = patchSide[i];

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

        std::string sidesetName = getSidesetName(patch.name);
        char name_padded[33];
        std::memset(name_padded, 0, 33);
        std::strncpy(name_padded, sidesetName.c_str(), 32);

        nc_inq_varid(ncid, "ss_names", &var_id);
        size_t start[2] = {i, 0};
        size_t count[2] = {1, 33};
        nc_put_vara_text(ncid, var_id, start, count, name_padded);

        std::cout << "Wrote sideset " << (i + 1) << ": " << sidesetName
                  << " with " << elem_list.size() << " sides ("
                  << patch.nFaces << " boundary faces)" << std::endl;
    }
}

// Template implementations (kept in header for inlining)
template <typename ReaderType>
void ExodusWriter::writeElementsImpl(const ReaderType& reader)
{
    // This is just a placeholder - actual implementation is in the concrete
    // overloads above
}

template <typename ReaderType>
void ExodusWriter::writeSideSetsImpl(const ReaderType& reader)
{
    // This is just a placeholder - actual implementation is in the concrete
    // overloads above
}
