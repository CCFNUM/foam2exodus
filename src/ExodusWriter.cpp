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
    // Node sets are not exported yet; the parameter keeps the signature
    // aligned with the Exodus initialisation call.
    (void)numNodeSets;
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

namespace
{

// OpenFOAM stores every face with its right-hand-rule normal pointing out of
// the face's owner cell. Exodus (like OpenFOAM's own cellShape models) wants
// the base face of an element wound the other way round, with its normal
// pointing into the element, so the loop is reversed for the owner and taken
// as stored for the neighbour. This is exact topology: it never looks at the
// geometry, so slivers and warped cells are oriented as reliably as cubes.
std::vector<int> inwardLoop(const std::vector<int>& facePoints,
                            int faceIdx,
                            int cellIdx,
                            const std::vector<int>& owner)
{
    std::vector<int> loop(facePoints);
    if (faceIdx >= 0 && faceIdx < (int)owner.size() &&
        owner[faceIdx] == cellIdx)
    {
        std::reverse(loop.begin(), loop.end());
    }
    return loop;
}

// Number of vertices of each face of a cell, plus how often each vertex is
// used. Returns false if any face index is out of range.
bool cellVertexUse(const Cell& cell,
                   const std::vector<Face>& faces,
                   std::map<int, int>& vertexUse,
                   int& nTri,
                   int& nQuad,
                   int& nOther)
{
    vertexUse.clear();
    nTri = nQuad = nOther = 0;
    for (int fi : cell.faceIndices)
    {
        if (fi < 0 || fi >= (int)faces.size())
            return false;
        const auto& fv = faces[fi].pointIndices;
        if (fv.size() == 3)
            ++nTri;
        else if (fv.size() == 4)
            ++nQuad;
        else
            ++nOther;
        for (int n : fv)
            ++vertexUse[n];
    }
    return true;
}

// Given a base edge (n0,n1) of a cell, return the vertex adjacent to n0 on the
// side face carrying that edge, i.e. the node directly "above" n0. -1 if no
// side face carries the edge, which means the cell is not the assumed shape.
int nodeAbove(int n0,
              int n1,
              int baseFace,
              const Cell& cell,
              const std::vector<Face>& faces,
              const std::set<int>& baseSet)
{
    for (int fi : cell.faceIndices)
    {
        if (fi == baseFace)
            continue;
        const auto& sf = faces[fi].pointIndices;
        const int nv = (int)sf.size();
        int a = -1, b = -1;
        for (int j = 0; j < nv; ++j)
        {
            if (sf[j] == n0)
                a = j;
            if (sf[j] == n1)
                b = j;
        }
        if (a < 0 || b < 0)
            continue;
        int cand = -1;
        if ((a + 1) % nv == b)
            cand = sf[(a + nv - 1) % nv];
        else if ((a + nv - 1) % nv == b)
            cand = sf[(a + 1) % nv];
        if (cand >= 0 && baseSet.count(cand) == 0)
            return cand;
    }
    return -1;
}

} // namespace

std::vector<int> ExodusWriter::orderTetNodes(int cellIdx,
                                             const Cell& cell,
                                             const std::vector<Face>& faces,
                                             const std::vector<int>& owner)
{
    // tetMatcher equivalent: four triangles, four vertices, each on 3 faces.
    if (cell.faceIndices.size() != 4)
        return {};
    std::map<int, int> vertexUse;
    int nTri = 0, nQuad = 0, nOther = 0;
    if (!cellVertexUse(cell, faces, vertexUse, nTri, nQuad, nOther))
        return {};
    if (nTri != 4 || nQuad != 0 || nOther != 0 || vertexUse.size() != 4)
        return {};
    for (const auto& vu : vertexUse)
        if (vu.second != 3)
            return {};

    const int baseFace = cell.faceIndices[0];
    std::vector<int> nodes =
        inwardLoop(faces[baseFace].pointIndices, baseFace, cellIdx, owner);
    int apex = -1;
    for (const auto& vu : vertexUse)
        if (std::find(nodes.begin(), nodes.end(), vu.first) == nodes.end())
            apex = vu.first;
    if (apex < 0)
        return {};
    nodes.push_back(apex);
    return nodes;
}

std::vector<int> ExodusWriter::orderPyramidNodes(int cellIdx,
                                                 const Cell& cell,
                                                 const std::vector<Face>& faces,
                                                 const std::vector<int>& owner)
{
    // pyrMatcher equivalent: one quad base, four triangles, five vertices with
    // the apex on four faces and each base vertex on three.
    if (cell.faceIndices.size() != 5)
        return {};
    std::map<int, int> vertexUse;
    int nTri = 0, nQuad = 0, nOther = 0;
    if (!cellVertexUse(cell, faces, vertexUse, nTri, nQuad, nOther))
        return {};
    if (nTri != 4 || nQuad != 1 || nOther != 0 || vertexUse.size() != 5)
        return {};

    int quadFace = -1;
    for (int fi : cell.faceIndices)
        if (faces[fi].pointIndices.size() == 4)
            quadFace = fi;

    std::vector<int> nodes =
        inwardLoop(faces[quadFace].pointIndices, quadFace, cellIdx, owner);
    int apex = -1;
    for (const auto& vu : vertexUse)
        if (std::find(nodes.begin(), nodes.end(), vu.first) == nodes.end())
            apex = vu.first;
    if (apex < 0 || vertexUse[apex] != 4)
        return {};
    nodes.push_back(apex);
    return nodes;
}

std::vector<int> ExodusWriter::orderWedgeNodes(int cellIdx,
                                               const Cell& cell,
                                               const std::vector<Face>& faces,
                                               const std::vector<int>& owner)
{
    // prismMatcher equivalent: two triangles, three quads, six vertices each
    // shared by three faces.
    if (cell.faceIndices.size() != 5)
        return {};
    std::map<int, int> vertexUse;
    int nTri = 0, nQuad = 0, nOther = 0;
    if (!cellVertexUse(cell, faces, vertexUse, nTri, nQuad, nOther))
        return {};
    if (nTri != 2 || nQuad != 3 || nOther != 0 || vertexUse.size() != 6)
        return {};
    for (const auto& vu : vertexUse)
        if (vu.second != 3)
            return {};

    int baseFace = -1;
    for (int fi : cell.faceIndices)
        if (faces[fi].pointIndices.size() == 3)
        {
            baseFace = fi;
            break;
        }

    std::vector<int> base =
        inwardLoop(faces[baseFace].pointIndices, baseFace, cellIdx, owner);
    std::set<int> baseSet(base.begin(), base.end());

    std::vector<int> nodes(6, -1);
    for (int i = 0; i < 3; ++i)
    {
        nodes[i] = base[i];
        int top = nodeAbove(
            base[i], base[(i + 1) % 3], baseFace, cell, faces, baseSet);
        if (top < 0)
            return {};
        nodes[3 + i] = top;
    }
    if (std::set<int>(nodes.begin(), nodes.end()).size() != 6)
        return {};
    return nodes;
}

std::vector<int> ExodusWriter::orderHexNodes(int cellIdx,
                                             const Cell& cell,
                                             const std::vector<Face>& faces,
                                             const std::vector<int>& owner)
{
    // hexMatcher equivalent: six quadrilaterals, eight vertices, every vertex
    // shared by exactly three faces. That set of conditions admits only the
    // hexahedron, so no geometric test is needed to recognise one.
    if (cell.faceIndices.size() != 6)
        return {};
    std::map<int, int> vertexUse;
    int nTri = 0, nQuad = 0, nOther = 0;
    if (!cellVertexUse(cell, faces, vertexUse, nTri, nQuad, nOther))
        return {};
    if (nQuad != 6 || nTri != 0 || nOther != 0 || vertexUse.size() != 8)
        return {};
    for (const auto& vu : vertexUse)
        if (vu.second != 3)
            return {};

    // Any face can serve as the base; the opposite face must share no vertex
    // with it, which is the remaining condition that rules out non-hexahedra.
    const int baseFace = cell.faceIndices[0];
    std::vector<int> base =
        inwardLoop(faces[baseFace].pointIndices, baseFace, cellIdx, owner);
    std::set<int> baseSet(base.begin(), base.end());

    bool hasOpposite = false;
    for (int fi : cell.faceIndices)
    {
        if (fi == baseFace)
            continue;
        bool shares = false;
        for (int n : faces[fi].pointIndices)
            if (baseSet.count(n))
            {
                shares = true;
                break;
            }
        if (!shares)
        {
            hasOpposite = true;
            break;
        }
    }
    if (!hasOpposite)
        return {};

    std::vector<int> nodes(8, -1);
    for (int i = 0; i < 4; ++i)
    {
        nodes[i] = base[i];
        int top = nodeAbove(
            base[i], base[(i + 1) % 4], baseFace, cell, faces, baseSet);
        if (top < 0)
            return {};
        nodes[4 + i] = top;
    }
    if (std::set<int>(nodes.begin(), nodes.end()).size() != 8)
        return {};
    return nodes;
}

std::vector<int>
ExodusWriter::orderStandardNodes(int cellIdx,
                                 const Cell& cell,
                                 const std::vector<Face>& faces,
                                 const std::vector<int>& owner)
{
    if (cell.type == "hex")
        return orderHexNodes(cellIdx, cell, faces, owner);
    if (cell.type == "tet")
        return orderTetNodes(cellIdx, cell, faces, owner);
    if (cell.type == "pyr")
        return orderPyramidNodes(cellIdx, cell, faces, owner);
    if (cell.type == "wedge")
        return orderWedgeNodes(cellIdx, cell, faces, owner);
    return {};
}

std::vector<std::string>
ExodusWriter::buildCellGroups(const std::vector<Cell>& cells,
                              const std::vector<CellZone>& cellZones) const
{
    // Mirror the standard-block grouping in writeElements: with >1 zone, cells
    // are grouped by zone name (or "unzoned"); otherwise a single "fluid"
    // group.
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

namespace
{

// Local topology tables in Exodus node order, shared by the element checker.
// Corner triads are the edge triples whose determinant is the Jacobian at that
// vertex (the Verdict/VTK definition, so signs match ParaView's Mesh Quality).
struct ElemTopo
{
    int nNodes;
    std::vector<std::array<int, 2>> edges;
    std::vector<std::array<int, 3>> corners;
    // Outward faces; a -1 in the fourth slot marks a triangle.
    std::vector<std::array<int, 4>> faces;
};

const ElemTopo& elemTopo(const std::string& type)
{
    static const ElemTopo hex{8,
                              {{0, 1},
                               {1, 2},
                               {2, 3},
                               {3, 0},
                               {4, 5},
                               {5, 6},
                               {6, 7},
                               {7, 4},
                               {0, 4},
                               {1, 5},
                               {2, 6},
                               {3, 7}},
                              {{1, 3, 4},
                               {2, 0, 5},
                               {3, 1, 6},
                               {0, 2, 7},
                               {7, 5, 0},
                               {4, 6, 1},
                               {5, 7, 2},
                               {6, 4, 3}},
                              {{0, 1, 5, 4},
                               {1, 2, 6, 5},
                               {2, 3, 7, 6},
                               {3, 0, 4, 7},
                               {0, 3, 2, 1},
                               {4, 5, 6, 7}}};
    static const ElemTopo tet{
        4,
        {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}},
        {{1, 2, 3}, {2, 0, 3}, {0, 1, 3}, {2, 1, 0}},
        {{0, 1, 3, -1}, {1, 2, 3, -1}, {2, 0, 3, -1}, {0, 2, 1, -1}}};
    static const ElemTopo pyr{
        5,
        {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 4}, {1, 4}, {2, 4}, {3, 4}},
        {{1, 3, 4}, {2, 0, 4}, {3, 1, 4}, {0, 2, 4}},
        {{0, 1, 4, -1},
         {1, 2, 4, -1},
         {2, 3, 4, -1},
         {3, 0, 4, -1},
         {0, 3, 2, 1}}};
    static const ElemTopo wedge{
        6,
        {{0, 1},
         {1, 2},
         {2, 0},
         {3, 4},
         {4, 5},
         {5, 3},
         {0, 3},
         {1, 4},
         {2, 5}},
        {{1, 2, 3}, {2, 0, 4}, {0, 1, 5}, {5, 4, 0}, {3, 5, 1}, {4, 3, 2}},
        {{0, 1, 4, 3},
         {1, 2, 5, 4},
         {2, 0, 3, 5},
         {0, 2, 1, -1},
         {3, 4, 5, -1}}};
    static const ElemTopo none{0, {}, {}, {}};
    if (type == "hex")
        return hex;
    if (type == "tet")
        return tet;
    if (type == "pyr")
        return pyr;
    if (type == "wedge")
        return wedge;
    return none;
}

double det3(const double a[3], const double b[3], const double c[3])
{
    return a[0] * (b[1] * c[2] - b[2] * c[1]) +
           a[1] * (b[2] * c[0] - b[0] * c[2]) +
           a[2] * (b[0] * c[1] - b[1] * c[0]);
}

// Trilinear volume of a HEX8 by 2x2x2 Gauss quadrature. The Jacobian
// determinant of a trilinear map is quadratic in each coordinate, so two
// points per direction integrate it exactly; this reproduces the volume
// ParaView reports. A volume built by fanning each quad face from one of its
// corners does not, and for warped high-aspect cells it can even change sign.
double hexTrilinearVolume(const std::vector<const Point*>& p)
{
    static const int sgn[8][3] = {{-1, -1, -1},
                                  {1, -1, -1},
                                  {1, 1, -1},
                                  {-1, 1, -1},
                                  {-1, -1, 1},
                                  {1, -1, 1},
                                  {1, 1, 1},
                                  {-1, 1, 1}};
    const double g = 1.0 / std::sqrt(3.0);
    double vol = 0.0;
    for (int q = 0; q < 8; ++q)
    {
        const double xi = sgn[q][0] * g, eta = sgn[q][1] * g,
                     ze = sgn[q][2] * g;
        double J[3][3] = {{0}};
        for (int i = 0; i < 8; ++i)
        {
            const double a = sgn[i][0], b = sgn[i][1], c = sgn[i][2];
            const double d[3] = {a * (1 + b * eta) * (1 + c * ze) / 8.0,
                                 b * (1 + a * xi) * (1 + c * ze) / 8.0,
                                 c * (1 + a * xi) * (1 + b * eta) / 8.0};
            const double x[3] = {p[i]->x, p[i]->y, p[i]->z};
            for (int r = 0; r < 3; ++r)
                for (int s = 0; s < 3; ++s)
                    J[r][s] += d[r] * x[s];
        }
        vol += det3(J[0], J[1], J[2]);
    }
    return vol;
}

// Volume from the outward faces via the divergence theorem, with quads split
// about their centroid so a warped face contributes symmetrically. Used for
// tets, pyramids and wedges, and matches how OpenFOAM computes cell volumes.
double faceFanVolume(const std::vector<const Point*>& p, const ElemTopo& topo)
{
    auto triDiv = [](const Point& A, const Point& B, const Point& C)
    {
        const double a[3] = {A.x, A.y, A.z};
        const double b[3] = {B.x, B.y, B.z};
        const double c[3] = {C.x, C.y, C.z};
        return det3(a, b, c);
    };
    double vol = 0.0;
    for (const auto& f : topo.faces)
    {
        if (f[3] < 0)
        {
            vol += triDiv(*p[f[0]], *p[f[1]], *p[f[2]]);
            continue;
        }
        Point ctr{0, 0, 0};
        for (int k = 0; k < 4; ++k)
        {
            ctr.x += p[f[k]]->x / 4.0;
            ctr.y += p[f[k]]->y / 4.0;
            ctr.z += p[f[k]]->z / 4.0;
        }
        for (int k = 0; k < 4; ++k)
            vol += triDiv(*p[f[k]], *p[f[(k + 1) % 4]], ctr);
    }
    return vol / 6.0;
}

} // namespace

ExodusWriter::ElemCheck
ExodusWriter::checkElement(const std::string& type,
                           const std::vector<int>& nodes,
                           const std::vector<Point>& points) const
{
    ElemCheck res;
    const ElemTopo& topo = elemTopo(type);
    if (topo.nNodes == 0)
    {
        res.reason = "unsupported topology '" + type + "'";
        return res;
    }
    if ((int)nodes.size() != topo.nNodes)
    {
        res.reason = "expected " + std::to_string(topo.nNodes) +
                     " nodes, got " + std::to_string(nodes.size());
        return res;
    }

    const int nBase = (int)points.size();
    std::vector<const Point*> p(topo.nNodes);
    for (int i = 0; i < topo.nNodes; ++i)
    {
        const int n = nodes[i];
        if (n < 0 || n >= nBase + (int)polyExtraPoints.size())
        {
            res.reason = "node index " + std::to_string(n) + " out of range";
            return res;
        }
        p[i] = (n < nBase) ? &points[n] : &polyExtraPoints[n - nBase];
    }

    if (std::set<int>(nodes.begin(), nodes.end()).size() != (size_t)topo.nNodes)
    {
        res.reason = "repeated node ID in connectivity";
        return res;
    }

    // Longest edge sets the length scale for the degeneracy tolerances.
    double maxEdge = 0.0;
    res.minEdge = 1e300;
    for (const auto& e : topo.edges)
    {
        const double dx = p[e[1]]->x - p[e[0]]->x;
        const double dy = p[e[1]]->y - p[e[0]]->y;
        const double dz = p[e[1]]->z - p[e[0]]->z;
        const double len = std::sqrt(dx * dx + dy * dy + dz * dz);
        res.minEdge = std::min(res.minEdge, len);
        maxEdge = std::max(maxEdge, len);
    }
    if (maxEdge <= 0.0)
    {
        res.reason = "all vertices coincide";
        return res;
    }
    for (int i = 0; i < topo.nNodes; ++i)
        for (int j = i + 1; j < topo.nNodes; ++j)
        {
            const double dx = p[j]->x - p[i]->x;
            const double dy = p[j]->y - p[i]->y;
            const double dz = p[j]->z - p[i]->z;
            if (std::sqrt(dx * dx + dy * dy + dz * dz) <= 1e-12 * maxEdge)
            {
                res.reason = "coincident local vertices " + std::to_string(i) +
                             " and " + std::to_string(j);
                return res;
            }
        }
    if (res.minEdge <= 1e-12 * maxEdge)
    {
        res.reason = "collapsed edge";
        return res;
    }

    // Jacobian at every corner; the minimum normalised value is the scaled
    // Jacobian, negative for an inverted or wrongly ordered element.
    res.scaledJacobian = 1e300;
    for (size_t ci = 0; ci < topo.corners.size(); ++ci)
    {
        // Triad entry k is the far end of the k-th edge leaving corner ci.
        const auto& c = topo.corners[ci];
        double v[3][3];
        double len[3];
        for (int k = 0; k < 3; ++k)
        {
            const Point& q = *p[c[k]];
            v[k][0] = q.x - p[ci]->x;
            v[k][1] = q.y - p[ci]->y;
            v[k][2] = q.z - p[ci]->z;
            len[k] = std::sqrt(v[k][0] * v[k][0] + v[k][1] * v[k][1] +
                               v[k][2] * v[k][2]);
        }
        if (len[0] <= 0.0 || len[1] <= 0.0 || len[2] <= 0.0)
        {
            res.reason = "degenerate corner";
            return res;
        }
        double sj = det3(v[0], v[1], v[2]) / (len[0] * len[1] * len[2]);
        if (type == "tet")
            sj *= std::sqrt(2.0);
        res.scaledJacobian = std::min(res.scaledJacobian, sj);
    }

    res.volume =
        (type == "hex") ? hexTrilinearVolume(p) : faceFanVolume(p, topo);

    if (res.scaledJacobian <= 0.0)
    {
        res.reason = "non-positive Jacobian (scaled Jacobian " +
                     std::to_string(res.scaledJacobian) + ")";
        return res;
    }
    if (res.volume <= 0.0)
    {
        res.reason = "non-positive volume (" + std::to_string(res.volume) + ")";
        return res;
    }

    res.valid = true;
    return res;
}

void ExodusWriter::validateElements(const std::vector<Point>& points,
                                    const std::vector<Cell>& cells) const
{
    std::vector<std::string> errors;
    double minSJ = 1e300, maxSJ = -1e300, minVol = 1e300;
    int nChecked = 0;

    auto record = [&](const ElemCheck& chk, const std::string& what)
    {
        if (!chk.valid)
        {
            if (errors.size() < 20)
                errors.push_back(what + ": " + chk.reason);
            return;
        }
        minSJ = std::min(minSJ, chk.scaledJacobian);
        maxSJ = std::max(maxSJ, chk.scaledJacobian);
        minVol = std::min(minVol, chk.volume);
        ++nChecked;
    };

    for (size_t c = 0; c < cells.size(); ++c)
    {
        if (cellDecomposed[c] || cellOrderedNodes[c].empty())
            continue;
        record(checkElement(cells[c].type, cellOrderedNodes[c], points),
               "OpenFOAM cell " + std::to_string(c) + " (" + cells[c].type +
                   ")");
    }

    for (size_t g = 0; g < polySubElems.size(); ++g)
    {
        const std::string type = (polySubElems[g].type == 'T') ? "tet" : "pyr";
        record(checkElement(type, polySubElems[g].nodes, points),
               "sub-element " + std::to_string(polySubElems[g].srcSub) +
                   " of OpenFOAM face " +
                   std::to_string(polySubElems[g].srcFace) +
                   " of OpenFOAM cell " + std::to_string(polySubElemCell[g]) +
                   " (" + type + ")");
    }

    if (!errors.empty())
    {
        std::string msg = "mesh validation failed for " +
                          std::to_string(errors.size()) +
                          "+ element(s); no Exodus file written:";
        for (const auto& e : errors)
            msg += "\n  " + e;
        throw std::runtime_error(msg);
    }

    std::cout << "Validated " << nChecked << " elements: scaled Jacobian in ["
              << minSJ << ", " << maxSJ << "], min volume " << minVol
              << std::endl;
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
    { return (idx < nBase) ? points[idx] : polyExtraPoints[idx - nBase]; };

    // Face centroids are shared per face so both adjacent poly cells fan an
    // n-gon identically, keeping the split conformal.
    std::unordered_map<int, int> faceCentroidNode;

    cellOrderedNodes.assign(cells.size(), {});
    int nRejected = 0;

    for (int c = 0; c < (int)cells.size(); ++c)
    {
        const auto& cell = cells[c];

        // A cell is written as a standard Exodus element when it matches that
        // topology and the resulting element is geometrically sound. Anything
        // else - a genuine polyhedron, or a standard cell whose geometry is
        // unusable - is reported and split into conformal tets/pyramids.
        if (cell.type != "unknown")
        {
            std::vector<int> on = orderStandardNodes(c, cell, faces, owner);
            if (!on.empty())
            {
                ElemCheck chk = checkElement(cell.type, on, points);
                if (chk.valid)
                {
                    cellOrderedNodes[c] = std::move(on);
                    continue;
                }
                std::cerr << "Warning: OpenFOAM cell " << c << " (" << cell.type
                          << ") rejected as a standard element: " << chk.reason
                          << "; decomposing it instead" << std::endl;
                ++nRejected;
            }
            else if (cell.type == "hex" || cell.type == "tet" ||
                     cell.type == "pyr" || cell.type == "wedge")
            {
                std::cerr << "Warning: OpenFOAM cell " << c
                          << " has the face counts of a " << cell.type
                          << " but does not match that topology; decomposing it"
                          << std::endl;
                ++nRejected;
            }
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
            // Wind the face so its normal points into this cell, which is the
            // side the centroid apex sits on. Doing it from the owner /
            // neighbour relation rather than from a signed volume keeps the
            // sub-elements correctly oriented even for warped faces.
            const std::vector<int> fv =
                inwardLoop(faces[fi].pointIndices, fi, c, owner);
            const int nv = (int)fv.size();
            const bool onBoundary = (fi >= boundaryStart) && (owner[fi] == c);

            if (nv == 3)
            {
                int gid = (int)polySubElems.size();
                polySubElems.push_back(
                    {'T', {fv[0], fv[1], fv[2], cNode}, fi, 0});
                if (onBoundary)
                    polyFaceToSubs[fi].push_back(
                        {gid, getTetFaceId(fv, polySubElems[gid].nodes)});
            }
            else if (nv == 4)
            {
                int gid = (int)polySubElems.size();
                polySubElems.push_back(
                    {'P', {fv[0], fv[1], fv[2], fv[3], cNode}, fi, 0});
                if (onBoundary)
                    polyFaceToSubs[fi].push_back(
                        {gid, getPyramidFaceId(fv, polySubElems[gid].nodes)});
            }
            else if (nv >= 5)
            {
                // An arbitrary polygon is fanned into nv triangles around a
                // per-face centroid node shared with the cell on the other
                // side of the face, so the split stays conformal. Each
                // triangle plus the cell centroid forms one tetrahedron.
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
                    const int a = fv[k], b = fv[(k + 1) % nv];
                    int gid = (int)polySubElems.size();
                    polySubElems.push_back({'T', {a, b, fNode, cNode}, fi, k});
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

    if (nRejected > 0)
    {
        std::cerr << "Warning: " << nRejected
                  << " cell(s) with standard face counts could not be written "
                     "as standard elements and were decomposed"
                  << std::endl;
    }
}

// Names of the per-element attributes that trace an Exodus element back to the
// OpenFOAM cell it came from. One-to-one elements carry the source cell ID
// only; sub-elements of a decomposed cell also carry the face they were built
// on and their index within that face.
const std::vector<std::string>& ExodusWriter::standardAttribNames()
{
    static const std::vector<std::string> names{"source_openfoam_cell_id"};
    return names;
}

const std::vector<std::string>& ExodusWriter::decomposedAttribNames()
{
    static const std::vector<std::string> names{"source_openfoam_cell_id",
                                                "source_openfoam_face_id",
                                                "source_sub_element_index"};
    return names;
}

void ExodusWriter::writeBlockAttributes(int blockId,
                                        bool decomposed,
                                        const std::vector<double>& values)
{
    const std::string blk_num = std::to_string(blockId);
    const std::vector<std::string>& names =
        decomposed ? decomposedAttribNames() : standardAttribNames();

    int var_id;
    checkError(nc_inq_varid(ncid, ("attrib" + blk_num).c_str(), &var_id),
               "Failed to get element attribute variable");
    if (!values.empty())
    {
        checkError(nc_put_var_double(ncid, var_id, values.data()),
                   "Failed to write element attributes");
    }

    checkError(nc_inq_varid(ncid, ("attrib_name" + blk_num).c_str(), &var_id),
               "Failed to get element attribute name variable");
    for (size_t i = 0; i < names.size(); ++i)
    {
        char buffer[33];
        std::memset(buffer, 0, sizeof(buffer));
        std::strncpy(buffer, names[i].c_str(), 32);
        size_t start[2] = {i, 0};
        size_t count[2] = {1, 33};
        checkError(nc_put_vara_text(ncid, var_id, start, count, buffer),
                   "Failed to write element attribute name");
    }
}

void ExodusWriter::writeElements(const OpenFOAMMeshReader& reader)
{
    const auto& cells = reader.getCells();
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
                if (cellIdx < (int)cells.size())
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

    int dim_len_name;
    checkError(nc_inq_dimid(ncid, "len_name", &dim_len_name),
               "Failed to get len_name dimension");

    int blockId = 1;
    for (const auto& block : blocks)
    {
        int numElemsInBlock = block.decomposed ? (int)block.subIndices.size()
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

        // Per-element attributes carrying the source-cell provenance. They are
        // static Exodus attributes rather than transient element variables, so
        // a pure mesh file stays free of time steps.
        const std::vector<std::string>& attNames =
            block.decomposed ? decomposedAttribNames() : standardAttribNames();
        int dim_num_att;
        status = nc_def_dim(ncid,
                            ("num_att_in_blk" + blk_num).c_str(),
                            attNames.size(),
                            &dim_num_att);
        checkError(status, "Failed to define element attribute dimension");

        int var_attrib, var_attrib_name;
        int dims_attrib[2] = {dim_num_el_in_blk, dim_num_att};
        status = nc_def_var(ncid,
                            ("attrib" + blk_num).c_str(),
                            NC_DOUBLE,
                            2,
                            dims_attrib,
                            &var_attrib);
        checkError(status, "Failed to define element attribute variable");

        int dims_attrib_name[2] = {dim_num_att, dim_len_name};
        status = nc_def_var(ncid,
                            ("attrib_name" + blk_num).c_str(),
                            NC_CHAR,
                            2,
                            dims_attrib_name,
                            &var_attrib_name);
        checkError(status, "Failed to define element attribute name variable");

        blockId++;
    }

    nc_enddef(ncid);

    cellToExodusElem.assign(cells.size(), 0);
    int nextExodusElem = 1;

    blockId = 1;
    for (const auto& block : blocks)
    {
        std::vector<int> connectivity;
        std::vector<double> attribSource;

        if (block.decomposed)
        {
            for (int gid : block.subIndices)
            {
                polySubElemExoId[gid] = nextExodusElem++;
                for (int n : polySubElems[gid].nodes)
                    connectivity.push_back(n + 1);
                attribSource.push_back(polySubElemCell[gid]);
                attribSource.push_back(polySubElems[gid].srcFace);
                attribSource.push_back(polySubElems[gid].srcSub);
            }
        }
        else
            for (int cellIdx : block.cellIndices)
            {
                cellToExodusElem[cellIdx] = nextExodusElem++;
                attribSource.push_back(cellIdx);

                // Connectivity was matched and validated up front, so it is
                // reused verbatim here and in the side sets.
                const std::vector<int>& nodes = cellOrderedNodes[cellIdx];

                int numNodesPerElem = (block.cellType == "tet")     ? 4
                                      : (block.cellType == "pyr")   ? 5
                                      : (block.cellType == "wedge") ? 6
                                                                    : 8;
                if ((int)nodes.size() != numNodesPerElem)
                {
                    throw std::runtime_error(
                        "internal error: OpenFOAM cell " +
                        std::to_string(cellIdx) + " (type=" + block.cellType +
                        ") has " + std::to_string(nodes.size()) +
                        " ordered nodes, expected " +
                        std::to_string(numNodesPerElem));
                }
                for (int i = 0; i < numNodesPerElem; ++i)
                {
                    connectivity.push_back(nodes[i] + 1);
                }
            }

        int var_id;
        nc_inq_varid(
            ncid, ("connect" + std::to_string(blockId)).c_str(), &var_id);
        nc_put_var_int(ncid, var_id, connectivity.data());

        writeBlockAttributes(blockId, block.decomposed, attribSource);

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
                  << " elements" << std::endl;

        blockId++;
    }
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

    if (patches.empty())
    {
        std::cout << "No boundary patches to write as sidesets" << std::endl;
        return;
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
            if (cellIdx < 0 || cellIdx >= (int)cells.size())
                continue;
            int exoId = cellToExodusElem[cellIdx];

            const auto& face = faces[faceIdx];
            const auto& cell = cells[cellIdx];
            int sideId = 1;
            const std::vector<int>& ordNodes = cellOrderedNodes[cellIdx];
            if (!ordNodes.empty())
            {
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

        // A patch with no faces gets no dimension or variables at all: the
        // classic netCDF model has no zero-length fixed dimension. It is still
        // named and counted, but its status is set to 0 below.
        if (patchElem[i].empty())
            continue;

        int dim_num_side_ss;
        checkError(nc_def_dim(ncid,
                              ("num_side_ss" + ss_num).c_str(),
                              patchElem[i].size(),
                              &dim_num_side_ss),
                   "Failed to define sideset dimension");

        int var_elem_ss, var_side_ss;
        checkError(nc_def_var(ncid,
                              ("elem_ss" + ss_num).c_str(),
                              NC_INT,
                              1,
                              &dim_num_side_ss,
                              &var_elem_ss),
                   "Failed to define sideset element variable");
        checkError(nc_def_var(ncid,
                              ("side_ss" + ss_num).c_str(),
                              NC_INT,
                              1,
                              &dim_num_side_ss,
                              &var_side_ss),
                   "Failed to define sideset side variable");

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
        if (!elem_list.empty())
        {
            nc_inq_varid(ncid, ("elem_ss" + ss_num).c_str(), &var_id);
            nc_put_var_int(ncid, var_id, elem_list.data());

            nc_inq_varid(ncid, ("side_ss" + ss_num).c_str(), &var_id);
            nc_put_var_int(ncid, var_id, side_list.data());
        }

        int ss_status = elem_list.empty() ? 0 : 1;
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
                  << " with " << elem_list.size() << " sides (" << patch.nFaces
                  << " boundary faces)" << std::endl;
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

    // Nothing is written until every element that would go into the file has
    // been checked, so a bad mesh fails loudly instead of producing a file
    // with inverted or collapsed elements in it.
    validateElements(reader.getPoints(), cells);

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
                if (cellIdx < (int)cells.size())
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

    // Nothing is written until every element that would go into the file has
    // been checked, so a bad mesh fails loudly instead of producing a file
    // with inverted or collapsed elements in it.
    validateElements(reader.getPoints(), cells);

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
                if (cellIdx < (int)cells.size())
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
                if (cellIdx < (int)cells.size())
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

    int dim_len_name;
    checkError(nc_inq_dimid(ncid, "len_name", &dim_len_name),
               "Failed to get len_name dimension");

    int blockId = 1;
    for (const auto& block : blocks)
    {
        int numElemsInBlock = block.decomposed ? (int)block.subIndices.size()
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

        // Per-element attributes carrying the source-cell provenance. They are
        // static Exodus attributes rather than transient element variables, so
        // a pure mesh file stays free of time steps.
        const std::vector<std::string>& attNames =
            block.decomposed ? decomposedAttribNames() : standardAttribNames();
        int dim_num_att;
        status = nc_def_dim(ncid,
                            ("num_att_in_blk" + blk_num).c_str(),
                            attNames.size(),
                            &dim_num_att);
        checkError(status, "Failed to define element attribute dimension");

        int var_attrib, var_attrib_name;
        int dims_attrib[2] = {dim_num_el_in_blk, dim_num_att};
        status = nc_def_var(ncid,
                            ("attrib" + blk_num).c_str(),
                            NC_DOUBLE,
                            2,
                            dims_attrib,
                            &var_attrib);
        checkError(status, "Failed to define element attribute variable");

        int dims_attrib_name[2] = {dim_num_att, dim_len_name};
        status = nc_def_var(ncid,
                            ("attrib_name" + blk_num).c_str(),
                            NC_CHAR,
                            2,
                            dims_attrib_name,
                            &var_attrib_name);
        checkError(status, "Failed to define element attribute name variable");

        blockId++;
    }

    nc_enddef(ncid);

    cellToExodusElem.assign(cells.size(), 0);
    int nextExodusElem = 1;

    blockId = 1;
    for (const auto& block : blocks)
    {
        std::vector<int> connectivity;
        std::vector<double> attribSource;

        if (block.decomposed)
        {
            for (int gid : block.subIndices)
            {
                polySubElemExoId[gid] = nextExodusElem++;
                for (int n : polySubElems[gid].nodes)
                    connectivity.push_back(n + 1);
                attribSource.push_back(polySubElemCell[gid]);
                attribSource.push_back(polySubElems[gid].srcFace);
                attribSource.push_back(polySubElems[gid].srcSub);
            }
        }
        else
            for (int cellIdx : block.cellIndices)
            {
                cellToExodusElem[cellIdx] = nextExodusElem++;
                attribSource.push_back(cellIdx);

                // Connectivity was matched and validated up front, so it is
                // reused verbatim here and in the side sets.
                const std::vector<int>& nodes = cellOrderedNodes[cellIdx];

                int numNodesPerElem = (block.cellType == "tet")     ? 4
                                      : (block.cellType == "pyr")   ? 5
                                      : (block.cellType == "wedge") ? 6
                                                                    : 8;
                if ((int)nodes.size() != numNodesPerElem)
                {
                    throw std::runtime_error(
                        "internal error: OpenFOAM cell " +
                        std::to_string(cellIdx) + " (type=" + block.cellType +
                        ") has " + std::to_string(nodes.size()) +
                        " ordered nodes, expected " +
                        std::to_string(numNodesPerElem));
                }
                for (int i = 0; i < numNodesPerElem; ++i)
                {
                    connectivity.push_back(nodes[i] + 1);
                }
            }

        int var_id;
        nc_inq_varid(
            ncid, ("connect" + std::to_string(blockId)).c_str(), &var_id);
        nc_put_var_int(ncid, var_id, connectivity.data());

        writeBlockAttributes(blockId, block.decomposed, attribSource);

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
                  << " elements" << std::endl;

        blockId++;
    }
}

void ExodusWriter::writeSideSets(const MergedMeshReader& reader)
{
    const auto& patches = reader.getBoundaryPatches();
    const auto& faces = reader.getFaces();
    const auto& owner = reader.getOwner();
    const auto& cells = reader.getCells();

    if (patches.empty())
    {
        std::cout << "No boundary patches to write as sidesets" << std::endl;
        return;
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
            if (cellIdx < 0 || cellIdx >= (int)cells.size())
                continue;
            int exoId = cellToExodusElem[cellIdx];

            const auto& face = faces[faceIdx];
            const auto& cell = cells[cellIdx];
            int sideId = 1;
            const std::vector<int>& ordNodes = cellOrderedNodes[cellIdx];
            if (!ordNodes.empty())
            {
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

        // A patch with no faces gets no dimension or variables at all: the
        // classic netCDF model has no zero-length fixed dimension. It is still
        // named and counted, but its status is set to 0 below.
        if (patchElem[i].empty())
            continue;

        int dim_num_side_ss;
        checkError(nc_def_dim(ncid,
                              ("num_side_ss" + ss_num).c_str(),
                              patchElem[i].size(),
                              &dim_num_side_ss),
                   "Failed to define sideset dimension");

        int var_elem_ss, var_side_ss;
        checkError(nc_def_var(ncid,
                              ("elem_ss" + ss_num).c_str(),
                              NC_INT,
                              1,
                              &dim_num_side_ss,
                              &var_elem_ss),
                   "Failed to define sideset element variable");
        checkError(nc_def_var(ncid,
                              ("side_ss" + ss_num).c_str(),
                              NC_INT,
                              1,
                              &dim_num_side_ss,
                              &var_side_ss),
                   "Failed to define sideset side variable");

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
        if (!elem_list.empty())
        {
            nc_inq_varid(ncid, ("elem_ss" + ss_num).c_str(), &var_id);
            nc_put_var_int(ncid, var_id, elem_list.data());

            nc_inq_varid(ncid, ("side_ss" + ss_num).c_str(), &var_id);
            nc_put_var_int(ncid, var_id, side_list.data());
        }

        int ss_status = elem_list.empty() ? 0 : 1;
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
                  << " with " << elem_list.size() << " sides (" << patch.nFaces
                  << " boundary faces)" << std::endl;
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
