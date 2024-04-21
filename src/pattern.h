///**
// * Fracture Pattern class.
//**/
//#pragma once
//
//#include <Windows.h>
//#include "clipper.h"
//
//class Pattern
//{
//public:
//    Pattern();
//    Pattern(AllCellVertices, AllCellFaces, AllCellEdges);
//    ~Pattern() {};
//
//    void createCellsfromVoro();
//    std::vector<spCell> getCells();
//    AllCellVertices getVertices();
//    AllCellFaces getFaces();
//    int numCells();
//
//    void setVertices(const AllCellVertices& verts);
//    void setFaces(const AllCellFaces& faces);
//
//private:
//    // original data outputed from Voronoi function of libigl
//    AllCellVertices o_cellVertices;
//    AllCellFaces o_cellFaces;
//    AllCellEdges o_cellEdges;
//    // data converted for mesh clipping and convex hull algorithm
//    std::vector<spCell> cells;
//    // possible data to be added 
//    // 1. the center of the pattern or the bounding box of Pattern 
//    // 2. transformation matrices to move the pattern around
//    // 3. scal
//
//};