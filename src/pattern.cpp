//#include "pattern.h"
//#include <iostream>
//
//Pattern::Pattern() {}
//
//Pattern::Pattern(AllCellVertices v, AllCellFaces f, AllCellEdges e):o_cellVertices(v), o_cellFaces(f), o_cellEdges(e)
//{}
//
//std::vector<spCell> Pattern::getCells() {
//    return cells;
//}
//
//AllCellVertices Pattern::getVertices() {
//    return o_cellVertices;
//}
//
//AllCellFaces Pattern::getFaces() {
//    return o_cellFaces;
//}
//
//int Pattern::numCells() {
//    return cells.size();
//}
//
//void Pattern::setVertices(const AllCellVertices& verts) {
//    for (int i = 0; i < verts.size(); i++) {
//        auto& cell = verts[i];
//        for (int j = 0; j < cell.size(); j++) {
//            o_cellVertices[i][j] = cell[j];
//        }
//    }
//}
//
//void Pattern::setFaces(const AllCellFaces& faces) {
//    for (int i = 0; i < faces.size(); i++) {
//        auto& cell = faces[i];
//        for (int j = 0; j < cell.size(); j++) {
//            auto& face = cell[j];
//            for (int k = 0; k < face.size(); k++) {
//                o_cellFaces[i][j][k] = face[k];
//            }
//        }
//    }
//}
//
//void Pattern::createCellsfromVoro() {
//    // for every cell create Cell and push to container
//    for (int i = 0; i < o_cellFaces.size(); ++i) {
//        auto curCell = o_cellFaces[i];
//        Surface_mesh sm; 
//        buildSMfromVF(o_cellVertices[i], o_cellFaces[i], sm);
//        spCell cellptr(new Cell{ i, std::vector<spConvex>{}, o_cellVertices[i], o_cellFaces[i], sm });
//        cells.push_back(cellptr);
//    }
//}
