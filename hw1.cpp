#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <list>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/number_utils.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

struct Face {
  unsigned int id;
  Kernel::Triangle_3 triangle;
  int material;
  std::set<unsigned int> neighbours;
  int block_id;
};

struct Block {
  std::vector<Kernel::Point_2> points;
  std::vector<Kernel::Point_2> convex_hull;
  Kernel::RT min, max;
};

typedef std::vector<Face>::iterator Faces_iterator;
typedef std::pair<Faces_iterator, Faces_iterator> Faces_intersection;
typedef CGAL::Box_intersection_d::Box_with_handle_d<Kernel::RT, 2, Faces_iterator> Box;

const std::string input_file = "/Users/jacksonye/GEO1004A1/10-282-562-obj/10-282-562-LoD13-3D.obj";
const std::string output_file = "/Users/jacksonye/GEO1004A1/10-282-562-obj/delft-LoD13-3D-2.obj";

struct Box_intersector {
  std::back_insert_iterator<std::vector<Faces_intersection>> back_inserter;

  Box_intersector(const std::back_insert_iterator<std::vector<Faces_intersection>> &bi) : back_inserter(bi) { }

  void operator()(const Box &a, const Box &b) {
    *back_inserter++ = Faces_intersection(a.handle(), b.handle());
  }
};

int main(int argc, const char *argv[]) {

  std::vector<Kernel::Point_3> input_vertices;
  std::vector<Face> input_faces;
  std::vector<Box> boxes;
  std::vector<Block> blocks;

  double expansion_threshold = 2;

  std::ifstream input_stream;
  input_stream.open(input_file);
  if (input_stream.is_open()) {
    std::string line;
    int current_material = -1;

    while (getline(input_stream, line)) {
      std::istringstream line_stream(line);
      std::string line_type;
      line_stream >> line_type;

      if (line_type == "v") {
        double x, y, z;
        line_stream >> x >> y >> z;
        input_vertices.emplace_back(Kernel::Point_3(x, y, z));
      }

      if (line_type == "f") {
        input_faces.emplace_back();
        input_faces.back().id = input_faces.size()-1;
        input_faces.back().material = current_material;
        input_faces.back().block_id = -1;
        std::vector<Kernel::Point_3> face_vertices;
        int v;
        while (!line_stream.eof()) {
          line_stream >> v;
          face_vertices.emplace_back(input_vertices[v-1]);
        }
        if (face_vertices.size() == 3) {
          // 1. Generating Triangle Geometry
          input_faces.back().triangle = Kernel::Triangle_3(
              face_vertices[0],
              face_vertices[1],
              face_vertices[2]);
        }
      }

      if (line_type == "usemtl") {
        line_stream >> current_material;
      }
    }
    input_stream.close();
  }

  for (Faces_iterator current_face = input_faces.begin(); current_face != input_faces.end(); ++current_face) {
    CGAL::Bbox_3 face_box = current_face->triangle.bbox();
    // 2. Computing 2D Expansion BBoxes
    Kernel::RT lo[2] = {
      face_box.xmin() - expansion_threshold,
      face_box.ymin() - expansion_threshold
    };
    Kernel::RT hi[2] = {
      face_box.xmax() + expansion_threshold,
      face_box.ymax() + expansion_threshold
    };
    boxes.emplace_back(lo, hi, current_face);   // 存入盒序列
  }

  std::vector<Faces_intersection> intersections;
  Box_intersector box_intersector(std::back_inserter(intersections));
  CGAL::box_self_intersection_d(boxes.begin(), boxes.end(), box_intersector);
  std::cout << intersections.size() << " intersections found" << std::endl;
  for (auto const &intersection: intersections) {
    // 3. Establishment of the adjacency table
    unsigned int id_a = intersection.first->id;
    unsigned int id_b = intersection.second->id;
    input_faces[id_a].neighbours.insert(id_b);
    input_faces[id_b].neighbours.insert(id_a);
  }

  int current_block = 0;
  std::list<unsigned int> faces_in_current_block;
  for (auto &face: input_faces) {
    if (face.block_id != -1) continue;
    face.block_id = current_block;
    faces_in_current_block.push_back(face.id);
    while (!faces_in_current_block.empty()) {
      for (auto const &neighbour: input_faces[faces_in_current_block.front()].neighbours) {
        if (input_faces[neighbour].block_id == -1) {
          input_faces[neighbour].block_id = current_block;
          faces_in_current_block.push_back(neighbour);
        }
      } faces_in_current_block.pop_front();
    } ++current_block;
  }
  std::cout << current_block << " blocks found" << std::endl;

  // 4. Create Block and fill with points & height
  blocks.resize(current_block);
  for (const auto &face : input_faces) {
    Block &blk = blocks[face.block_id];
    for (int i = 0; i < 3; ++i) {
      const Kernel::Point_3 &p = face.triangle[i];
      blk.points.emplace_back(p.x(), p.y());
      if (blk.points.size() == 1) {
        blk.min = blk.max = p.z();
      } else {
        blk.min = std::min(blk.min, p.z());
        blk.max = std::max(blk.max, p.z());
      }
    }
  }

  for (auto &block: blocks) {
    CGAL::convex_hull_2(block.points.begin(), block.points.end(), std::back_inserter(block.convex_hull));
  }

  // 5. Generation of output triangles
  std::vector<Kernel::Point_3> output_vertices;
  std::vector<Kernel::Triangle_3> output_faces;

  for (const auto &block : blocks) {
    const auto &h = block.convex_hull;
    if (h.size() < 3) continue; // 忽略退化块

    // roof (z = max)
    Kernel::Point_3 roof_origin(h[0].x(), h[0].y(), block.max);
    for (std::size_t i = 1; i + 1 < h.size(); ++i) {
      Kernel::Point_3 p1(h[i].x(),   h[i].y(),   block.max);
      Kernel::Point_3 p2(h[i+1].x(), h[i+1].y(), block.max);
      output_faces.emplace_back(roof_origin, p1, p2); // CCW 法向量向上
    }

    // floors (z = min)
    Kernel::Point_3 floor_origin(h[0].x(), h[0].y(), block.min);
    for (std::size_t i = 1; i + 1 < h.size(); ++i) {
      Kernel::Point_3 p1(h[i+1].x(), h[i+1].y(), block.min);
      Kernel::Point_3 p2(h[i].x(),   h[i].y(),   block.min);
      output_faces.emplace_back(floor_origin, p1, p2); // CW 法向量向下
    }

    // façade
    for (std::size_t i = 0; i < h.size(); ++i) {
      std::size_t j = (i + 1) % h.size();
      Kernel::Point_3 a(h[i].x(), h[i].y(), block.max);
      Kernel::Point_3 b(h[j].x(), h[j].y(), block.max);
      Kernel::Point_3 c(h[j].x(), h[j].y(), block.min);
      Kernel::Point_3 d(h[i].x(), h[i].y(), block.min);
      // (a,b,c) + (a,c,d)
      output_faces.emplace_back(a, b, c);
      output_faces.emplace_back(a, c, d);
    }
  }

  std::ofstream output_stream;
  output_stream.open(output_file);
  // 6. write .OBJ
  for (const auto &tri : output_faces) {
    for (int i = 0; i < 3; ++i) {
      const auto &p = tri[i];
      output_vertices.push_back(p);
      output_stream << "v "
                    << CGAL::to_double(p.x()) << " "
                    << CGAL::to_double(p.y()) << " "
                    << CGAL::to_double(p.z()) << "\n";
    }
  }

  for (std::size_t i = 0; i < output_vertices.size(); i += 3) {
    output_stream << "f " << i + 1 << " " << i + 2 << " " << i + 3 << "\n";
  }

  output_stream.close();
  std::cout << "Wrote " << output_faces.size() << " triangles to " << output_file << std::endl;
  return 0;
}
