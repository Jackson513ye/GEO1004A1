#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <set>
#include <list>    // For BFS queue in block assignment
#include <map>     // For vertex map in OBJ output
#include <algorithm> // For std::sort, std::unique, std::min, std::max
#include <limits>    // For std::numeric_limits
#include <iomanip>   // For std::fixed, std::setprecision

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Bbox_2.h> // Required for CGAL::Bbox_2
// CGAL::Bbox_3 is implicitly available via Kernel::Triangle_3::bbox()

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

struct Face {
    unsigned int id;
    Kernel::Triangle_3 triangle;
    int material;
    std::set<unsigned int> neighbours;
    int block_id;
};

struct Block {
    std::vector<Kernel::Point_2> points; // 2D points (x,y) of all vertices in the block
    std::vector<Kernel::Point_2> convex_hull; // 2D convex hull of the points
    Kernel::RT min_z = std::numeric_limits<Kernel::RT>::max();   // Min z-coordinate in the block
    Kernel::RT max_z = std::numeric_limits<Kernel::RT>::lowest(); // Max z-coordinate in the block
};

typedef std::vector<Face>::iterator Faces_iterator;
typedef std::pair<Faces_iterator, Faces_iterator> Faces_intersection;
// Box type for CGAL's box intersection
typedef CGAL::Box_intersection_d::Box_with_handle_d<Kernel::RT, 2, Faces_iterator> Box;

// Input and output file paths (can be modified or passed as arguments)
const std::string input_file = "/Users/ken/Downloads/10-282-562-obj/10-282-562-LoD12-3D.obj"; // Please update this path
const std::string output_file = "/Users/ken/Downloads/delft_simplified.obj"; // Please update this path

// Callback functor for box intersection
struct Box_intersector {
    std::back_insert_iterator<std::vector<Faces_intersection>> back_inserter;

    Box_intersector(const std::back_insert_iterator<std::vector<Faces_intersection>>& bi) : back_inserter(bi) {}

    // This operator is called for each pair of intersecting boxes
    void operator()(const Box& a, const Box& b) {
        // Store the pair of iterators to the Face objects whose boxes intersect
        *back_inserter++ = Faces_intersection(a.handle(), b.handle());
    }
};

int main(int argc, const char* argv[]) {

    std::vector<Kernel::Point_3> input_vertices_master_list; // Stores all vertices read from OBJ
    std::vector<Face> input_faces; // Stores all triangular faces
    std::vector<Box> boxes; // Stores 2D bounding boxes for each face
    std::vector<Block> blocks; // Stores block information

    double expansion_threshold = 1.0; // Threshold for expanding bounding boxes (in model units)

    // Read file
    std::ifstream input_stream;
    input_stream.open(input_file);
    if (input_stream.is_open()) {
        std::cout << "Reading OBJ file: " << input_file << std::endl;
        std::string line;
        int current_material = -1; // Tracks current material from "usemtl"

        // Parse line by line
        while (getline(input_stream, line)) {
            std::istringstream line_stream(line);
            std::string line_type;
            line_stream >> line_type;

            // Vertex
            if (line_type == "v") {
                double x, y, z;
                line_stream >> x >> y >> z;
                input_vertices_master_list.emplace_back(Kernel::Point_3(x, y, z));
            }

            // Face
            if (line_type == "f") {
                std::vector<Kernel::Point_3> face_polygon_vertices;
                std::string v_str;
                // Read all vertex specifications for this face line
                while (line_stream >> v_str) {
                    std::stringstream v_ss(v_str);
                    int v_idx;
                    v_ss >> v_idx; // Parses the first integer (vertex index) from "v", "v/vt", or "v/vt/vn"

                    // OBJ indices are 1-based, vector indices are 0-based
                    if (v_idx > 0 && (size_t)v_idx <= input_vertices_master_list.size()) {
                        face_polygon_vertices.emplace_back(input_vertices_master_list[v_idx - 1]);
                    }
                    else if (v_idx < 0 && (size_t)(-v_idx) <= input_vertices_master_list.size()) { // Handle relative indices
                        face_polygon_vertices.emplace_back(input_vertices_master_list[input_vertices_master_list.size() + v_idx]);
                    }
                    else {
                        std::cerr << "Warning: Invalid vertex index '" << v_str << "' in face definition. Skipping vertex." << std::endl;
                    }
                }

                // Triangulate if polygon has more than 3 vertices (simple fan triangulation)
                if (face_polygon_vertices.size() >= 3) {
                    for (size_t i = 0; i < face_polygon_vertices.size() - 2; ++i) {
                        input_faces.emplace_back(); // Create a new Face object for each triangle
                        Face& current_new_face = input_faces.back();
                        current_new_face.id = input_faces.size() - 1; // ID is its index in input_faces
                        current_new_face.material = current_material;
                        current_new_face.block_id = -1; // Initialize block_id to unassigned
                        // Create triangle using vertices (v0, v_i+1, v_i+2) for fan
                        current_new_face.triangle = Kernel::Triangle_3(face_polygon_vertices[0], face_polygon_vertices[i + 1], face_polygon_vertices[i + 2]);
                    }
                }
                else if (face_polygon_vertices.size() > 0) { // Not enough vertices for a triangle
                    std::cerr << "Warning: Face with " << face_polygon_vertices.size() << " vertices found after parsing. Skipping face." << std::endl;
                }
            }

            // Material
            if (line_type == "usemtl") {
                // For simplicity, we parse it as an int. If material names are strings, adapt this.
                // If it's a string, you might map it to an int or store the string.
                std::string material_name;
                line_stream >> material_name;
                // Example: try to convert to int if it's numeric, otherwise hash or map.
                try {
                    current_material = std::stoi(material_name);
                }
                catch (const std::invalid_argument& ia) {
                    // Handle non-integer material names if necessary, e.g., by assigning a default or using a map
                    // For now, let's assign a default if conversion fails
                    // std::cerr << "Warning: Non-integer material name '" << material_name << "'. Using default -1." << std::endl;
                    current_material = -1; // Or some other indicator
                }
            }
        }
        input_stream.close();
        std::cout << "Finished reading OBJ. " << input_vertices_master_list.size() << " vertices, " << input_faces.size() << " triangles." << std::endl;
    }
    else {
        std::cerr << "Error: Could not open input file " << input_file << std::endl;
        return 1;
    }

    // Generate expanded 2D bounding boxes for each triangle
    std::cout << "Generating expanded 2D bounding boxes..." << std::endl;
    for (Faces_iterator current_face_iter = input_faces.begin(); current_face_iter != input_faces.end(); ++current_face_iter) {
        // Project triangle vertices to 2D (XY plane)
        Kernel::Point_2 p0_2d(current_face_iter->triangle.vertex(0).x(), current_face_iter->triangle.vertex(0).y());
        Kernel::Point_2 p1_2d(current_face_iter->triangle.vertex(1).x(), current_face_iter->triangle.vertex(1).y());
        Kernel::Point_2 p2_2d(current_face_iter->triangle.vertex(2).x(), current_face_iter->triangle.vertex(2).y());

        // Compute 2D bounding box of the projected triangle
        CGAL::Bbox_2 face_box_2d = p0_2d.bbox() + p1_2d.bbox() + p2_2d.bbox();

        // Expand the 2D bounding box
        CGAL::Bbox_2 expanded_box_2d(face_box_2d.xmin() - expansion_threshold,
            face_box_2d.ymin() - expansion_threshold,
            face_box_2d.xmax() + expansion_threshold,
            face_box_2d.ymax() + expansion_threshold);

        boxes.emplace_back(expanded_box_2d, current_face_iter); // Store box with iterator to the Face
    }

    // Compute neighbours using box intersections
    std::cout << "Computing triangle neighbours based on box intersections..." << std::endl;
    std::vector<Faces_intersection> intersections;
    Box_intersector box_intersector_instance(std::back_inserter(intersections));
    CGAL::box_self_intersection_d(boxes.begin(), boxes.end(), box_intersector_instance);
    std::cout << intersections.size() << " potential intersections found." << std::endl;

    for (const auto& intersection_pair : intersections) {
        Faces_iterator face_iter_a = intersection_pair.first;
        Faces_iterator face_iter_b = intersection_pair.second;

        // Add each face to the other's set of neighbours
        // (IDs are used for storage in the set)
        if (face_iter_a->id != face_iter_b->id) { // Ensure they are not the same face
            face_iter_a->neighbours.insert(face_iter_b->id);
            face_iter_b->neighbours.insert(face_iter_a->id);
        }
    }

    // Assign blocks to faces (Connected Components Algorithm using BFS)
    std::cout << "Assigning faces to blocks..." << std::endl;
    int current_block_id_counter = 0;
    for (auto& face : input_faces) {
        if (face.block_id == -1) { // If face not yet assigned to a block
            std::list<unsigned int> faces_to_visit_in_block; // Queue for BFS

            face.block_id = current_block_id_counter;
            faces_to_visit_in_block.push_back(face.id);

            while (!faces_to_visit_in_block.empty()) {
                unsigned int current_face_id = faces_to_visit_in_block.front();
                faces_to_visit_in_block.pop_front();

                for (unsigned int neighbour_id : input_faces[current_face_id].neighbours) {
                    if (input_faces[neighbour_id].block_id == -1) { // If neighbour not yet assigned
                        input_faces[neighbour_id].block_id = current_block_id_counter;
                        faces_to_visit_in_block.push_back(neighbour_id);
                    }
                }
            }
            current_block_id_counter++; // Move to the next block ID
        }
    }
    std::cout << current_block_id_counter << " blocks found." << std::endl;

    // Create Block objects and populate them with points and min/max Z values
    if (current_block_id_counter > 0) {
        blocks.resize(current_block_id_counter);
    }

    for (const auto& face : input_faces) {
        if (face.block_id != -1 && (size_t)face.block_id < blocks.size()) {
            Block& current_block_ref = blocks[face.block_id];
            for (int i = 0; i < 3; ++i) {
                const Kernel::Point_3& p3d = face.triangle.vertex(i);
                current_block_ref.points.emplace_back(p3d.x(), p3d.y()); // Store 2D point

                // Update min/max Z for the block
                current_block_ref.min_z = std::min(current_block_ref.min_z, p3d.z());
                current_block_ref.max_z = std::max(current_block_ref.max_z, p3d.z());
            }
        }
    }

    // Compute 2D convex hulls for each block
    std::cout << "Computing convex hulls for blocks..." << std::endl;
    for (auto& block : blocks) {
        if (block.points.empty()) continue;

        // Remove duplicate points to ensure robust convex hull computation
        std::sort(block.points.begin(), block.points.end());
        block.points.erase(std::unique(block.points.begin(), block.points.end()), block.points.end());

        if (block.points.size() >= 3) {
            CGAL::convex_hull_2(block.points.begin(), block.points.end(), std::back_inserter(block.convex_hull));
        }
        else {
            // If less than 3 unique points, the "convex hull" is the points themselves.
            // Output generation will skip blocks with less than 3 hull points.
            block.convex_hull = block.points;
        }
    }

    // Create output triangles (floor, roof, walls) for each block
    std::cout << "Generating output geometry..." << std::endl;
    std::vector<Kernel::Triangle_3> output_triangles; // Stores all generated triangles for output

    for (const auto& block : blocks) {
        if (block.convex_hull.size() < 3) {
            // Need at least 3 points in the convex hull to form a polygon for roof/floor
            // std::cerr << "Warning: Block with < 3 points in convex hull. Skipping output for this block." << std::endl;
            continue;
        }

        std::vector<Kernel::Point_3> floor_3d_points;
        for (const auto& p2d : block.convex_hull) {
            floor_3d_points.emplace_back(p2d.x(), p2d.y(), block.min_z);
        }

        std::vector<Kernel::Point_3> roof_3d_points;
        for (const auto& p2d : block.convex_hull) {
            roof_3d_points.emplace_back(p2d.x(), p2d.y(), block.max_z);
        }

        // Create floor faces (triangulate the polygon - fan triangulation)
        // CGAL convex_hull_2 outputs points in counter-clockwise (CCW) order.
        // For floor (viewed from above), we want clockwise (CW) to have normal pointing downwards.
        for (size_t i = 1; i < floor_3d_points.size() - 1; ++i) {
            output_triangles.emplace_back(floor_3d_points[0], floor_3d_points[i + 1], floor_3d_points[i]);
        }

        // Create roof faces (triangulate the polygon - fan triangulation)
        // For roof (viewed from above), CCW order gives normal pointing upwards.
        for (size_t i = 1; i < roof_3d_points.size() - 1; ++i) {
            output_triangles.emplace_back(roof_3d_points[0], roof_3d_points[i], roof_3d_points[i + 1]);
        }

        // Create wall faces
        // Each edge of the convex hull polygon forms a quad wall panel, split into two triangles.
        size_t num_hull_points = block.convex_hull.size();
        for (size_t i = 0; i < num_hull_points; ++i) {
            size_t next_i = (i + 1) % num_hull_points; // Next point, wraps around

            const auto& p0_floor = floor_3d_points[i];
            const auto& p1_floor = floor_3d_points[next_i];
            const auto& p0_roof = roof_3d_points[i];
            const auto& p1_roof = roof_3d_points[next_i];

            // Wall quad: (p0_floor, p1_floor, p1_roof, p0_roof)
            // Triangle 1 (CCW for outward normal): (p0_floor, p1_floor, p1_roof)
            output_triangles.emplace_back(p0_floor, p1_floor, p1_roof);
            // Triangle 2 (CCW for outward normal): (p0_floor, p1_roof, p0_roof)
            output_triangles.emplace_back(p0_floor, p1_roof, p0_roof);
        }
    }

    // Write output OBJ file
    std::cout << "Writing simplified model to OBJ file: " << output_file << std::endl;
    std::ofstream output_stream_file; // Renamed to avoid conflict with std::ofstream
    output_stream_file.open(output_file);
    if (output_stream_file.is_open()) {
        std::vector<Kernel::Point_3> final_output_vertex_list;
        std::map<Kernel::Point_3, int> vertex_to_output_idx_map; // For unique vertices and their 1-based OBJ index

        // Helper lambda to add a vertex to the list if unique, and return its 1-based index
        auto get_vertex_output_idx =
            [&](const Kernel::Point_3& p) -> int {
            auto map_iterator = vertex_to_output_idx_map.find(p);
            if (map_iterator != vertex_to_output_idx_map.end()) {
                return map_iterator->second; // Vertex already exists, return its index
            }
            // New vertex, add to list and map
            final_output_vertex_list.push_back(p);
            int new_idx = final_output_vertex_list.size(); // New 1-based index
            vertex_to_output_idx_map[p] = new_idx;
            return new_idx;
            };

        // First pass: populate unique vertex list and map from all output triangles
        for (const auto& tri : output_triangles) {
            get_vertex_output_idx(tri.vertex(0));
            get_vertex_output_idx(tri.vertex(1));
            get_vertex_output_idx(tri.vertex(2));
        }

        // Set precision for floating point output
        output_stream_file << std::fixed << std::setprecision(6);

        // Write vertices to OBJ file
        for (const auto& v : final_output_vertex_list) {
            output_stream_file << "v " << v.x() << " " << v.y() << " " << v.z() << std::endl;
        }

        // Write faces to OBJ file
        for (const auto& tri : output_triangles) {
            // Retrieve 1-based indices for each vertex of the triangle
            int idx0 = get_vertex_output_idx(tri.vertex(0));
            int idx1 = get_vertex_output_idx(tri.vertex(1));
            int idx2 = get_vertex_output_idx(tri.vertex(2));
            output_stream_file << "f " << idx0 << " " << idx1 << " " << idx2 << std::endl;
        }

        output_stream_file.close();
        std::cout << "Successfully written simplified model." << std::endl;
        std::cout << "Output model stats: " << final_output_vertex_list.size() << " vertices, " << output_triangles.size() << " faces." << std::endl;

    }
    else {
        std::cerr << "Error: Could not open output file " << output_file << std::endl;
        return 1;
    }

    return 0;
}