/*
 *
 * This file contains the complete source code of the graph analysis library. While having one
 *   single file is not the most clear and organized choice for source codes, this format has
 *   been chosen to faciliate portability of the code.
 *
 * Compilation example: g++ -O3 my_code.cpp
 * 
 * Author:  Antoine Allard
 * WWW:     antoineallard.info
 * Date:    April 2018
 * 
 * Version: 1.0.1
 * 
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * 
 */

// Standard Template Library
#include <algorithm>         // std::min, std::fill
#include <cstdlib>           // std::terminate
#include <fstream>           // std::fstream
#include <iomanip>           // std::setw
#include <iostream>          // std::cerr
#include <map>               // std::map
#include <set>               // std::set, std::multiset
#include <sstream>           // std::stringstream
#include <string>            // std::string, std::getline
#include <utility>           // std::pair, std::make_pair, std::swap
#include <vector>            // std::vector


class undirected_graph_t
{
  // =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
  // Internal storage of the graph object.
  public:
    // Name to ID conversion. Every vertex is assigned a numerical ID in [0, |V|).
    std::map<std::string, int> Name2ID;
    // Edgelist.
    std::set< std::pair<int, int> > edgelist;
    // Adjacency list.
    std::vector<std::vector<int> > adjacency_list;
  // =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
  // Typedefs and accessors.
  public:
    typedef std::set< std::pair<int, int> >::iterator edgelist_iterator;
    edgelist_iterator edgelist_begin() { return edgelist.begin(); }
    edgelist_iterator edgelist_end()   { return edgelist.end();   }
  // =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
  // Objects/functions related to outputs.
  public:
    // Custom types to indicate how to identify vertices in outputs.
    enum vID_t { vID_none, vID_num, vID_name };
    static const bool header_true = true;
    static const bool header_false = false;
    // ID to name conversion.
    std::vector< std::string > ID2Name;
    // Build the ID2Name vector.
    void build_ID2Name();
  private:
    // Default width of columns.
    static const int default_column_width = 15;
    // Available vertex properties.
    std::set< std::string > available_vertex_prop;
    // Available vertex integer properties.
    std::set< std::string > available_vertex_integer_prop;
    // Default headers associated with vertex properties.
    std::map< std::string, std::string > v_prop_header;
  // =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
  // Properties.
  public:
    // Properties of the graph.
    std::map< std::string, double> g_prop;
    // Properties of vertices.
    std::map< std::string, std::vector<double> > v_prop;
    // List of all triangles.
    std::vector< std::vector<int> > triangles;
    // Shortest path lengths.
    int shortest_path_length(int v1, int v2);
    std::vector<double> average_shortest_path_lengths;
    std::vector< std::vector<int> > shortest_path_lengths;
    // Diameters of the components.
    std::vector<double> diameters;
    // Unique integer properties of vertices.
    std::map< std::string, std::set<int> > v_class_1p;
    std::map< std::pair<std::string, std::string>, std::set< std::pair<int, int> > > v_class_2p;
    // Number of vertices in each unique integer properties of vertices.
    std::map< std::string, std::map<int, int> > v_class_count_1p;
    std::map< std::pair<std::string, std::string>, std::map<std::pair<int, int>, int> > v_class_count_2p;
    // Number of edges in each unique integer properties of edges.
    std::map< std::string, std::map< std::pair<int, int>, int> > e_class_count_1p;
    std::map< std::pair< std::string, std::string >, std::map< std::multiset< std::pair<int, int> >, int> > e_class_count_2p;
  private:
    // Function verifying wether a vertex property exists.
    void is_vertex_property(std::string prop);
    // Function verifying wether a vertex property qualifies as a "integer" property.
    void is_vertex_integer_property(std::string prop);
  // =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
  // Functions loading/saving edgelists and vertices/edges properties.
  public:
    // Loads the graph structure from an edgelist in a file.
    void load_edgelist(std::string edgelist_filename);
    // Outputing the vertices properties (the last 3 inputs can be omitted and/or put in any order).
    void save_vertices_properties(std::string filename, std::vector<std::string> props_id, vID_t vID = vID_name, int width = default_column_width, bool header = header_true);
    void save_vertices_properties(std::string filename, std::vector<std::string> props_id, vID_t vID,            bool header,                      int width = default_column_width)                { save_vertices_properties(filename, props_id, vID, width, header); };
    void save_vertices_properties(std::string filename, std::vector<std::string> props_id, bool header,          vID_t vID = vID_name,             int width = default_column_width)                { save_vertices_properties(filename, props_id, vID, width, header); };
    void save_vertices_properties(std::string filename, std::vector<std::string> props_id, bool header,          int width,                        vID_t vID = vID_name)                            { save_vertices_properties(filename, props_id, vID, width, header); };
    void save_vertices_properties(std::string filename, std::vector<std::string> props_id, int width,            vID_t vID = vID_name,             bool header = header_true)                       { save_vertices_properties(filename, props_id, vID, width, header); };
    void save_vertices_properties(std::string filename, std::vector<std::string> props_id, int width,            bool header,                      vID_t vID = vID_name)                            { save_vertices_properties(filename, props_id, vID, width, header); };
    // Outputing the number of vertices in each vertex class.
    void save_number_vertices_in_classes(std::string filename, std::string prop, bool header = header_true, int width = default_column_width);
    void save_number_vertices_in_classes(std::string filename, std::string prop, int width,                 bool header = header_true)         { save_number_vertices_in_classes(filename, prop, header, width); };
    void save_number_vertices_in_classes(std::string filename, std::pair<std::string, std::string> props, bool header = header_true, int width = default_column_width);
    void save_number_vertices_in_classes(std::string filename, std::pair<std::string, std::string> props, int width,                 bool header = header_true)         { save_number_vertices_in_classes(filename, props, header, width); };
    // Outputing the number of edges in each vertex class.
    void save_number_edges_in_classes(std::string filename, std::string prop, bool header = header_true, int width = default_column_width);
    void save_number_edges_in_classes(std::string filename, std::string prop, int width, bool header = header_true)                         { save_number_edges_in_classes(filename, prop, header, width); };
    void save_number_edges_in_classes(std::string filename, std::pair<std::string, std::string> props, bool header = header_true, int width = default_column_width);
    void save_number_edges_in_classes(std::string filename, std::pair<std::string, std::string> props, int width, bool header = header_true)                         { save_number_edges_in_classes(filename, props, header, width); };
  // =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
  // Functions building secondary objects related to the graph.
  private:
    // Builds the adjacency list.
    void build_adjacency_list();
  // =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
  // Functions extracting properties of graph/vertices/edges.
  public:
    // Computes the average degree of neighbors. It can be normalized or not.
    void average_degree_of_neighbors(bool normalized = false);
    // Computes the average value of a vertex property.
    void average_vertex_prop(std::string prop);
    // Computes the average value of a vertex property of neighbors.
    void average_vertex_prop_of_neighbors(std::string prop);
    // Computes the closeness centrality.
    void closeness_centrality(bool correct_for_component_size = true);
    // Computes the degree correlation coefficient.
    void degree_correlation_coefficient();
    // Extracts the degree of vertices.
    void degrees();
    // Computes the harmonic centrality.
    void harmonic_centrality();
    // Extracts the k-core decomposition.
    void kcore_decomposition();
    // Extracts the local clustering coefficients.
    void local_clustering_coefficients();
    // Counts the number of triangles around every vertex.
    void number_of_triangles();
    // Extracts the onion decomposition (includes the k-core decomposition as a by-product).
    void onion_decomposition();
    // Identifies the connected components to which vertices belong.
    void survey_connected_components();
    // Finds the length of all shortest paths (equal to the number of vertices if no path exists).
    void survey_shortest_path_lengths();
    // Compiles a list of all triangles in the graph.
    void survey_triangles();
    // Computes the average shortest path lengtha and the diameters.
    void topological_dimensions();
    // Identifies unique integer properties of vertices (ex. degree or degree-layer).
    void identify_classes_of_vertices(std::string prop);
    void identify_classes_of_vertices(std::pair<std::string, std::string> props);
    // Counts the number of vertices in each class of vertices.
    void count_vertices_in_classes(std::string prop);
    void count_vertices_in_classes(std::pair<std::string, std::string> props);
    // Counts the number of edges in each class of edges.
    void count_edges_in_classes(std::string prop);
    void count_edges_in_classes(std::pair<std::string, std::string> props);
  private:
    // Counts the number of triangles around a given vertex.
    int count_triangles_around_vertex(int v1);
    // Functions associated to the extraction of the components.
    int get_root(int i, std::vector<int> &clust_id);
    void merge_clusters(std::vector<int> &size, std::vector<int> &clust_id, int nb_vertices);
  // =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
  // Constructors (and related functions).
  private:
    // Function setting default values (to avoid requiring the C++11 standard).
    void initialization();
  public:
    // Empty constructor.
    undirected_graph_t() { initialization(); };
    // Constructor with edgelist.
    undirected_graph_t(std::string edgelist_filename) { initialization(); load_edgelist(edgelist_filename); };
};










// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::initialization()
{
  // Available vertex properties.
  available_vertex_prop.insert("degrees");
  available_vertex_prop.insert("local_clust");
  available_vertex_prop.insert("kcore");
  available_vertex_prop.insert("od_layer");
  available_vertex_prop.insert("nb_triangles");
  available_vertex_prop.insert("conn_comp");
  available_vertex_prop.insert("close_cent");
  available_vertex_prop.insert("corr_close_cent");
  available_vertex_prop.insert("harmo_cent");
  available_vertex_prop.insert("avg_neigh_degrees");
  available_vertex_prop.insert("avg_neigh_degrees_norm");
  available_vertex_prop.insert("avg_neigh_local_clust");
  available_vertex_prop.insert("avg_neigh_kcore");
  available_vertex_prop.insert("avg_neigh_od_layer");
  available_vertex_prop.insert("avg_neigh_nb_triangles");
  // Available vertex properties.
  available_vertex_integer_prop.insert("degrees");
  available_vertex_integer_prop.insert("kcore");
  available_vertex_integer_prop.insert("od_layer");
  available_vertex_integer_prop.insert("nb_triangles");
  available_vertex_integer_prop.insert("conn_comp");
  // Headers for properties.
  v_prop_header["degrees"]                = "Degree";
  v_prop_header["local_clust"]            = "Clust.";
  v_prop_header["kcore"]                  = "Coreness";
  v_prop_header["od_layer"]               = "ODlayer";
  v_prop_header["nb_triangles"]           = "NbTriang.";
  v_prop_header["conn_comp"]              = "Conn.Comp.";
  v_prop_header["close_cent"]             = "Close.Cent.";
  v_prop_header["corr_close_cent"]        = "CClose.Cent.";
  v_prop_header["harmo_cent"]             = "Harmo.Cent.";
  v_prop_header["avg_neigh_degrees"]      = "ANDegree";
  v_prop_header["avg_neigh_degrees_norm"] = "ANDegreeNorm";
  v_prop_header["avg_neigh_local_clust"]  = "ANClust.";
  v_prop_header["avg_neigh_kcore"]        = "ANCoreness";
  v_prop_header["avg_neigh_od_layer"]     = "ANODlayer";
  v_prop_header["avg_neigh_nb_triangles"] = "ANNbTriang.";
  // // Sets some default values.
  // header_true = true;
  // header_false = false;
  // default_column_width = 15;
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::average_degree_of_neighbors(bool normalized)
{
  // Computes the average degree of neighbors.
  average_vertex_prop_of_neighbors("degrees");

  if(normalized)
  {
    // =============================================================================================
    // Initializes relevant objects of the class.
    int nb_vertices = g_prop["nb_vertices"];

    std::vector<double>& Vertex2Degree = v_prop["degrees"];

    v_prop["avg_neigh_degrees_norm"].clear();
    std::vector<double>& Vertex2ANDegree = v_prop["avg_neigh_degrees_norm"];
    Vertex2ANDegree.resize(nb_vertices, 0);

    std::vector<double>& Vertex2ADegree = v_prop["avg_neigh_degrees"];
    // =============================================================================================

    // Extracts the degrees of the vertices.
    if(Vertex2Degree.size() != nb_vertices)
    {
      degrees();
    }

    // Computes the two first "moments" of the degree distribution.
    double m1 = 0;
    double m2 = 0;
    double d;
    for(int v(0); v<nb_vertices; ++v)
    {
      d = Vertex2Degree[v];
      m1 += d;
      m2 += d * d;
    }

    // Compute the normalized values of the average degree of neighbors.
    for(int v(0); v<nb_vertices; ++v)
    {
      Vertex2ANDegree[v] = Vertex2ADegree[v] * m1 / m2;
    }
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::average_vertex_prop(std::string prop)
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  std::vector<double>& Vertex2Prop = v_prop[prop];
  // ===============================================================================================

  // Checks if the property has been extracted/computed already.
  if(v_prop[prop].size() != nb_vertices)
  {
    std::cerr << "ERROR: The property " << prop << " has not been extracted/computed." << std::endl;
    std::terminate();
  }

  // Computes the average degree of vertices.
  double sum_of_values = 0;
  for(int v(0); v<nb_vertices; ++v)
  {
    sum_of_values += Vertex2Prop[v];
  }

  // ===============================================================================================
  // Updates the properties of the graph.
  std::string avg_v_prop = "avg_" + prop;
  g_prop[avg_v_prop] = sum_of_values / nb_vertices;
  // ===============================================================================================
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::average_vertex_prop_of_neighbors(std::string prop)
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];

  std::vector<double>& Vertex2Prop = v_prop[prop];

  std::string avg_neigh_prop = "avg_neigh_" + prop;
  v_prop[avg_neigh_prop].clear();
  std::vector<double>& Vertex2ANProp = v_prop[avg_neigh_prop];
  Vertex2ANProp.resize(nb_vertices, 0);
  // ===============================================================================================

  // Ensures that the adjacency list has been built.
  if(adjacency_list.size() != nb_vertices)
  {
    build_adjacency_list();
  }

  // Computes the average vertex property of neighbors.
  int d1, v2;
  double value;
  for(int v1(0); v1<nb_vertices; ++v1)
  {
    value = 0;
    d1 = adjacency_list[v1].size();
    for(int s(0); s<d1; ++s)
    {
      v2 = adjacency_list[v1][s];
      value += Vertex2Prop[v2];
    }
    Vertex2ANProp[v1] = value / d1;
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::build_adjacency_list()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  adjacency_list.clear();
  adjacency_list.resize(nb_vertices);
  // ===============================================================================================


  // Loops over all edges.
  int v1, v2;
  std::set< std::pair<int, int> >::iterator it = edgelist.begin();
  std::set< std::pair<int, int> >::iterator end = edgelist.end();
  for(; it!=end; ++it)
  {
    // Identifies the vertices.
    v1 = it->first;
    v2 = it->second;
    // Adds the percolated edge.
    adjacency_list[v1].push_back(v2);
    adjacency_list[v2].push_back(v1);
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::build_ID2Name()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  ID2Name.clear();
  ID2Name.resize(nb_vertices);
  // ===============================================================================================


  // Loops over all names.
  std::map<std::string, int>::iterator it = Name2ID.begin();
  std::map<std::string, int>::iterator end = Name2ID.end();
  for(; it!=end; ++it)
  {
    ID2Name[it->second] = it->first;
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::closeness_centrality(bool correct_for_component_size)
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  std::string v_prop_name;
  if(correct_for_component_size)
  {
    v_prop_name = "corr_close_cent";
  }
  else
  {
    v_prop_name = "close_cent";
  }
  std::vector<double>& Vertex2Cent = v_prop[v_prop_name];
  Vertex2Cent.clear();
  Vertex2Cent.resize(nb_vertices, 0);
  // ===============================================================================================


  // Ensures that the shortest path lengths have been computed.
  if(shortest_path_lengths.size() != nb_vertices)
  {
    survey_shortest_path_lengths();
  }

  // Compiles the closeness centrality of every vertex (ignores paths that does not exist).
  int length, cnt;
  double value;
  for(int v1(0); v1<nb_vertices; ++v1)
  {
    cnt = 0;
    value = 0;
    for(int v2(0); v2<nb_vertices; ++v2)
    {
      length = shortest_path_length(v1, v2);
      if(length != nb_vertices)
      {
        value += length;
        ++cnt;
      }
    }
    if(correct_for_component_size)
    {
      // There is also the Wasserman and Faust correction.
      Vertex2Cent[v1] = ((cnt - 1) / ((double) (nb_vertices - 1))) * ((cnt - 1) / value);
    }
    else
    {
      Vertex2Cent[v1] = (cnt - 1) / value;
    }
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::count_edges_in_classes(std::string prop)
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  std::vector<double>& Vertex2Prop = v_prop[prop];

  e_class_count_1p[prop].clear();
  std::map< std::pair<int, int>, int>& PropCount = e_class_count_1p[prop];
  // ===============================================================================================

  // Loops over every edges and counts the number of edges in every unique class of edges.
  int p1, p2;
  std::pair<int, int> p;
  std::set< std::pair<int, int> >::iterator it = edgelist.begin();
  std::set< std::pair<int, int> >::iterator end = edgelist.end();
  for(; it!=end; ++it)
  {
    p1 = Vertex2Prop[it->first];
    p2 = Vertex2Prop[it->second];
    if(p1 < p2)
    {
      p = std::make_pair(p1, p2);
    }
    else
    {
      p = std::make_pair(p2, p1);
    }
    if(PropCount.find(p) == PropCount.end())
    {
      PropCount[p] = 1;
    }
    else
    {
      PropCount[p] += 1;
    }
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::count_edges_in_classes(std::pair<std::string, std::string> props)
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  std::vector<double>& Vertex2Prop1 = v_prop[props.first];
  std::vector<double>& Vertex2Prop2 = v_prop[props.second];

  e_class_count_2p[props].clear();
  std::map< std::multiset< std::pair<int, int> >, int>& PropCount = e_class_count_2p[props];
  // ===============================================================================================

  // Loops over every edges and counts the number of edges in every unique class of edges.
  int p1, p2;
  std::multiset< std::pair<int, int> > p;
  std::multiset< std::pair<int, int> >::iterator it = edgelist.begin();
  std::multiset< std::pair<int, int> >::iterator end = edgelist.end();
  for(; it!=end; ++it)
  {
    p.clear();
    p1 = Vertex2Prop1[it->first];
    p2 = Vertex2Prop2[it->first];
    p.insert(std::make_pair(p1, p2));
    p1 = Vertex2Prop1[it->second];
    p2 = Vertex2Prop2[it->second];
    p.insert(std::make_pair(p1, p2));
    if(PropCount.find(p) == PropCount.end())
    {
      PropCount[p] = 1;
    }
    else
    {
      PropCount[p] += 1;
    }
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
int undirected_graph_t::count_triangles_around_vertex(int v1)
{
  // Variables.
  int v2, v3, d1, d2;
  int nb_triangles = 0;
  // Vector objects.
  std::vector<int> intersection;
  // Set objects.
  std::set<int> neighbours_v1, neighbours_v2;
  // Iterator objects.
  std::vector<int>::iterator it;
  // Degree of vertex v1.
  d1 = adjacency_list[v1].size();
  // Performs the calculation only if d1>1.
  if( d1 > 1 )
  {
    // Builds an ordered list of the neighbourhood of v1
    neighbours_v1.clear();
    neighbours_v1.insert(adjacency_list[v1].begin(), adjacency_list[v1].end());
    // Loops over the neighbours of vertex v1.
    for(int n1(0); n1<d1; ++n1)
    {
      // Identity and degree of vertex 2.
      v2 = adjacency_list[v1][n1];
      d2 = adjacency_list[v2].size();
      // Performs the calculation only if d2>1.
      if( d2 > 1 )
      {        
        // Builds an ordered list of the neighbourhood of v2
        neighbours_v2.clear();
        for(int n2(0); n2<d2; ++n2)
        {
          // Identifies the neighbors.
          v3 = adjacency_list[v2][n2];
          if(v2 < v3) // Ensures that triangles will be counted only once.
          {
            neighbours_v2.insert(v3);
          }
        }
        // Identifies the triangles.
        d2 = neighbours_v2.size();
        intersection.clear();
        intersection.resize(std::min(d1, d2));
        it = std::set_intersection(neighbours_v1.begin(), neighbours_v1.end(), neighbours_v2.begin(), neighbours_v2.end(), intersection.begin());
        // intersection.resize(it-intersection.begin());
        // Computes the local clustering coeficient.
        // nb_triangles += intersection.size();
        nb_triangles += it-intersection.begin();
      }
    }
    // Returns the number of triangles.
    return nb_triangles;
  }
  else
  {
    return 0;
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::count_vertices_in_classes(std::string prop)
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];

  std::vector<double>& Vertex2Prop = v_prop[prop];

  std::set<int>& PropClass = v_class_1p[prop];

  v_class_count_1p[prop].clear();
  std::map<int, int>& PropCount = v_class_count_1p[prop];
  // ===============================================================================================

  // Ensures that the relevant classes of vertices have been identified.
  if(PropClass.size() == 0)
  {
    identify_classes_of_vertices(prop);
  }

  // Initializes the count structure.
  std::set<int>::iterator it = PropClass.begin();
  std::set<int>::iterator end = PropClass.end();
  for(; it!=end; ++it)
  {
    PropCount[*it] = 0;
  }

  // Loops over every vertices and counts the number of vertices in every unique class of vertex.
  int p;
  for(int v(0); v<nb_vertices; ++v)
  {
    p = Vertex2Prop[v];
    PropCount[p] += 1;
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::count_vertices_in_classes(std::pair<std::string, std::string> props)
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];

  std::vector<double>& Vertex2Prop1 = v_prop[props.first];
  std::vector<double>& Vertex2Prop2 = v_prop[props.second];

  v_class_count_2p[props].clear();
  std::map<std::pair<int, int>, int>& PropCount = v_class_count_2p[props];

  std::set< std::pair<int, int> >& PropClass = v_class_2p[props];
  // ===============================================================================================

  // Ensures that the relevant classes of vertices have been identified.
  if(PropClass.size() == 0)
  {
    identify_classes_of_vertices(props);
  }

  // Initializes the count structure.
  std::set< std::pair<int, int> >::iterator it = PropClass.begin();
  std::set< std::pair<int, int> >::iterator end = PropClass.end();
  for(; it!=end; ++it)
  {
    PropCount[*it] = 0;
  }

  // Loops over every vertices and counts the number of vertices in every unique class of vertex.
  int p1, p2;
  for(int v(0); v<nb_vertices; ++v)
  {
    p1 = Vertex2Prop1[v];
    p2 = Vertex2Prop2[v];
    PropCount[std::make_pair(p1, p2)] += 1;
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::degree_correlation_coefficient()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  int nb_edges = g_prop["nb_edges"];
  std::set<int>& DegreeClasses = v_class_1p["degrees"];
  // ===============================================================================================

  // Extracts the degrees of the vertices.
  if(v_prop["degrees"].size() != nb_vertices)
  {
    degrees();
  }

  // Computes the average degree of vertices.
  average_vertex_prop("degrees");

  // Compiles the degree distribution.
  identify_classes_of_vertices("degrees");
  count_vertices_in_classes("degrees");

  // Computes the excess degree distribution and its variance.
  int k;
  double avg_degree = g_prop["avg_degrees"];
  double var_q, value;
  double m1 = 0;
  double m2 = 0;
  std::map<int, double> q;
  // std::map<int, double> e;
  std::map<int,int>::iterator it1 = v_class_count_1p["degrees"].begin();
  std::map<int,int>::iterator end1 = v_class_count_1p["degrees"].end();
  for(; it1!=end1; ++it1)
  {
    k = it1->first;
    value = k * it1->second / avg_degree / nb_vertices;
    q[k - 1] = value;
    // e[k - 1] = 0;
    m2 += (k - 1) * (k - 1) * value;
    m1 += (k - 1) * value;
  }
  var_q = m2 - (m1 * m1);

  // Extracts the degree-degree correlation matrix.
  count_edges_in_classes("degrees");

  // Computes the degree correlation coefficient.
  double k1, k2;
  double corr_coef = 0;
  std::map< std::pair<int, int>, int>::iterator it2 = e_class_count_1p["degrees"].begin();
  std::map< std::pair<int, int>, int>::iterator end2 = e_class_count_1p["degrees"].end();
  for(; it2!=end2; ++it2)
  {
    k1 = it2->first.first;
    k2 = it2->first.second;
    corr_coef += (k1 - 1) * (k2 - 1) * it2->second / nb_edges;
  }
  std::map<int, double>::iterator it3 = q.begin();
  std::map<int, double>::iterator end3 = q.end();
  std::map<int, double>::iterator it4;
  std::map<int, double>::iterator end4 = q.end();
  for(; it3!=end3; ++it3)
  {
    it4 = it3;
    corr_coef -= it3->first * it4->first * it3->second * it4->second;
    for(++it4; it4!=end4; ++it4)
    {
      corr_coef -= 2 * it3->first * it4->first * it3->second * it4->second;
    }
  }

  // ===============================================================================================
  // Updates the properties of the graph.
  g_prop["degree_correlation_coefficient"] = corr_coef / var_q;
  // ===============================================================================================
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::degrees()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  v_prop["degrees"].clear();
  std::vector<double>& Vertex2Degree = v_prop["degrees"];
  Vertex2Degree.resize(nb_vertices, 0);
  // ===============================================================================================


  // Extracts the degrees from the adjacency list (if already generated).
  if(adjacency_list.size() == nb_vertices)
  {
    // Loops over all vertices.
    for(int v(0); v<nb_vertices; ++v)
    {
      Vertex2Degree[v] = adjacency_list[v].size();
    }
  }
  // Otherwise extracts the degrees from the edgelist.
  else
  {
    // Loops over all edges.
    int v1, v2;
    std::set< std::pair<int, int> >::iterator it = edgelist.begin();
    std::set< std::pair<int, int> >::iterator end = edgelist.end();
    for(; it!=end; ++it)
    {
      // Identifies the vertices.
      v1 = it->first;
      v2 = it->second;
      // Adds the contribution to the degree of vertices.
      Vertex2Degree[v1] += 1;
      Vertex2Degree[v2] += 1;
    }
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
int undirected_graph_t::get_root(int i, std::vector<int> &clust_id)
{
  while(i != clust_id[i])
  {
    clust_id[i] = clust_id[clust_id[i]];
    i = clust_id[i];
  }
  return i;
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::harmonic_centrality()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  std::vector<double>& Vertex2Cent = v_prop["harmo_cent"];
  Vertex2Cent.clear();
  Vertex2Cent.resize(nb_vertices, 0);
  // ===============================================================================================


  // Ensures that the shortest path lengths have been computed.
  if(shortest_path_lengths.size() != nb_vertices)
  {
    survey_shortest_path_lengths();
  }

  // Compiles the harmonic centrality of every vertex.
  int length;
  double value;
  for(int v1(0); v1<nb_vertices; ++v1)
  {
    value = 0;
    for(int v2(0); v2<nb_vertices; ++v2)
    {
      length = shortest_path_length(v1, v2);
      if(length != nb_vertices && length != 0)
      {
        value += 1 / (double) length;
      }
    }
    Vertex2Cent[v1] = value / (nb_vertices - 1);
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::identify_classes_of_vertices(std::string prop)
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];

  std::vector<double>& Vertex2Prop = v_prop[prop];

  v_class_1p[prop].clear();
  std::set<int>& PropClass = v_class_1p[prop];
  // ===============================================================================================

  is_vertex_property(prop);
  is_vertex_integer_property(prop);
  if(v_prop[prop].size() != nb_vertices)
  {
    std::cerr << "ERROR: The property " << prop << " has not been extracted/computed." << std::endl;
    std::terminate();
  }

  // Loops over every vertices and identifies all unique class of vertex.
  int p;
  for(int v(0); v<nb_vertices; ++v)
  {
    p = Vertex2Prop[v];
    PropClass.insert(p);
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::identify_classes_of_vertices(std::pair<std::string, std::string> props)
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];

  std::vector<double>& Vertex2Prop1 = v_prop[props.first];
  std::vector<double>& Vertex2Prop2 = v_prop[props.second];

  v_class_2p[props].clear();
  std::set< std::pair<int, int> >& PropClass = v_class_2p[props];
  // ===============================================================================================

  is_vertex_property(props.first);
  is_vertex_integer_property(props.first);
  if(v_prop[props.first].size() != nb_vertices)
  {
    std::cerr << "ERROR: The property " << props.first << " has not been extracted/computed." << std::endl;
    std::terminate();
  }

  is_vertex_property(props.second);
  is_vertex_integer_property(props.second);
  if(v_prop[props.second].size() != nb_vertices)
  {
    std::cerr << "ERROR: The property " << props.second << " has not been extracted/computed." << std::endl;
    std::terminate();
  }

  // Loops over every vertices and identifies all unique class of vertex.
  int p1, p2;
  for(int v(0); v<nb_vertices; ++v)
  {
    p1 = Vertex2Prop1[v];
    p2 = Vertex2Prop2[v];
    PropClass.insert(std::make_pair(p1, p2));
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::is_vertex_integer_property(std::string prop)
{
  if(available_vertex_integer_prop.find(prop) == available_vertex_integer_prop.end())
  {
    std::cerr << "ERROR: " << prop << " is not a valid vertex property." << std::endl;
    std::terminate();
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::is_vertex_property(std::string prop)
{
  if(available_vertex_prop.find(prop) == available_vertex_prop.end())
  {
    std::cerr << "ERROR: " << prop << " is not a valid vertex property." << std::endl;
    std::terminate();
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::kcore_decomposition()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];

  v_prop["kcore"].clear();
  std::vector<double>& Vertex2Coreness = v_prop["kcore"];
  Vertex2Coreness.resize(nb_vertices, 0);
  
  std::vector<double>& Vertex2Degree = v_prop["degrees"];
  // ===============================================================================================

  // Ensures that the adjacency list has been built.
  if(adjacency_list.size() != nb_vertices)
  {
    build_adjacency_list();
  }

  // Extracts the degree (if not already generated).
  if(Vertex2Degree.size() != nb_vertices)
  {
    degrees();
  }

  // Builds two lists (std::vector, std::set) of the degree of the vertices.
  std::vector<int> DegreeVec(nb_vertices);
  std::set<std::pair<int, int> > DegreeSet;
  for(int v(0); v<nb_vertices; ++v)
  {
    DegreeSet.insert(std::make_pair(Vertex2Degree[v], v));
    DegreeVec[v] = Vertex2Degree[v];
  }

  // Determines the coreness and the layer based on the algorithm of Batagelj and Zaversnik.
  int v1, v2, d1, d2;
  std::set< std::pair<int, int> >::iterator m_it;
  while(DegreeSet.size() > 0)
  {
    // Sets the coreness of the first vertex in the list and then removes it.
    m_it = DegreeSet.begin();
    d1 = m_it->first;
    v1 = m_it->second;
    Vertex2Coreness[v1] = d1;
    DegreeSet.erase(m_it);
    // Reduces the "effective" degree of its neighbours.
    for(int n(0), nn(adjacency_list[v1].size()); n<nn; ++n)
    {
      // Identifies the neighbor.
      v2 = adjacency_list[v1][n];
      d2 = DegreeVec[v2];
      // Finds the neighbor in the list "effective" degrees.
      m_it = DegreeSet.find(std::make_pair(d2, v2));
      if(m_it != DegreeSet.end())
      {
        if(d2 > d1)
        {
          DegreeVec[v2] = d2 - 1;
          DegreeSet.erase(m_it);
          DegreeSet.insert(std::make_pair(d2 - 1, v2));
        }
      }
    }
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::load_edgelist(std::string edgelist_filename)
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  Name2ID.clear();
  edgelist.clear();
  // ===============================================================================================


  // Stream objects.
  std::fstream edgelist_file;
  std::stringstream one_line;
  // Variables.
  int v1, v2;
  // Initializes the number of vertices.
  int nb_vertices = 0;
  // String objects.
  std::string full_line, name1_str, name2_str;
  // Iterators.
  std::map< std::string, int >::iterator name_it;


  // Opens the stream and terminates if the operation did not succeed.
  edgelist_file.open(edgelist_filename.c_str(), std::fstream::in);
  if( !edgelist_file.is_open() )
  {
    std::cerr << "ERROR: Could not open file: " << edgelist_filename << "." << std::endl;
    std::terminate();
  }
  else
  {
    // Reads the OD file line by line.
    while( !edgelist_file.eof() )
    {
      // Reads a line of the file.
      std::getline(edgelist_file, full_line);
      edgelist_file >> std::ws;
      one_line.str(full_line);
      one_line >> std::ws;
      one_line >> name1_str >> std::ws;
      // Skips a line of comment.
      if(name1_str == "#")
      {
        one_line.clear();
        continue;
      }
      one_line >> name2_str >> std::ws;
      one_line.clear();

      // Ignores self-loops.
      if(name1_str != name2_str)
      {

        // Is name1 new?
        name_it = Name2ID.find(name1_str);
        if(name_it == Name2ID.end())
        {
          Name2ID[name1_str] = nb_vertices;
          v1 = nb_vertices;
          ++nb_vertices;
        }
        else
        {
          v1 = name_it->second;
        }

        // Is name2 new?
        name_it = Name2ID.find(name2_str);
        if(name_it == Name2ID.end())
        {
          Name2ID[name2_str] = nb_vertices;
          v2 = nb_vertices;
          ++nb_vertices;
        }
        else
        {
          v2 = name_it->second;
        }

        // Adds the edge (multiedges are ignored).
        if(v1 > v2)
        {
          std::swap(v1, v2);
        }
        edgelist.insert(std::make_pair(v1, v2));

      }
    }
  }
  // Closes the stream.
  edgelist_file.close();

  // ===============================================================================================
  // Updates the properties of the graph.
  g_prop["nb_vertices"] = nb_vertices;
  g_prop["nb_edges"] = edgelist.size();
  // ===============================================================================================
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::local_clustering_coefficients()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];

  v_prop["local_clust"].clear();
  std::vector<double>& Vertex2Prop = v_prop["local_clust"];
  Vertex2Prop.resize(nb_vertices, 0);
  // ===============================================================================================

  // Ensures that the adjacency list has been built.
  if(adjacency_list.size() != nb_vertices)
  {
    build_adjacency_list();
  }

  // Computes the intersection for the in- and out- neighbourhoods of each node.
  int d1;
  double nb_triangles, tmp_value;
  double average_local_clustering_coefficient = 0;
  for(int v1(0); v1<nb_vertices; ++v1)
  {
    // The degree of the vertex.
    d1 = adjacency_list[v1].size();
    // Counts the number of triangle.
    nb_triangles = count_triangles_around_vertex(v1);
    tmp_value = 0;
    // Computes the coefficient of clustering.
    if(d1 > 1)
    {
      tmp_value = 2.0 * nb_triangles / d1 / (d1 - 1);
      Vertex2Prop[v1] = tmp_value;
      average_local_clustering_coefficient += tmp_value;
    }
    else
    {
      Vertex2Prop[v1] = 0;
    }
  }
  // ===============================================================================================
  // Updates the properties of the graph.
  g_prop["avg_local_clust"] = average_local_clustering_coefficient / g_prop["nb_vertices"];
  // ===============================================================================================
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::merge_clusters(std::vector<int> &size, std::vector<int> &clust_id, int nb_vertices)
{
  // Variables.
  int v1, v2, v3, v4;
  // Loops over the vertices.
  for(int i(0); i<nb_vertices; ++i)
  {
    // Loops over the neighbors.
    for(int j(0), jj(adjacency_list[i].size()); j<jj; ++j)
    {
      if(get_root(i, clust_id) != get_root(adjacency_list[i][j], clust_id))
      {
        // Adjust the root of vertices.
        v1 = i;
        v2 = adjacency_list[i][j];
        if(size[v2] > size[v1])
          std::swap(v1, v2);
        v3 = get_root(v1, clust_id);
        v4 = get_root(v2, clust_id);
        clust_id[v4] = v3;
        size[v3] += size[v4];
      }
    }
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::number_of_triangles()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  v_prop["nb_triangles"].clear();
  std::vector<double>& Vertex2Prop = v_prop["nb_triangles"];
  Vertex2Prop.resize(nb_vertices, 0);
  // ===============================================================================================


  // Ensures that the adjacency list has been built.
  if(adjacency_list.size() != nb_vertices)
  {
    build_adjacency_list();
  }

  // Counts the number of triangles to which each vertex participates.
  int tmp_value;
  int nb_triangles = 0;
  for(int v1(0); v1<nb_vertices; ++v1)
  {
    tmp_value = count_triangles_around_vertex(v1);
    Vertex2Prop[v1] = tmp_value;
    nb_triangles += tmp_value;
  }

  // ===============================================================================================
  // Updates the properties of the graph.
  g_prop["nb_triangles"] = nb_triangles;
  // ===============================================================================================
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::onion_decomposition()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];

  v_prop["kcore"].clear();
  std::vector<double>& Vertex2Coreness = v_prop["kcore"];
  Vertex2Coreness.resize(nb_vertices, 0);

  v_prop["od_layer"].clear();
  std::vector<double>& Vertex2Layer = v_prop["od_layer"];
  Vertex2Layer.resize(nb_vertices, 0);
  
  std::vector<double>& Vertex2Degree = v_prop["degrees"];
  // ===============================================================================================

  // Ensures that the adjacency list has been built.
  if(adjacency_list.size() != nb_vertices)
  {
    build_adjacency_list();
  }

  // Extracts the degree (if not already generated).
  if(Vertex2Degree.size() != nb_vertices)
  {
    degrees();
  }

  // Builds two lists (std::vector, std::set) of the degree of the vertices.
  std::vector<int> DegreeVec(nb_vertices);
  std::set<std::pair<int, int> > DegreeSet;
  for(int v(0); v<nb_vertices; ++v)
  {
    DegreeSet.insert(std::make_pair(Vertex2Degree[v], v));
    DegreeVec[v] = Vertex2Degree[v];
  }


  // Determines the coreness and the layer based on the modified algorithm of Batagelj and
  //   Zaversnik by Hébert-Dufresne, Grochow and Allard.
  int v1, v2, d1, d2;
  int layer = 0;
  std::set< std::pair<int, int> > LayerSet;
  std::set< std::pair<int, int> >::iterator m_it;
  while(DegreeSet.size() > 0)
  {
    // Increases the layer id.
    layer += 1;
    // Populates the set containing the vertices belonging to the same layer.
    m_it = DegreeSet.begin();
    d1 = m_it->first;
    // Sets the coreness and the layer the vertices with the same degree.
    while(m_it->first == d1 && m_it != DegreeSet.end())
    {
      // Sets the coreness and the layer.
      v1 = m_it->second;
      Vertex2Coreness[v1] = d1;
      Vertex2Layer[v1] = layer;
      // Looks at the next vertex.
      ++m_it;
    }
    // Adds the vertices of the layer to the set.
    LayerSet.insert(DegreeSet.begin(), m_it);
    // Removes the vertices of the current layer.
    DegreeSet.erase(DegreeSet.begin(), m_it);
    // Modifies the "effective" degree of the neighbors of the vertices in the layer.
    while(LayerSet.size() > 0)
    {
      // Gets information about the next vertex of the layer.
      v1 = LayerSet.begin()->second;
      // Reduces the "effective" degree of its neighbours.
      for(int n(0), nn(adjacency_list[v1].size()); n<nn; ++n)
      {
        // Identifies the neighbor.
        v2 = adjacency_list[v1][n];
        d2 = DegreeVec[v2];
        // Finds the neighbor in the list "effective" degrees.
        m_it = DegreeSet.find(std::make_pair(d2, v2));
        if(m_it != DegreeSet.end())
        {
          if(d2 > d1)
          {
            DegreeVec[v2] = d2 - 1;
            DegreeSet.erase(m_it);
            DegreeSet.insert(std::make_pair(d2 - 1, v2));
          }
        }
      }
      // Removes the vertices from the LayerSet.
      LayerSet.erase(LayerSet.begin());
    }
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::save_number_edges_in_classes(std::string filename, std::string prop, bool header, int width)
{
  // Stream objects.
  std::fstream output_file;
  // Opens the stream and terminates if the operation did not succeed.
  output_file.open(filename.c_str(), std::fstream::out);
  if( !output_file.is_open() )
  {
    std::cerr << "ERROR: Could not open file: " << filename << "." << std::endl;
    std::terminate();
  }


  // ===============================================================================================
  // Initializes relevant objects of the class.
  std::map< std::pair<int, int>, int>& PropCount = e_class_count_1p[prop];
  // ===============================================================================================


  if(PropCount.size() == 0)
  {
    std::cerr << "ERROR: the number of vertices in classes <" << prop << "> is not available." << std::endl;
    std::terminate();
  }

  // Prints the header (if required).
  if(header)
  {
    output_file << "#" << std::setw(width - 2) << v_prop_header[prop] << "1 ";
    output_file << std::setw(width - 1) << v_prop_header[prop] << "2 ";
    output_file << std::setw(width) << "Nb. edges" << " ";
    output_file << std::endl;
  }

  // Prints the properties.
  std::map< std::pair<int, int>, int>::iterator it = PropCount.begin();
  std::map< std::pair<int, int>, int>::iterator end = PropCount.end();
  for(; it!=end; ++it)
  {
    output_file << std::setw(width) << it->first.first << " ";
    output_file << std::setw(width) << it->first.second << " ";
    output_file << std::setw(width) << it->second << " ";
    output_file << std::endl;
  }

  // Closes the stream.
  output_file.close();
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::save_number_edges_in_classes(std::string filename, std::pair<std::string, std::string> props, bool header, int width)
{
  // Stream objects.
  std::fstream output_file;
  // Opens the stream and terminates if the operation did not succeed.
  output_file.open(filename.c_str(), std::fstream::out);
  if( !output_file.is_open() )
  {
    std::cerr << "ERROR: Could not open file: " << filename << "." << std::endl;
    std::terminate();
  }


  // ===============================================================================================
  // Initializes relevant objects of the class.
  std::map< std::multiset< std::pair<int, int> >, int>& PropCount = e_class_count_2p[props];
  // ===============================================================================================


  if(PropCount.size() == 0)
  {
    std::cerr << "ERROR: the number of vertices in classes <" << props.first << "," << props.second << "> is not available." << std::endl;
    std::terminate();
  }

  // Prints the header (if required).
  if(header)
  {
    output_file << "#" << std::setw(width - 2) << v_prop_header[props.first] << "1 ";
    output_file << std::setw(width - 1) << v_prop_header[props.second] << "1 ";
    output_file << std::setw(width - 1) << v_prop_header[props.first] << "2 ";
    output_file << std::setw(width - 1) << v_prop_header[props.second] << "2 ";
    output_file << std::setw(width) << "Nb. edges" << " ";
    output_file << std::endl;
  }

  // Prints the properties.
  std::map< std::multiset< std::pair< int, int > >, int>::iterator it = PropCount.begin();
  std::map< std::multiset< std::pair< int, int > >, int>::iterator end = PropCount.end();
  for(; it!=end; ++it)
  {
    output_file << std::setw(width) << (it->first).begin()->first << " ";
    output_file << std::setw(width) << (it->first).begin()->second << " ";
    output_file << std::setw(width) << (++(it->first.begin()))->first << " ";
    output_file << std::setw(width) << (++(it->first.begin()))->second << " ";
    output_file << std::setw(width) << it->second << " ";
    output_file << std::endl;
  }

  // Closes the stream.
  output_file.close();
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::save_number_vertices_in_classes(std::string filename, std::string prop, bool header, int width)
{
  // Stream objects.
  std::fstream output_file;
  // Opens the stream and terminates if the operation did not succeed.
  output_file.open(filename.c_str(), std::fstream::out);
  if( !output_file.is_open() )
  {
    std::cerr << "ERROR: Could not open file: " << filename << "." << std::endl;
    std::terminate();
  }


  // ===============================================================================================
  // Initializes relevant objects of the class.
  std::map<int, int>& PropCount = v_class_count_1p[prop];
  // ===============================================================================================


  if(PropCount.size() == 0)
  {
    std::cerr << "ERROR: the number of vertices in classes <" << prop << "> is not available." << std::endl;
    std::terminate();
  }

  // Prints the header (if required).
  if(header)
  {
    output_file << "#" << std::setw(width - 1) << v_prop_header[prop] << " ";
    output_file << std::setw(width) << "Nb. vertices" << " ";
    output_file << std::endl;
  }

  // Prints the properties.
  std::map<int, int>::iterator it = PropCount.begin();
  std::map<int, int>::iterator end = PropCount.end();
  for(; it!=end; ++it)
  {
    output_file << std::setw(width) << it->first << " ";
    output_file << std::setw(width) << it->second << " ";
    output_file << std::endl;
  }

  // Closes the stream.
  output_file.close();
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::save_number_vertices_in_classes(std::string filename, std::pair<std::string, std::string> props, bool header, int width)
{
  // Stream objects.
  std::fstream output_file;
  // Opens the stream and terminates if the operation did not succeed.
  output_file.open(filename.c_str(), std::fstream::out);
  if( !output_file.is_open() )
  {
    std::cerr << "ERROR: Could not open file: " << filename << "." << std::endl;
    std::terminate();
  }


  // ===============================================================================================
  // Initializes relevant objects of the class.
  std::map< std::pair<int, int>, int>& PropCount = v_class_count_2p[props];
  // ===============================================================================================


  if(PropCount.size() == 0)
  {
    std::cerr << "ERROR: the number of vertices in classes <" << props.first << "," << props.second << "> is not available." << std::endl;
    std::terminate();
  }

  // Prints the header (if required).
  if(header)
  {
    output_file << "#" << std::setw(width - 1) << v_prop_header[props.first] << " ";
    output_file << std::setw(width) << v_prop_header[props.second] << " ";
    output_file << std::setw(width) << "Nb. vertices" << " ";
    output_file << std::endl;
  }

  // Prints the properties.
  std::map< std::pair<int, int>, int>::iterator it = PropCount.begin();
  std::map< std::pair<int, int>, int>::iterator end = PropCount.end();
  for(; it!=end; ++it)
  {
    output_file << std::setw(width) << it->first.first << " ";
    output_file << std::setw(width) << it->first.second << " ";
    output_file << std::setw(width) << it->second << " ";
    output_file << std::endl;
  }

  // Closes the stream.
  output_file.close();
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::save_vertices_properties(std::string filename, std::vector<std::string> props_id, vID_t vID, int width, bool header)
{
  // Stream objects.
  std::fstream output_file;
  // Opens the stream and terminates if the operation did not succeed.
  output_file.open(filename.c_str(), std::fstream::out);
  if( !output_file.is_open() )
  {
    std::cerr << "ERROR: Could not open file: " << filename << "." << std::endl;
    std::terminate();
  }


  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  // ===============================================================================================


  // Number of properties.
  int nb_props = props_id.size();

  // Checks whether all properties have been extracted/computed.
  for(int i(0); i<nb_props; ++i)
  {
    is_vertex_property(props_id[i]);
    if(v_prop[props_id[i]].size() != nb_vertices)
    {
      std::cerr << "ERROR: The property " << props_id[i] << " has not been extracted/computed." << std::endl;
      std::terminate();
    }
  }

  // Ensures that the ID2Name has been built. Otherwise, build it.
  if(ID2Name.size() != nb_vertices)
  {
    build_ID2Name();
  }


  // Prints the header (if required).
  if(header)
  {
    int i = 0;
    if(vID == vID_name || vID == vID_num)
    {
      output_file << "#" << std::setw(width - 1) << "Vertex" << " ";
      output_file << std::setw(width) << v_prop_header[props_id[i]] << " ";
      ++i;
    }
    else if(vID == vID_none)
    {
      output_file << "#" << std::setw(width - 1) << v_prop_header[props_id[i]] << " ";
      ++i;
    }
    else
    {
      std::cerr << "ERROR: Unknown vertex identifier type." << std::endl;
      std::terminate();
    }
    for(; i<nb_props; ++i)
    {
      output_file << std::setw(width) << v_prop_header[props_id[i]] << " ";
    }
    output_file << std::endl;
  }

  // Prints the properties.
  for(int v(0); v<nb_vertices; ++v)
  {
    if(vID == vID_name)
    {
      output_file << std::setw(width) << ID2Name[v] << " ";
    }
    else if(vID == vID_num)
    {
      output_file << std::setw(width) << v << " ";
    }
    else if(vID == vID_none) { }
    else
    {
      std::cerr << "ERROR: Unknown vertex identifier type." << std::endl;
      std::terminate();
    }
    for(int i(0); i<nb_props; ++i)
    {
      output_file << std::setw(width) << v_prop[props_id[i]][v] << " ";
    }
    output_file << std::endl;
  }

  // Closes the stream.
  output_file.close();
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
int undirected_graph_t::shortest_path_length(int v1, int v2)
{
  // Ensures that the shortest path lengths have been computed.
  if(shortest_path_lengths.size() != g_prop["nb_vertices"])
  {
    survey_shortest_path_lengths();
  }
  if(v2 > v1)
  {
    return shortest_path_lengths[v2][v1];
  }
  else
  {
    return shortest_path_lengths[v1][v2];
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::survey_connected_components()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];

  v_prop["conn_comp"].clear();
  std::vector<double>& Vertex2Prop = v_prop["conn_comp"];
  Vertex2Prop.resize(nb_vertices, 0);
  // ===============================================================================================

  // Ensures that the adjacency list has been built.
  if(adjacency_list.size() != nb_vertices)
  {
    build_adjacency_list();
  }

  // Starts with every vertex as an isolated cluster.
  std::vector<int> clust_id(nb_vertices);
  std::vector<int> clust_size(nb_vertices, 1);
  for(int v(0); v<nb_vertices; ++v)
  {
    clust_id[v] = v;
  }
  // Merges clusters until the minimal set is obtained.
  merge_clusters(clust_size, clust_id, nb_vertices);
  clust_size.clear();
  // Identifies the connected component to which each vertex belongs.
  int nb_conn_comp = 0;
  int comp_id;
  std::map<int, int> CompID;
  for(int v(0); v<nb_vertices; ++v)
  {
    comp_id = get_root(v, clust_id);
    if(CompID.find(comp_id) == CompID.end())
    {
      CompID[comp_id] = nb_conn_comp;
      ++nb_conn_comp;
    }
    Vertex2Prop[v] = CompID[comp_id];
  }

  // ===============================================================================================
  // Updates the properties of the graph.
  g_prop["nb_conn_comp"] = nb_conn_comp;
  // ===============================================================================================
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::survey_shortest_path_lengths()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  std::vector<double>& Vertex2Degree = v_prop["degrees"];
  shortest_path_lengths.clear();
  // ===============================================================================================


  // Ensures that the adjacency list has been built.
  if(adjacency_list.size() != nb_vertices)
  {
    build_adjacency_list();
  }

  // Extracts the degree (if not already generated).
  if(Vertex2Degree.size() != nb_vertices)
  {
    degrees();
  }

  // Resizes the container for the shortest path lengths.
  shortest_path_lengths.resize(nb_vertices);
  for(int v(0); v<nb_vertices; ++v)
  {
    shortest_path_lengths[v].resize(v+1, nb_vertices);
  }

  // Finds the length of shortest paths via a breadth first search.
  int v1, v2;
  std::vector<int> to_visit;
  to_visit.resize(nb_vertices);
  int next_current, last_current, last_next, current_distance;
  bool keep_going;
  std::vector<bool> is_new(nb_vertices);
  for(int v(0); v<nb_vertices; ++v)
  {
    // No vertex has been visited yet, except "v".
    std::fill(is_new.begin(), is_new.end(), true);
    is_new[v] = false;
    to_visit[0] = v;

    next_current = 0;
    last_current = 0;
    last_next = 0;
    current_distance = 0;
    keep_going = true;
    while(keep_going)
    {
      // Sets the shortest path distance lengths for vertices in the current layer and populates the
      //   next layer.
      for(; next_current <= last_current; ++next_current)
      {
        // Identifies the vertex that has just been reached.
        v1 = to_visit[next_current];
  
        // Sets the shortest path length.
        if(v > v1)
          shortest_path_lengths[v][v1] = current_distance;
        else
          shortest_path_lengths[v1][v] = current_distance;
  
        // Loops over the neighbors of the vertex that has just been reached.
        for(int s(0), ss(Vertex2Degree[v1]); s<ss; ++s)
        {
          // Identifies the neighbor.
          v2 = adjacency_list[v1][s];
          // Adds the vertex to the next layer if it has not been reached yet.
          if(is_new[v2])
          {
            is_new[v2] = false;
            ++last_next;
            to_visit[last_next] = v2;
          }
        }
      }

      // Checks if all vertices have been reached.
      if(last_current == last_next)
      {
        keep_going = false;
      }
      else
      {
        last_current = last_next;
        ++current_distance;
      }
    }
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::survey_triangles()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  triangles.clear();
  // ===============================================================================================


  // Variables.
  int v2, v3, d1, d2;
  // Vector objects.
  std::vector<int> intersection;
  // Set objects.
  std::set<int> neighbours_v1, neighbours_v2;
  // Iterator objects.
  std::vector<int>::iterator it;
  // Newly found triangle.
  std::vector<int> triangle(3);

  // Ensures that the adjacency list has been built.
  if(adjacency_list.size() != nb_vertices)
  {
    build_adjacency_list();
  }

  // Finds all the triangles.
  for(int v1(0); v1<nb_vertices; ++v1)
  {
    // Degree of vertex v1.
    d1 = adjacency_list[v1].size();
    // Performs the calculation only if d1>1.
    if( d1 > 1 )
    {
      // Builds an ordered list of the neighbourhood of v1
      neighbours_v1.clear();
      neighbours_v1.insert(adjacency_list[v1].begin(), adjacency_list[v1].end());
      // Loops over the neighbours of vertex v1.
      for(int n1(0); n1<d1; ++n1)
      {
        // Identity and degree of vertex 2.
        v2 = adjacency_list[v1][n1];
        d2 = adjacency_list[v2].size();
        // Performs the calculation only if d2>1 and if v2>v1 (ensures that each triangle is counted once).
        if( v1 < v2 && d2 > 1 )
        {        
          // Builds an ordered list of the neighbourhood of v2
          neighbours_v2.clear();
          for(int n2(0); n2<d2; ++n2)
          {
            // Identifies the neighbors.
            v3 = adjacency_list[v2][n2];
            if(v2 < v3) // Ensures that triangles will be counted only once.
            {
              neighbours_v2.insert(v3);
            }
          }
          // Identifies the triangles.
          d2 = neighbours_v2.size();
          intersection.clear();
          intersection.resize(std::min(d1, d2));
          it = std::set_intersection(neighbours_v1.begin(), neighbours_v1.end(), neighbours_v2.begin(), neighbours_v2.end(), intersection.begin());
          intersection.resize(it-intersection.begin());
          // Loops over the common neighbours of vertices v1 and v2.
          for(int n(0), nn(intersection.size()); n<nn; ++n)
          {
            // Adds the triangle to the list.
            triangle[0] = v1;
            triangle[1] = v2;
            triangle[2] = intersection[n];
            triangles.push_back(triangle);
          }
        }
      }
    }
  }

  // ===============================================================================================
  // Updates the properties of the graph.
  g_prop["nb_triangles"] = triangles.size();
  // ===============================================================================================
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::topological_dimensions()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  int nb_conn_comp = g_prop["nb_conn_comp"];
  std::vector<double>& Vertex2Comp = v_prop["conn_comp"];
  // ===============================================================================================


  // Ensures that the shortest path lengths have been computed.
  if(shortest_path_lengths.size() != nb_vertices)
  {
    survey_shortest_path_lengths();
  }

  // Ensures that the connected components have been identified.
  if(Vertex2Comp.size() != nb_vertices)
  {
    survey_connected_components();
    nb_conn_comp = g_prop["nb_conn_comp"];
  }

  // Resizes the containers.
  average_shortest_path_lengths.resize(nb_conn_comp, 0);
  diameters.resize(nb_conn_comp, 0);

  // Computes the average length.
  int c;
  double l;
  std::vector<int> comp_size(nb_conn_comp, 0);
  double avg_value = 0;
  for(int i(0), ii(shortest_path_lengths.size()); i<ii; ++i)
  {
    c = Vertex2Comp[i];
    for(int j(0), jj(shortest_path_lengths[i].size()-1); j<jj; ++j)
    {
      l = shortest_path_lengths[i][j];
      if(l < nb_vertices)
      {
        comp_size[c] += 1;
        average_shortest_path_lengths[c] += l;
        if(l > diameters[c])
        {
          diameters[c] = l;
        }
      }
    }
  }
  // Completes the calculation of the average values.
  for(int s(0); s<nb_conn_comp; ++s)
  {
    average_shortest_path_lengths[s] /= comp_size[s];
  }
}
