/*
## Copyright (C) <2017>  <bioMerieux, Universite Claude Bernard Lyon 1,
## Centre National de la Recherche Scientifique>

## 1. This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU Affero General Public License as published
## by the Free Software Foundation version 3 of the  License and under the
## terms of article 2 below.
## 2. This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
## or FITNESS FOR A PARTICULAR PURPOSE. See below the GNU Affero General
## Public License for more details.
## You should have received a copy of the GNU Affero General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
## 3. Communication to the public by any means, in particular in the form of
## a scientific paper, a poster, a slideshow, an internet page, or a patent,
## of a result obtained directly or indirectly by running this program must
## cite the following paper :
##  Magali Jaillard, Maud Tournoud, Leandro Lima, Vincent Lacroix,
##  Jean-Baptiste Veyrieras and Laurent Jacob, "Representing Genetic
##  Determinants in Bacterial GWAS with Compacted De Bruijn Graphs", 2017,
##  Cold Spring Harbor Labs Journals, doi:10.1101/113563.
##  (url: http://www.biorxiv.org/content/early/2017/03/03/113563)
## -------------------------------------------------------------------------

## Authors (alphabetically): Jacob L., Jaillard M., Lima L.
*/

#include "generate_output.h"
#include "Utils.h"
#include "global.h"
#include <cstdlib>
#include <ctime>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/copy.hpp>
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS

using namespace std;

//variables regarding the size of the nodes
int minSize = 50;
int maxSize = 400;



// The Tool constructor allows to give a name to our tool.
// This name appears when one gets help in the command line or in the final output
generate_output::generate_output ()  : Tool ("generate_output") //give a name to our tool
{
  populateParser(this);
}

void generateCytoscapeOutput(const graph_t &graph, const vector<int> &nodes, const string &typeOfGraph, int i,
                             const string &outputFolder, const vector<int> &selectedUnitigs, int nbPheno0, int nbPheno1) {
  //gets the maxCoverage of the nodes in this component
  int maxCoverage=-1;
  for (const auto &node : nodes) {
    MyVertex v = vertex(node, graph);
    maxCoverage = max(maxCoverage, graph[v].phenoCounter.getTotal());
  }


  //declares the stringstream which will store the elements to be printed
  stringstream elementsSS;
  elementsSS << scientific;

  //goes through the nodes and print them
  set<MyVertex> verticesInThisComponent;
  for (const auto &node : nodes) {
    MyVertex v = vertex(node, graph);
    verticesInThisComponent.insert(v);

    //TODO: if the graph file is too large, a good idea is to try to remove the name (node sequence)
    elementsSS << "{data: {id: 'n" << graph[v].id << "', name: '" << graph[v].name << "'" <<
        ", info: '" << graph[v].id << "'" <<
        ", total: '" << graph[v].phenoCounter.getTotal() << "'" <<
        ", pheno0: '" << graph[v].phenoCounter.getPheno0() << "/" << nbPheno0 << "'" <<
        ", pheno1: '" << graph[v].phenoCounter.getPheno1() << "/" << nbPheno1 << "'" <<
        ", NA: '" << graph[v].phenoCounter.getNA() << "'" <<
        ", qValue: '" << graph[v].unitigStats.getQValueAsStr() << "'" <<
        ", weight: '" << graph[v].unitigStats.getWeightAsStr() << "'" <<
        ", waldStatistic: '" << graph[v].unitigStats.getWaldStatisticAsStr() << "'" <<
        "}, style: {'background-color': rgbToHex(" << graph[v].unitigStats.getRGB() << ")" <<
        ", 'width': " << (minSize + (((double)graph[v].phenoCounter.getTotal())/maxCoverage*(maxSize-minSize))) <<
        ", 'height': " << (minSize + (((double)graph[v].phenoCounter.getTotal())/maxCoverage*(maxSize-minSize))) <<
        //if it is not a selected unitig, then it becomes a little bit transparent
        (find(selectedUnitigs.begin(), selectedUnitigs.end(), graph[v].id)==selectedUnitigs.end() ? ", 'background-opacity': 0.3" : "") <<
        "}}, ";
  }

  //goes through the edges and print them
  for (auto ep = edges(graph); ep.first != ep.second; ++ep.first) {
    MyEdge e = *ep.first;
    if (verticesInThisComponent.find(source(e, graph)) != verticesInThisComponent.end() &&
        verticesInThisComponent.find(target(e, graph)) != verticesInThisComponent.end()) {
      elementsSS << "{data: {id: 'e" << graph[e].id << "', source: 'n" << graph[source(e, graph)].id
      << "', target: 'n" << graph[target(e, graph)].id << "'}}, ";
    }
  }
  string elements = elementsSS.str();

  //read template file
  string templatePath = pathToExecParent + string("/cytoscape_template.html");
  string cytoscapeOutput = readFileAsString(templatePath.c_str());

  //put the info in the template file
  cytoscapeOutput.replace(cytoscapeOutput.find("<elementsTag>"), ((string) "<elementsTag>").size(), elements);

  //output the file
  string outfilename;
  {
    stringstream ss;
    ss << outputFolder << "/visualisations/components/" << typeOfGraph << "_" << i << ".html";
    outfilename=ss.str();
  }
  ofstream outFile;
  openFileForWriting(outfilename, outFile);
  outFile << cytoscapeOutput;
  outFile.close();

  //copy the lib folder, if it is not already copied
  string fromLibPath = pathToExecParent + string("/lib");
  string toLibPath = outputFolder + string("/visualisations/components/lib");
  if (!boost::filesystem::exists(toLibPath))
    copyDirectoryRecursively(fromLibPath, toLibPath);
}

void generate_output::execute () {
  checkParametersGenerateOutput(this);
  int neighbourhood = getInput()->getInt(STR_MAX_NEIGHBOURHOOD);
  //string outputFolder = getInput()->getStr(STR_OUTPUT);
  string outputFolder("output");
  //TODO: use several cores here?
  //int nbCores = getInput()->getInt(STR_NBCORES);

  //get the nbContigs
  int nbContigs = getNbLinesInFile(outputFolder+string("/graph.nodes"));



  cerr << "[Getting the significant unitigs from the patterns...]" << endl;

  //read the patterns
  auto sortedPatterns = PatternFromStats::readFile(outputFolder + "/patterns.txt", true);

  //associates each unitig to its pattern stats
  vector<const PatternFromStats*> unitigToPatternStats(nbContigs, NULL);
  {
    //since not all patterns have stats, this maps the id of the pattern to its stats, so we know easily which pattern has stats
    map<int, const PatternFromStats *> patternsToStats;
    for (const auto &pattern : sortedPatterns)
      patternsToStats[pattern.pattern] = &pattern;

    //associates each unitig to its Pattern stats
    ifstream unitigToPatternFilestream;
    openFileForReading(outputFolder + string("/gemma_input.unitig_to_pattern.binary"), unitigToPatternFilestream);

    int unitig, pattern;
    while (unitigToPatternFilestream >> unitig >> pattern) {
      if (patternsToStats.find(pattern) != patternsToStats.end()) {
        //the pattern has a stats
        unitigToPatternStats[unitig] = patternsToStats[pattern];
      }
    }
    unitigToPatternFilestream.close();
  }

  //read unitigToWeightCorrection
  vector<int> unitigToWeightCorrection(nbContigs);
  ifstream unitigToWeightCorrectionFilestream;
  openFileForReading(outputFolder+string("/weight_correction"), unitigToWeightCorrectionFilestream);
  for (int i=0;i<unitigToWeightCorrection.size();i++)
    unitigToWeightCorrectionFilestream >> unitigToWeightCorrection[i];
  unitigToWeightCorrectionFilestream.close();

  //get the significant patterns according to SFF
  vector<PatternFromStats> significantPatterns;
  {
    boost::apply_visitor(GetSignificantPatterns(sortedPatterns, significantPatterns), SFF);
    cerr << "Selected " << significantPatterns.size() << " most significant patterns..." << endl;
  }

  //transform the significant patterns into significant unitigs
  {
    ofstream significantUnitigsFile;
    openFileForWriting(outputFolder + string("/significant_unitigs.txt"), significantUnitigsFile);
    ifstream pattern2UnitigsFile;
    openFileForReading(outputFolder + string("/bugwas_input.unique_rows_to_all_rows.binary"), pattern2UnitigsFile);

    //read pattern2UnitigsFile
    string line;
    vector<vector<int> > pattern2Unitigs;
    while (getline(pattern2UnitigsFile, line)) {
      stringstream ss(line);
      int unitig;
      vector<int> unitigs;
      while (ss >> unitig)
        unitigs.push_back(unitig);
      if (unitigs.size() > 0)
        pattern2Unitigs.push_back(unitigs);
    }

    //process significantPatternsFile and produces significantUnitigsFile
    for (const auto &significantPattern : significantPatterns) {
      for (const auto &unitig : pattern2Unitigs[significantPattern.pattern]) {
        significantUnitigsFile << unitig << endl;
      }
    }

    pattern2UnitigsFile.close();
    significantUnitigsFile.close();
  }
  cerr << "[Getting significant unitigs from the patterns...] - Done!" << endl;

  cerr << "[Reading input and creating BOOST graph...]" << endl;

  //declare the graph
  graph_t graph(nbContigs);

  //get the unitigs selected for visualization
  vector<int> selectedUnitigs;
  {
    ifstream selectedUnitigsStream;
    openFileForReading(outputFolder + string("/significant_unitigs.txt"), selectedUnitigsStream);
    int selectedUnitig;
    while (selectedUnitigsStream >> selectedUnitig)
      selectedUnitigs.push_back(selectedUnitig);
    selectedUnitigsStream.close();
  }

  //read the frequency file
  vector< PhenoCounter > phenotypeCounters;
  {
    string frequencyFilename(outputFolder+string("/frequency_unitig_to_total_pheno0_pheno1_NA_count"));
    ifstream frequencyFile;
    openFileForReading(frequencyFilename, frequencyFile);

    int temp1, temp2, temp3, temp4;
    while (frequencyFile >> temp1) {
      frequencyFile >> temp2 >> temp3 >> temp4;
      phenotypeCounters.push_back(PhenoCounter(temp2, temp3, temp4));
    }

    frequencyFile.close();
  }

  //read the total number of pheno0 and pheno1 strains
  int nbPheno0, nbPheno1;
  {
    ifstream totalNbOfStrainsInEachPheno;
    openFileForReading(outputFolder + string("/total_nb_of_strains_in_each_pheno"), totalNbOfStrainsInEachPheno);
    totalNbOfStrainsInEachPheno >> nbPheno0 >> nbPheno1;
    totalNbOfStrainsInEachPheno.close();
  }




  //create the nodes of the graph
  //assume the graph is undirected, instead of directed
  string nodesFile = outputFolder+string("/graph.nodes");
  {
    ifstream nodesFileReader;
    openFileForReading(nodesFile, nodesFileReader);
    int id;
    string seq;
    int index = 0;
    while (nodesFileReader >> id >> seq) {
      MyVertex vF = vertex(index, graph);
      //TODO: sequence should be stored as a 2-bit array
      graph[vF].name = seq;
      graph[vF].id = id;
      graph[vF].strand = 'F';
      graph[vF].phenoCounter = phenotypeCounters[id];
      graph[vF].unitigStats = UnitigStats(unitigToPatternStats[id], unitigToWeightCorrection[id]);
      index++;
    }
    nodesFileReader.close();
  }

  //create the edges of the graph
  string edgesFile = outputFolder+string("/graph.edges.dbg");
  {
    ifstream edgesFileReader;
    openFileForReading(edgesFile, edgesFileReader);
    int from, to;
    string label;
    int index = 0;
    while (edgesFileReader >> from >> to >> label) {
        pair<MyEdge, bool> return_from_add_edge = add_edge(vertex(from, graph),
                                                           vertex(to, graph),
                                                           graph);
        if (return_from_add_edge.second) {
          graph[return_from_add_edge.first].id = index;
          graph[return_from_add_edge.first].weight = 1;
          index++;
        }
    }
    edgesFileReader.close();
  }
  cerr << "[Reading input and creating BOOST graph...] - Done!" << endl;
  //DEBUG - print the graph to a file
  //GraphWriter::writeGraphToFile(graph, outputFolder+string("/graph.debug.cytoscape");
  //DEBUG - print the graph to a file


  cerr << "[Computing nodes' neighbourhoods...]" << endl;
  //compute the neighbourhood
  vector<int> distances(num_vertices(graph));
  set<MyVertex> verticesInTheNeighbourhood;
  map <int, set<MyVertex> > verticesInTheNeighbourhoodOfThisUnitig;

  for (auto unitig : selectedUnitigs) {
    MyVertex selectedVertex = vertex(unitig, graph);
    try {
      fill(distances.begin(), distances.end(), 0);
      dijkstra_shortest_paths(graph, selectedVertex, weight_map(get(&EdgeInfo::weight, graph)).
          distance_map(make_iterator_property_map(distances.begin(),
                                                  boost::get(boost::vertex_index, graph))).
          visitor(StopWhenVeryDistantFromSourceDijkstraVisitor(distances, neighbourhood,
                                                               verticesInTheNeighbourhoodOfThisUnitig[unitig])));
    } catch (TooDistant &e) { }

    //debug
    /*
    cout << "Vertices in the neighbourhood of " << unitig << endl;
    for (auto v : verticesInTheNeighbourhoodOfThisUnitig)
      cout << v << ",";
    cout << endl;
     */
    //debug

    //update verticesInTheNeighbourhood
    verticesInTheNeighbourhood.insert(verticesInTheNeighbourhoodOfThisUnitig[unitig].begin(), verticesInTheNeighbourhoodOfThisUnitig[unitig].end());
  }
  cerr << "[Computing nodes' neighbourhoods...] - Done!" << endl;

  //debug
  /*
  cout << "Vertices in the neighbourhood: " << endl;
  for (auto v : verticesInTheNeighbourhood)
    cout << v << " ";
  cout << endl;
   */

  cerr << "[Generating the visualisation files...]" << endl;
  //print #1:
  //print the whole graph
/*
  {
    vector<int> vectorVertex(verticesInTheNeighbourhood.begin(), verticesInTheNeighbourhood.end());
    generateCytoscapeOutput(*filteredGraph, graph, vectorVertex, "whole", 0);
  }
*/
  //print #2:
  //print one graph per component
  //create a subgraph containing only the nodes in the neighbourhoods
  graph_t& newGraph = graph.create_subgraph();
  vector<vector<int> > nodesInComponent; //care: the nodes in this variable are the nodes in newGraph, not in the normal graph
  int numberOfComponents=0;
  {
    for (auto vp = vertices(graph); vp.first != vp.second; ++vp.first) {
      MyVertex v = *vp.first;
      if (verticesInTheNeighbourhood.find(v)!=verticesInTheNeighbourhood.end()) {
        add_vertex(v, newGraph);
      }
    }

    //get the components of this subgraphs
    vector<int> componentOfThisNode = vector<int>(num_vertices(newGraph));
    int num = connected_components(newGraph, &componentOfThisNode[0]);
    nodesInComponent = vector<vector<int> >(num);
    for (int i = 0; i != componentOfThisNode.size(); ++i) {
      nodesInComponent[componentOfThisNode[i]].push_back(i);
    }

    for (int i = 0; i < nodesInComponent.size(); i++) {
      generateCytoscapeOutput(newGraph, nodesInComponent[i], "comp", i, outputFolder, selectedUnitigs, nbPheno0, nbPheno1);
    }
    numberOfComponents = nodesInComponent.size();
  }

  //print #3:
  //print one graph per selected unitig
  /*
  {
    for (auto unitig : selectedUnitigs) {
      vector<int> vectorVertex(verticesInTheNeighbourhoodOfThisUnitig[unitig].begin(), verticesInTheNeighbourhoodOfThisUnitig[unitig].end());
      generateCytoscapeOutput(*filteredGraph, graph, vectorVertex, "unitig", unitig);
    }
  }*/
  cerr << "[Generating the visualisation files...] - Done!" << endl;


  //create the index
  cerr << "[Creating index file...]" << endl;
  //create the thumbnails
  for (int i=0;i<numberOfComponents;i++) {
    string HTMLFile(boost::filesystem::canonical(outputFolder+"/visualisations/components/comp_"+std::to_string(i)+".html").string());
    string PNGFile = HTMLFile+".png";
    executeCommand("./phantomjs render_graph.js " + HTMLFile + " " + PNGFile, false);
  }

  //create the index
  //create the object previews for each component
  vector<ObjectPreview> previews;
  string preview="<table style=\"float: left; margin-left: 10px; margin-bottom: 10px\">"
      "        <tr> <td>"
      "          <center><strong>Comp_<id></strong></center>"
      "        </td> </tr>"
      "        <tr> <td>"
      "          q-value: <q-value>"
      "        </td> </tr>"
      "        <tr> <td>"
      "          Annotations: <annotations>"
      "        </td> </tr>"
      "        <tr> <td>"
      "          Preview (click to open) <br/>"
      "          <center>"
      "          <img id=\"img_<id>\" onclick=\"window.open(\\'components/comp_<id>.html\\',\\'_blank\\')\" src=\"components/comp_<id>.html.png\"/>"
      "          </center>"
      "        </td> </tr>"
      "      </table>";
  for (int i=0;i<numberOfComponents;i++) {
    string idString = std::to_string(i);
    string annotations = string("TODO");
    string thisPreview(preview);

    //get the q-value
    long double smallerQValue = std::numeric_limits<long double>::max();
    for(const auto &node : nodesInComponent[i]) {
      MyVertex v = vertex(node, newGraph);
      int id = newGraph[v].id;

      //check if the patterns exists
      if (unitigToPatternStats[id])
        smallerQValue=min(smallerQValue, unitigToPatternStats[id]->qValue);
    }

    //fix this preview
    boost::replace_all(thisPreview, "<id>", idString);
    boost::replace_all(thisPreview, "<annotations>", annotations);
    string qValueAsStr;
    {
      stringstream ss;
      ss << scientific;
      ss << smallerQValue;
      ss >> qValueAsStr;
    }
    boost::replace_all(thisPreview, "<q-value>", qValueAsStr);

    previews.push_back(ObjectPreview(smallerQValue, annotations, thisPreview));
  }


  //output the object previews
  stringstream ssPreview;
  ssPreview << "[";
  for (const auto &preview : previews)
    ssPreview << preview.toJSObject() << ", ";
  ssPreview << "]";

  //create the index file
  //read template file
  string templatePath = pathToExecParent + string("/index_template.html");
  string indexOutput = readFileAsString(templatePath.c_str());

  //put the info in the template file
  boost::replace_all(indexOutput, "<previews>", ssPreview.str());

  //output the file
  ofstream indexFile;
  openFileForWriting(outputFolder+string("/visualisations/index.html"), indexFile);
  indexFile << indexOutput;
  indexFile.close();
  cerr << "[Creating index file...] - Done!" << endl;


  //tell we are done
  cout << endl << endl <<
      "******************************************************************************" << endl <<
      "We are done. The output can be found at " << outputFolder << "/visualisations/index.html" << endl <<
      "******************************************************************************" << endl << endl;
  cout.flush();
}