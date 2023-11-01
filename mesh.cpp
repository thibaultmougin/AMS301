#include "headers.hpp"

extern int myRank;
extern int nbTasks;

//================================================================================
// Read the mesh from a gmsh-file (.msh) and store in a mesh-structure 'mesh'
//================================================================================

void readMsh(Mesh& mesh, string fileName)
{
  if(myRank == 0)
    cout << "== read mesh file" << endl;
  
  // Open mesh file
  ifstream meshfile(fileName.c_str());
  if(!meshfile){
    cout << "ERROR: Mesh '" << fileName.c_str() << "' not opened." << endl;
    exit(EXIT_FAILURE);
  }
  
  // Read mesh file
  double dummy;
  while(meshfile.good()){
    
    string line;
    getline(meshfile,line);
    
    // Read NODES
    if(line.compare("$Nodes") == 0){
      
      // Read the number of nodes
      meshfile >> mesh.nbOfNodes;
      
      // Read the node coordinates
      mesh.coords.resize(mesh.nbOfNodes,3);
      for(int nGlo=0; nGlo<mesh.nbOfNodes; nGlo++)
        meshfile >> dummy >> mesh.coords(nGlo,0) >> mesh.coords(nGlo,1) >> mesh.coords(nGlo,2);
      
      if(myRank == 0){
        cout << "   -> " << mesh.nbOfNodes << " nodes" << endl;
      }
    }
    
    // Read ELEMENTS
    if(line.compare("$Elements") == 0){
      
      // Save the total number of triangles
      meshfile >> mesh.nbOfTri;
      if(myRank == 0)
        cout << "   -> " <<  mesh.nbOfTri << " triangles" << endl;
      if(mesh.nbOfTri < 1){
        cout << "ERROR: No elements (mesh.nbOfTri = " << mesh.nbOfTri << ")." << endl;
        exit(EXIT_FAILURE);
      }
      
      // Resize arrays for elements infos/nodes for each type of element
      mesh.triNum.resize(mesh.nbOfTri);
      mesh.triNodes.resize(mesh.nbOfTri,3);
      mesh.triPart.resize(mesh.nbOfTri);
      
      // Temporary arrays for elements infos/nodes
      IntMatrix elemNodes(mesh.nbOfTri,3);      // associated nodes (max 3. for triangle)
      IntVector elemNum(mesh.nbOfTri);          // gmsh numerbing
      IntVector elemPart(mesh.nbOfTri);         // partition tag
      
      // Read infos/nodes for all the elements
      for(int iTri=0; iTri<mesh.nbOfTri; iTri++){
        getline(meshfile, line);
        
        // Save gmsh number for each triangle
        meshfile >> mesh.triNum(iTri);
        
        // Check the type of element
        int elemType;
        meshfile >> elemType;
        if (elemType != 2){
          cout << "ERROR: Element type '" << elemType << "' not supported." << endl;
          exit(EXIT_FAILURE);
          break;
        }
        
        // Save partition number for each triangle
        int infos;
        meshfile >> infos >> dummy >> dummy;
        int elemPart = 1;
        if(infos > 2)
          meshfile >> dummy >> elemPart;    // partition tag
        for(int j=5; j<=infos; j++)
          meshfile >> dummy;                // useless infos
        mesh.triPart(iTri) = (elemPart-1) % nbTasks; // (gmsh is 1-index, here is 0-index)
        
        // Save nodes for each triangle
        IntVector elemNodes(3);
        meshfile >> elemNodes(0);
        meshfile >> elemNodes(1);
        meshfile >> elemNodes(2);
        mesh.triNodes(iTri,0) = elemNodes(0)-1; // (gmsh is 1-index, here is 0-index)
        mesh.triNodes(iTri,1) = elemNodes(1)-1;
        mesh.triNodes(iTri,2) = elemNodes(2)-1;
      }
    }
  }
  meshfile.close();
}

//================================================================================
// Export a field contained in 'vec' in a gmsh-file (.msh)
//================================================================================

void exportFieldMsh(ScaVector& vec, Mesh& mesh, string viewName, string fileName)
{
  if(nbTasks > 1){
    ostringstream ss;
    ss << fileName << "_" << myRank;
    fileName = ss.str();
  }

  ofstream posFile(fileName.c_str());
  posFile << "$MeshFormat" << endl;
  posFile << "2.2 0 0" << endl;
  posFile << "$EndMeshFormat" << endl;
  
  posFile << "$ElementNodeData" << endl;
  posFile << "2" << endl;
  posFile << "\"" << viewName.c_str() << "\"" << endl;  // name of the view
  posFile << "" << endl;
  posFile << "1" << endl;
  posFile << "0" << endl; // ("Time")
  posFile << "4" << endl;
  posFile << "0" << endl; // ("timeStep")
  posFile << "1" << endl; // ("numComp")
  posFile << mesh.nbOfTri << endl;   // total number of elementNodeData in this file
  posFile << myRank << endl;
  for(int iTriLoc=0; iTriLoc<mesh.nbOfTri; iTriLoc++){
    posFile << mesh.triNum[iTriLoc] << " " << 3;
    for(int n=0; n<3; n++){
      int nLoc = mesh.triNodes(iTriLoc,n);
      posFile << " " << vec(nLoc);
    }
    posFile << endl;
  }
  posFile << "$EndElementNodeData" << endl;
  posFile.close();
}