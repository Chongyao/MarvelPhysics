#include <src/readjson.h>
#include <iostream>
using namespace std;
using namespace boost::property_tree;

int read_cmd_by_json(const int &argc,  char* argv[],  boost::property_tree::ptree &pt){

  const string jsonfile_path = argv[1];
  cout << jsonfile_path << endl;
  const size_t ext = jsonfile_path.rfind(".json");
  if (ext != std::string::npos){
    read_json(jsonfile_path, pt);
    return 0;    
  }

  
  else
    return -1;
  
  
}
