#ifndef READJSON_H
#define READJSON_H
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

int read_cmd_by_json(const int &argc, const char* argv[], boost::property_tree::ptree &pt);

#endif
