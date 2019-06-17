#ifndef VECTOR_IO_H
#define VECTOR_IO_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <iostream>

template<typename INT, typename DOUBLE> 
struct vec {
  INT len;
  DOUBLE *val;
};


template<class INT, class DOUBLE>
vec<INT, DOUBLE> new_vec(const int len)
{
  vec<INT, DOUBLE> v;
  v.len = len;
  v.val = (DOUBLE*)malloc(sizeof(DOUBLE)*len);
  return v;
}

template<class INT, class DOUBLE>
vec<INT, DOUBLE> print_vec(vec<INT, DOUBLE> v)
{
  for (int i=0; i<v.len; ++i)
  {
    std::cout << v.val[i] << std::endl;
  }
}

template<class INT, class DOUBLE>
DOUBLE l2_norm(vec<INT, DOUBLE> v1,
                  vec<INT, DOUBLE> v2)
{
  if (v1.len != v2.len)
  {
    std::cout << "Different vector length make no sense." << std::endl;
    exit(1);
  }

  DOUBLE norm = 0.0;

  for (int i=0; i<v1.len; ++i)
  {
    DOUBLE d = v1.val[i] - v2.val[i];
    norm = norm + (d*d);
  }

  return sqrt(norm);
}

template<class INT, class DOUBLE>
vec<INT, DOUBLE> read_vec(const char *filename)
{
  FILE *fd = fopen(filename, "r");
    
  if (fd == NULL){
    std::cout << "Unable to open file: " << filename << std::endl;
    exit(1);
  }
  fclose(fd);
  std::cout << "Reading vector from file: " << filename << std::endl;
  
  INT count = 0;
  std::ifstream ifs(filename);
  std::string line;
  getline(ifs, line);
  std::istringstream iss_num(line);
  int num;
  iss_num >> num;
  std::cout << line << std::endl;
  count = (INT)num;

  std::cout << "Vector length: " << count << std::endl;

  DOUBLE *val = (DOUBLE*)malloc(sizeof(DOUBLE)*count);

  for (int i=0; i<count; ++i)
  {
    getline(ifs, line);
    std::istringstream iss_val(line);
    double v;
    iss_val >> v;
    val[i] = (DOUBLE)v;
  }
  
  vec<INT, DOUBLE> v;
  v.len = count;
  v.val = val;
  return v;
}

template <class INT, class DOUBLE>
void delete_vec(vec<INT, DOUBLE> v)
{
  free(v.val);
}

#endif
