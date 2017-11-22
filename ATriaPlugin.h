#ifndef ATRIAPLUGIN_H
#define ATRIAPLUGIN_H

#include "Plugin.h"
#include "PluginProxy.h"
#include "csv_parser/csv_parser.cpp"
#include <string>
#include <vector>
//#define GSIZE 126
//#define GSIZE 241

class ATriaPlugin : public Plugin 
{
  public:    
  ~ATriaPlugin();
  //std::string toString(){return "ATria";}
  void input(std::string file);
  void run();
  void output(std::string file);

  private:
     float* OrigGraph;
     //std::string bacteria[GSIZE];
     std::string* bacteria;//[GSIZE];
     float* H_G;
     std::vector<float> U;   
     float* H_pay;
     int GSIZE;
     void _CPU_Floyd(float*, int);
};


#endif
