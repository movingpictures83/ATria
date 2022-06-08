#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include "PluginManager.h"
#include "ATriaPlugin.h"
#include "PluMA.h"

ATriaPlugin::~ATriaPlugin() {
   if (OrigGraph) free(OrigGraph);
   if (bacteria) delete bacteria;
   if (H_G) delete H_G;
   if (H_pay) delete H_pay;
}

void ATriaPlugin::input(std::string file) {
                const char field_terminator = ',';
                const char line_terminator  = '\n';
                const char enclosure_char   = '"';
   // File is in CSV format
  csv_parser file_parser;
                file_parser.set_skip_lines(1);
                file_parser.init(file.c_str());
                file_parser.set_enclosed_char(enclosure_char, ENCLOSURE_OPTIONAL);
                file_parser.set_field_term_char(field_terminator);
                file_parser.set_line_term_char(line_terminator);

                GSIZE = 0;
                while (file_parser.has_more_rows()) {
                   file_parser.get_row();
                   GSIZE++;
                }
                const int NumBytes=(GSIZE*2)*(GSIZE*2)*sizeof(float);
                OrigGraph=(float *)malloc(NumBytes);//will be original Adjancency matrix, will NOT be changed
                bacteria = new std::string[GSIZE];                               
 
                file_parser.init(file.c_str());
                unsigned int row_count = 0;
                while(file_parser.has_more_rows())
                {
                        unsigned int i = 0;

                        csv_row row = file_parser.get_row();

                        bacteria[row_count] = row[0];

                        for (i = 1; i < row.size(); i++) {
                              int bac1 = row_count;
                              int bac2 = i-1;
                              float weight = atof(row[i].c_str());
                              if (bac1 != bac2) {
                                 if (weight > 0) {
                                    OrigGraph[(bac1*2)*(2*GSIZE)+(bac2*2)] = weight;
                                    OrigGraph[(bac1*2+1)*(2*GSIZE)+(bac2*2+1)] = weight;
                                    OrigGraph[(bac1*2+1)*(2*GSIZE)+(bac2*2)] = 0;
                                    OrigGraph[(bac1*2)*(2*GSIZE)+(bac2*2+1)] = 0;
                                 }
                                 else if (weight < 0) {
                                    OrigGraph[(bac1*2+1)*(2*GSIZE)+(bac2*2)] = weight;
                                    OrigGraph[(bac1*2)*(2*GSIZE)+(bac2*2+1)] = weight;
                                    OrigGraph[(bac1*2)*(2*GSIZE)+(bac2*2)] = 0;
                                    OrigGraph[(bac1*2+1)*(2*GSIZE)+(bac2*2+1)] = 0;
                                 }
                                else {
                                    OrigGraph[(bac1*2)*(2*GSIZE)+(bac2*2)] = 0;
                                    OrigGraph[(bac1*2+1)*(2*GSIZE)+(bac2*2+1)] = 0;
                                    OrigGraph[(bac1*2+1)*(2*GSIZE)+(bac2*2)] = 0;
                                    OrigGraph[(bac1*2)*(2*GSIZE)+(bac2*2+1)] = 0;
                                }
                              }
                              else {
                                    OrigGraph[(bac1*2)*(2*GSIZE)+(bac2*2)] = 1; // Start these at 1, because they are starting verts.
                                    OrigGraph[(bac1*2+1)*(2*GSIZE)+(bac2*2+1)] = 1;
                                    OrigGraph[(bac1*2+1)*(2*GSIZE)+(bac2*2)] = 0;
                                    OrigGraph[(bac1*2)*(2*GSIZE)+(bac2*2+1)] = 0;
                              }
                        }
                        row_count++;
                }
  

}

void ATriaPlugin::run() {
   PluginManager::log("I am running ATria");
const int NumBytes=(GSIZE*2)*(GSIZE*2)*sizeof(float);
   H_G=(float *)malloc(NumBytes);
   U.resize(GSIZE, 0.0);

const int NumBytesPay=(GSIZE)*sizeof(float);   // N by 2N.  2 paths from every i to every j
                H_pay = (float *)malloc(NumBytesPay);

		int currentrank = 1;
for (int a = 0; a < GSIZE; a++) {
                //cout<<"Successfully created random highly connected graph in adjacency Matrix form with "<<RANDOM_GSIZE*RANDOM_GSIZE<< " elements.\n";
                for(int i=0;i<GSIZE*2*GSIZE*2;i++){//copy for use in computation
                        H_G[i]=OrigGraph[i];//copy for use in computation
                }

                //cout<<"\nFloyd-Warshall on CPU underway:\n";

                _CPU_Floyd(H_G,/*H_Gpath,*/GSIZE*2);//find shortest paths (with path construction) on serial CPU (Intel i7 3770 3.9 ghz)

                for (int i = 0; i < GSIZE; i++)
                {
                   float pay = 0;
                   for (int j = 0; j < 2*GSIZE; j++)
                      pay += H_G[i*2*(2*GSIZE)+j];
                   pay--;
                   H_pay[i] = pay;
                }
		vector<int> maxnodes;
                int mnode = -1;
                float maxpay = -1;
                for (int i = 0; i < GSIZE; i++) {
                   //cout << "Pay for " << bacteria[i] << ": " << H_pay[i] << endl;
                   if (fabs(H_pay[i]) > maxpay) {
                      mnode = i;
                      maxpay = fabs(H_pay[i]);
                   }
                }
                if (maxpay == 0)
                   break;
		maxnodes.push_back(mnode);
		for (int i = 0; i < GSIZE; i++) {
                   if ((i != mnode) && fabs(H_pay[i]) == maxpay) {
                      maxnodes.push_back(i);     
	           }
		}
		for (int w = 0; w < maxnodes.size(); w++) {
			int maxnode = maxnodes[w];
                PluginManager::log(std::string("Node with highest pay: "+bacteria[maxnode]+": "+std::to_string(H_pay[maxnode])));
                U[maxnode] = currentrank;//H_pay[maxnode];
                // Non-GPU Triad Removal
                for (int i = 0; i < GSIZE*2; i++) {
                   if ((i/2 != maxnode) &&
                       (OrigGraph[maxnode*2*GSIZE*2+i] != 0 || OrigGraph[(maxnode*2+1)*GSIZE*2+i] != 0)){
                      for (int j = i+1; j < GSIZE*2; j++) {
                         if ((j/2 != maxnode) && (OrigGraph[maxnode*2*GSIZE*2+j] != 0 || OrigGraph[(maxnode*2+1)*GSIZE*2+j] != 0) &&
                             OrigGraph[i*GSIZE*2+j] != 0) {
                            OrigGraph[i*GSIZE*2+j] = 2;
                            OrigGraph[j*GSIZE*2+i] = 2;
                         }
                      }
                      if (OrigGraph[maxnode*2*GSIZE*2+i] != 0) {
                         OrigGraph[maxnode*2*GSIZE*2+i] = 2;
                         OrigGraph[i*GSIZE*2+(maxnode*2)] = 2;
                      }
                      if (OrigGraph[(maxnode*2+1)*GSIZE*2+i] != 0) {
                         OrigGraph[(maxnode*2+1)*GSIZE*2+i] = 2;
                         OrigGraph[i*GSIZE*2+(maxnode*2+1)] = 2;
                      }
                   }
                }
                // Now sweep through
                for (int i = 0; i < GSIZE*2*GSIZE*2; i++)
                   if (OrigGraph[i] == 2) OrigGraph[i] = 0;

                //_generate_result_file( bool(same_adj_Matrix==0 && same_path_Matrix==0),cpu_time,gpu_time,RANDOM_GSIZE);
           }
		currentrank += maxnodes.size();
}

}




void ATriaPlugin::output(std::string file) {
for (int i = GSIZE-1; i >= 0; i--)
           for (int j = 0; j < i; j++) {
              if (fabs(U[j]) > fabs(U[j+1])) {
                 float tmp = U[j];
                 U[j] = U[j+1];
                 U[j+1] = tmp;
                 string tmp2 = bacteria[j];
                 bacteria[j] = bacteria[j+1];
                 bacteria[j+1] = tmp2;
              }
           }

        std::ofstream noafile(file.c_str(), std::ios::out);
        //noafile << "Name\tCentrality\tRank" << endl;
        noafile << "Name\tRank" << endl;
        float min = 0;
        float max = 0;
        for (int i = 0; i < GSIZE; i++) {
           /*U[i] = fabs(U[i]);
           if (fabs(U[i]) > max)
              max = fabs(U[i]);
           if (fabs(U[i]) < min)
              min = fabs(U[i]);*/
           //noafile << bacteria[i] << "\t" << U[i] << "\t\t" << GSIZE-i << endl;
	   if (U[i] != 0)
           noafile << bacteria[i] << "\t" << U[i] << endl;// << "\t\t" << GSIZE-i << endl;
        }


}


void ATriaPlugin::_CPU_Floyd(float *G,/*int *Gpath,*/int N){//standard N^3 algo
        for(int k=0;k<N;++k)
        {
          for(int i=0;i<N;++i){
          // Changed to be += 2 because we only care about the shortest path STARTING from + vertices
          //for(int i=0;i<N;i+=2){
           for(int j=0;j<N;++j){
             // Important: i CAN equal k.  The shortest path from i to j can be through itself.
             if (i != j && j != k) {
                int curloc=i*N+j,loca=i*N+k,locb=k*N+j;
                //int pathloc=(i/2)*N+j;
                //if(G[curloc]>(G[loca]+G[locb])){
                int evenodd = i+j;
                if(  (evenodd % 2 == 0 && G[curloc]<(G[loca]*G[locb])) ||
                     (evenodd % 2 == 1 && G[curloc]>(G[loca]*G[locb])) ) {
                        G[curloc]=(G[loca]*G[locb]);
                        //Gpath[pathloc] = k;
                        //Gpath[curloc]=k;
                }
             }
           }
         }
        }
}
PluginProxy<ATriaPlugin> ATriaPluginProxy = PluginProxy<ATriaPlugin>("ATria", PluginManager::getInstance());
