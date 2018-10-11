#include "reactions.h"
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

void write_crds_complex(int c1, Complex *ind_com, Fullmol *bases) {

  cout <<"complex: "<<c1<<" size: "<<ind_com[c1].mysize<<" crds: "<<ind_com[c1].xcom<<' '<<ind_com[c1].ycom<<' '<<ind_com[c1].zcom<<endl;
  for(int i=0;i<ind_com[c1].mysize;i++){
    int mp=ind_com[c1].plist[i];
    cout <<"protein: "<<i<<" is: "<<mp<<endl;
    write_crds(bases, mp);
  }
  cout <<"Complex rad: "<<ind_com[c1].radR<<" Drot: "<<ind_com[c1].Drx<<' '<<ind_com[c1].Dry<<' '<<ind_com[c1].Drz<<endl;
  cout <<" Dtrans: "<<ind_com[c1].Dx<<' '<<ind_com[c1].Dy<<' '<<ind_com[c1].Dz<<endl;
  
}
