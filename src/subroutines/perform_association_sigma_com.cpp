#include "reactions.h"



void perform_association_sigma_com(int p1, Fullmol *bases, Complex *ind_com, Parms &plist, int **crosspart, int **cross_rxn, int **crossint, int **Rlist, int *i_home, double **probvec, int *Ncoup, int **mycoupled, int *p_home, int ci1, int ci2, int it, int *movestat, int *ncross, int bindrad, int *ncrosscom, double **traj, double **trajR, int *nBoundPairs)
{
  int p2 = crosspart[p1][ci1];//?????? WHY P2 AS FUNCTION OF P1
  int rxn1 = cross_rxn[p1][ci1];
  
  /*Associate proteins, move them to contact and update their free and bound lists//BASILIO: What are the free and bound list?
    For binding to lipids (they have D=0), move to contact but align vertically.
   */
  
  /*test for being in same complex to avoid moving binding interfaces too far*/
  //flag2=same_complex_test(p1, p2, bases, crossint, i_home, bindrad, rxn1, ci1, ci2);
  cout << "Associate between proteins: p1 " << p1 << ' ' << p2 << " interfaces: " << crossint[p1][ci1] << ' ' << crossint[p2][ci2] << " reaction: " << rxn1 << " pact; " << probvec[p1][ci1] <<  " iter: " << it << endl;
  
  //BASILIO: crossint[p1][ci1]  interface ci1 of protein p1...
  /*for self binding, use freeleg for clathrin type proteins.*/
  //BASILIO REVIEW IT (OSAMN SAID DO IT!!!)
  int cancel;
  int cind=bases[p1].mycomplex;
  int cind2=bases[p2].mycomplex;
  //cout <<"complex 1: "<<cind<<" size: "<<ind_com[cind].mysize<<" crds: "<<ind_com[cind].xcom<<' '<<ind_com[cind].ycom<<' '<<ind_com[cind].zcom<<endl;
  // for(int i=0;i<ind_com[cind].mysize;i++){
  //   int mp=ind_com[cind].plist[i];
  //   cout <<"protein: "<<i<<" is: "<<mp<<endl;
  //   write_crds(bases, mp);
  // }
  //  cout <<"complex 2: "<<cind2<<" size: "<<ind_com[cind2].mysize<<" crds: "<<ind_com[cind2].xcom<<' '<<ind_com[cind2].ycom<<' '<<ind_com[cind2].zcom<<endl;
  // for(int i=0;i<ind_com[cind2].mysize;i++){
  //   int mp=ind_com[cind2].plist[i];
  //   cout <<"protein: "<<i<<" is: "<<mp<<endl;
  //   write_crds(bases, mp);
  // }
  ////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////// MODIFY BY BASILIO ////////////////////////////////////
  if (bases[p1].protype == plist.pclath && bases[p2].protype == plist.pclath){
    //updating the member clatcount in complex ind_com[bases[p1].mycomplex].clatcount
    if(ind_com[bases[p1].mycomplex].Dz!=0 && ind_com[bases[p2].mycomplex].Dz!=0) {
      cancel=associate_freeleg(p1, p2, rxn1, crossint[p1][ci1], crossint[p2][ci2], bases, Rlist, i_home, ind_com, plist, bindrad, ncrosscom);// 3d versus 3d clat clat
      //      cout <<"associate freeleg "<<endl;
    } else if (ind_com[bases[p1].mycomplex].Dz==0 && ind_com[bases[p2].mycomplex].Dz==0) {
      cancel=associate_freeleg(p1, p2, rxn1, crossint[p1][ci1], crossint[p2][ci2], bases, Rlist, i_home, ind_com, plist, bindrad, ncrosscom);// 2d versus 2d clat clat
      //cout <<"associate freeleg "<<endl;
    } else {
      cancel=associate_2dclat_3dclat_lipid_sigma(p1, p2, rxn1, crossint[p1][ci1], crossint[p2][ci2], bases, Rlist, i_home, ind_com, plist, bindrad, ncrosscom); // 3d vs 2d and 2d vs 3d clat clat
      //cout <<"associate 2dclat_3dclat "<<endl;
    }
  }else if (bases[p1].protype == plist.pclath || bases[p2].protype == plist.pclath) { //special ap2 clath
    
    cancel=associate_ap2_clath_sigma(p1, p2, rxn1, crossint[p1][ci1], crossint[p2][ci2], bases, Rlist, i_home, ind_com, plist, bindrad, ncrosscom);
    //cout <<"associate ap2_clath"<<endl;
    //After this line, all interaction is without clathrin
    
  }else if(bases[p1].Dz==0 || bases[p2].Dz==0){
    //at least one partner is in 2D always.//??? This is just for PI2 ? Clath can be in the membrane too...
    if(bases[p1].Dz==0 && bases[p2].Dz==0){
      /*If both proteins are 2D, just place them to contact.*/
      cancel=associate_zsigma(p1, p2, rxn1, crossint[p1][ci1], crossint[p2][ci2], bases, Rlist, i_home, ind_com, plist, bindrad, ncrosscom);
    }else
      cancel=associate_lipid_sigmaAP2(p1, p2, rxn1, crossint[p1][ci1], crossint[p2][ci2], bases, Rlist, i_home, ind_com, plist, bindrad, ncrosscom);//this comes from associate_freelegPBCCELL
    //cout <<"associate lipid_sigmaAP2 "<<endl;
    // 		  associate_zsigmaPBCCELL(p1, p2, rxn1, crossint[p1][ci1], crossint[p2][ci2], bases, Rlist, i_home, ind_com, plist, bindrad, ncrosscom);//2D association
    // 	  }else{
    // 		  associate_lipid_sigmaPBCCELLAP2(p1, p2, rxn1, crossint[p1][ci1], crossint[p2][ci2], bases, Rlist, i_home, ind_com, plist, bindrad, ncrosscom);//this comes from associate_freelegPBCCELL
    // 	  }
    // do you need the next section?
  }else{
    cancel=associate_zsigma(p1, p2, rxn1, crossint[p1][ci1], crossint[p2][ci2], bases, Rlist, i_home, ind_com, plist, bindrad, ncrosscom); // when I modified by associate_ap2_clath_sigmaPBCCEL ...
    //cout <<"associate zsigm "<<endl;
  }
    // ap2 interact not bertically with pip2
  ///////////////////// MODIFICATION END HERE /////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////
  if(cancel==0){//this means you continued with association.
    cind=bases[p1].mycomplex;
    double *M=new double[9];
    ap2_coupledrxn_add(rxn1, Ncoup, mycoupled, bases, Rlist, i_home, p_home, p1, p2);
    reflect_complex_rad_rot(p1, bases, ind_com, plist);
    update_Nboundpairs(bases[p1].protype, bases[p2].protype, 1, plist, nBoundPairs);
    
    /*Update the copied proteins trajectory*/
    traj[cind2][0]=traj[plist.ntotalcomplex][0];
    traj[cind2][1]=traj[plist.ntotalcomplex][1];
    traj[cind2][2]=traj[plist.ntotalcomplex][2];

    trajR[cind2][0]=trajR[plist.ntotalcomplex][0];
    trajR[cind2][1]=trajR[plist.ntotalcomplex][1];
    trajR[cind2][2]=trajR[plist.ntotalcomplex][2];
    /*Remove p1 and p2 from the list of potential reaction partners so they don't try to
      associate again in this turn. They will also not try to avoid overlap anymore, and proteins
      close by will not try to avoid overlapping them this turn. 
    */
    //  remove_reaction_all_ncom(p1, p2, ncross, crosspart, probvec,  cross_rxn, crossint, ncrosscom, bases);
    /*
      To allow proteins that just associated to continue avoiding overlap, set probabilities to zero, and trajectories to zero.
    */
    traj[bases[p1].mycomplex][0]=0;
    traj[bases[p1].mycomplex][1]=0;
    traj[bases[p1].mycomplex][2]=0;

    trajR[bases[p1].mycomplex][0]=0;
    trajR[bases[p1].mycomplex][1]=0;
    trajR[bases[p1].mycomplex][2]=0;

    
    /*Set probability to zero */
    remove_one_prob_all(p1, ncross, crosspart, probvec, cross_rxn, crossint);
    /*Set probability to zero */
    remove_one_prob_all(p2, ncross, crosspart, probvec, cross_rxn, crossint);
    ncross[p1]=-1;
    ncross[p2]=-1;
    /*These should both be the same, now that p1 and p2 are in the same complex.*/
    ncrosscom[bases[p1].mycomplex]=-1;//you won't avoid them, but they should avoid you. you won't move.
    ncrosscom[bases[p2].mycomplex]=-1;
    
    /*Since these proteins have moved to associate and taken their complex with them,
      don't allow any proteins in their complex to move again.
    */
    set_movestat_zero(p1, bases, ind_com, movestat);
    //    cout <<"FINAL POSITIONS "<<endl;
    /*If cind was Ntotalcomplex, it will be copied now to c2*/
    //cout <<"complex 1: "<<cind<<" size: "<<ind_com[cind].mysize<<" crds: "<<ind_com[cind].xcom<<' '<<ind_com[cind].ycom<<' '<<ind_com[cind].zcom<<endl;
    // for(int i=0;i<ind_com[cind].mysize;i++){
    //   int mp=ind_com[cind].plist[i];
    //   cout <<"protein: "<<i<<" is: "<<mp<<endl;
    //   write_crds(bases, mp);
    // }
    delete[] M;
  }else
    cout <<"Canceled reaction due to overlap crash of protein structure "<<endl;
}
