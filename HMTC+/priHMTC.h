#ifndef PRIHMTC_H
#define PRIHMTC_H
#include"lib/pairing_3.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <set>
#include <tr1/unordered_map>
#include <tr1/unordered_set>
#include <cstring>
using namespace std;

struct PP
{
    G1 g;
    G1 h;
    G2 w_g;
    G2 w_h;
};

struct MSK
{
    Big xta;
};

struct SKCS
{
    Big x1;
    Big x2;
    Big y;
    Big z1;
    Big z2;
    vector<Big> w; // len=nk+1
};

struct PKCS
{
    G2 w_X1;
    G2 w_X2;
    G2 w_Y;
    G2 w_Z1;
    G2 w_Z2;
    vector<G1> W; // len=nk+1
    vector<G2> w_W; // len=nk+1
};

struct SKMU
{
    Big xu;
    G1 S0;
    G1 S1;
    G1 S2;
    vector<G1> A; // len=nk+1
    vector<G1> Btp; // len=ntp+1,Btp[0]=null
};

struct PKMU
{
    G1 Yu;
};

struct Pi1
{
    G1 h_Yu;
    Big o_xu;
};

struct Pi2
{
    G1 h_Yu;
    G1 h_Yu_p;
    Big o_xu;
    Big o_xu_p;
    vector<G1> h_Rtp_list;
};


struct Pi3{
    G2 h_T1;
    G2 h_T2;
    G1 h_T3;
    Big o_xu;
    Big o_t;
};

struct Sigma
{
    G1 sigma1;
    G1 sigma2;
    G1 sigma3;
    G2 w_T1;
    G2 w_T2;
    G1 T3;
    Pi3 pi3;
};


struct h_Sigma
{
    G1 sigma1;
    G1 sigma2;
    G1 sigma3;
    G2 w_T1;
    G2 w_T2;
    G1 T3;
    G2 w_D1;
    G2 w_D2;
    Pi3 pi3;
};



struct Task_tp{
    vector<int> task_id;
}; // len<nk, 0<val<=nk

struct AS_tp{
    vector<bool> ass_id;
}; // len=nk+1, ass_id[0]=null

struct Task_all{
    vector<Task_tp> task_all;
}; // len=ntp+1,task_all[0]=null

struct AS{
    vector<AS_tp> as;
}; // len=ntp+1,task_all[0]=null


struct Dataset{
    vector<char*> data;
};  // len=len(Task_tp)

struct AggData{
    vector<Dataset> datasets;
    Task_tp disclosure_policy;
}; // len=ns

struct AggSigma{
    vector<h_Sigma> signatures;
}; // len=ns


class Pri_HMTC
{

public:
    PFC *pfc;
    int nk;
    int ntp;
    int nu=0; // number of users;
    PP pp;
    PKCS pkcs;
    vector<PKMU> mulist; //mulist[i]--> i is ID_MU
    vector <G1> Lds;

    Pri_HMTC(PFC *p);
    ~Pri_HMTC();
    int Setup(int ntp,PP& pp,MSK& msk);
    int CSKeygen(int nk,SKCS& skcs,PKCS& pkcs);
    int MUKeygen(Task_all assignment,SKCS skcs,SKMU& skmu,PKMU& pkmu);
    int Update(Task_all assignment_pre,Task_all assignment_new,SKCS skcs,SKMU& skmu_pre,PKMU& pkmu_pre,SKMU& skmu_new,PKMU& pkmu_new);
    int Sign(SKMU skmu,int tp,Task_tp tasktp,Dataset dataset,char* msg,Sigma& signature);
    int Sign_t(SKMU skmu,int tp,Task_tp tasktp,Dataset dataset,char* msg,Sigma& signature, Big &ot);
    int Verify(int tp,Task_tp tasktp,Dataset dataset,char* msg,Sigma& signature);
    int Redact(SKCS skcs,Task_tp tasktp,Dataset dataset,char* msg,Task_tp disclosure,Sigma& signature,h_Sigma& redsig);
    int RedactVerify(Dataset dataset,Task_tp disclosure,h_Sigma redsig);
    // int AggVerify(AggData aggdata,AggSigma aggsigma);
    int Reveal(MSK msk,h_Sigma redsig,int& idmu,PKMU& pkmu);

    //print task assignments
    AS_tp Map(Task_tp& tasktp);
    void show_Tasktp(Task_tp tasktp);
    void show_AStp(AS_tp astp);
    void show_Taskall(Task_all taskall);
    void show_AS(AS as);

    //measure the time cost of MUKeygen
    int MUKeygen_t(Task_all assignment,SKCS skcs,SKMU& skmu,PKMU& pkmu,double& time_mu,double& time_cs);
};

#endif // PRIHMTC_H
