#include "priHMTC.h"
#include "lib/bn_transfer.h"
#include "lib/bn_struct.h"
#include <iomanip>
#include <algorithm>
#define TEST_TIME 1
Pri_HMTC::Pri_HMTC(PFC *p)
{
    pfc=p;
}
Pri_HMTC::~Pri_HMTC()
{

}
int Pri_HMTC::Setup(int n_tp,PP& rpp,MSK& msk){
    ntp=n_tp;

    Big xta;

    pfc->random(rpp.g);
    pfc->random(rpp.w_g);
    pfc->random(xta);

    rpp.h=pfc->mult(rpp.g,xta);
    rpp.w_h=pfc->mult(rpp.w_g,xta);

    msk.xta=xta;
    pp=rpp;
    return 0;
}


int Pri_HMTC::CSKeygen(int n_k,SKCS& skcs,PKCS& rpkcs){
    nk=n_k;
    pfc->random(skcs.x1);
    pfc->random(skcs.x2);
    pfc->random(skcs.y);
    pfc->random(skcs.z1);
    pfc->random(skcs.z2);

    rpkcs.w_X1=pfc->mult(pp.w_g,skcs.x1);
    rpkcs.w_X2=pfc->mult(pp.w_g,skcs.x2);
    rpkcs.w_Y=pfc->mult(pp.w_g,skcs.y);
    rpkcs.w_Z1=pfc->mult(pp.w_g,skcs.z1);
    rpkcs.w_Z2=pfc->mult(pp.w_g,skcs.z2);

    for(int i=0;i<=nk;i++){
        Big w;
        pfc->random(w);
        skcs.w.push_back(w);

        G1 W=pfc->mult(pp.g,w);
        G2 w_W=pfc->mult(pp.w_g,w);
        rpkcs.W.push_back(W);
        rpkcs.w_W.push_back(w_W);

    }

    pkcs=rpkcs;

    return 0;
}
int Pri_HMTC::MUKeygen(Task_all assignment,SKCS skcs,SKMU& skmu,PKMU& pkmu){

    //Step-1 MU
    Big xu;
    pfc->random(xu);
    G1 Yu;
    Yu=pfc->mult(pp.g,xu);

    //generate Pi1
    Big h_xu;
    pfc->random(h_xu);
    G1 h_Yu;
    h_Yu=pfc->mult(pp.g,h_xu);
    pfc->start_hash();
    pfc->add_to_hash(Yu);
    pfc->add_to_hash(h_Yu);
    Big c=pfc->finish_hash_to_group();
    Big o_xu=h_xu-c*xu;

    Pi1 pi1;
    pi1.h_Yu=h_Yu;
    pi1.o_xu=o_xu;

    pkmu.Yu=Yu;

    //Step-2 CS
    //verify pi1
    pfc->start_hash();
    pfc->add_to_hash(Yu);
    pfc->add_to_hash(h_Yu);
    Big c2=pfc->finish_hash_to_group();
    if(pi1.h_Yu!=pfc->mult(Yu,c2)+pfc->mult(pp.g,o_xu)){
        printf("verification of Pi_1 failed\n");
        return -1;
    }

    //turn assignment to AS
    AS as;
    AS_tp astp0;
    as.as.push_back(astp0);
    for(int i=1;i<=ntp;i++){
        as.as.push_back(Map(assignment.task_all[i]));
    }

    // cout<<"~~~~~~AS:"<<endl;
    // show_AS(as);

    //key generation
    Big s;
    pfc->random(s);
    skmu.S0=pfc->mult(pp.h,s);
    skmu.S1=pfc->mult(pp.g,s);
    skmu.S2=pfc->mult(pfc->mult(pp.g,skcs.x1)+pfc->mult(Yu,skcs.y),s);

    for(int i=0;i<=nk;i++){
        G1 Ai=pfc->mult(pkcs.W[i],s);
        skmu.A.push_back(Ai);
    }

    G1 btmp;
    skmu.Btp.push_back(btmp);
    for(int tp=1;tp<=ntp;tp++){

        //map binary AS[tp] to char*, implicity require that nk<100
        char astp[101]={0};
        for(int i=1;i<=nk;i++) 
            astp[i-1]=as.as[tp].ass_id[i]; 
        //get H(AS[tp])
        pfc->start_hash();
        pfc->add_to_hash(astp);
        Big hastp=pfc->finish_hash_to_group();

        G1 Btp=pfc->mult(pfc->mult(pp.g,skcs.x2+skcs.z1*tp+skcs.z2*hastp)+pfc->mult(Yu,skcs.y),s);
        skmu.Btp.push_back(Btp);
    }

    mulist.push_back(pkmu);

    //Step-3 MU
    GT p1,p2,p3,p4;
    skmu.xu=xu;
    G2 sigW=pkcs.w_W[0];
    for(int i=1;i<=nk;i++)
        sigW=sigW+pkcs.w_W[i];
    p1=pfc->pairing(pkcs.w_X1+pfc->mult(pkcs.w_Y,skmu.xu)+pp.w_h+sigW,skmu.S1);


    G1 sigA=skmu.A[0];
    for(int i=1;i<=nk;i++)
        sigA=sigA+skmu.A[i];
    p2=pfc->pairing(pp.w_g,skmu.S0+skmu.S2+sigA);

    if(p1!=p2){
        printf("verification of tb-SK_MU failed 1\n");
        return -1;
    }

    G2 sigZ=pfc->mult(pkcs.w_Z1,Big(1));
    for(int tp=1;tp<=ntp;tp++){

        //map binary AS[tp] to char*, implicity require that nk<=100
        char astp[101]={0};
        for(int i=1;i<=nk;i++) 
            astp[i-1]=as.as[tp].ass_id[i]; 
        //get H(AS[tp])
        pfc->start_hash();
        pfc->add_to_hash(astp);
        Big hastp=pfc->finish_hash_to_group();

        if(tp==1)
            sigZ=sigZ+pfc->mult(pkcs.w_Z2,hastp);
        else
            sigZ=sigZ+pfc->mult(pkcs.w_Z1,Big(tp))+pfc->mult(pkcs.w_Z2,hastp);
    }
    p3=pfc->pairing(pfc->mult(pkcs.w_X2+pfc->mult(pkcs.w_Y,skmu.xu),Big(ntp))+sigZ,skmu.S1);

    G1 sigB=skmu.Btp[1];
    for(int tp=2;tp<=ntp;tp++){
        sigB=sigB+skmu.Btp[tp];
    }
    p4=pfc->pairing(pp.w_g,sigB);

    if(p3!=p4){
        printf("verification of tb-SK_MU failed 2\n");
        return -1;
    }
    

    // printf("verification of tb-SK_MU successful\n");
    return 0;
}





int Pri_HMTC::MUKeygen_t(Task_all assignment,SKCS skcs,SKMU& skmu,PKMU& pkmu,double& time_mu,double& time_cs){
    double sum1=0,sum2=0;
    clock_t start,finish;
    start=clock();


    //Step-1 MU
    Big xu;
    pfc->random(xu);
    G1 Yu;
    Yu=pfc->mult(pp.g,xu);

    //generate Pi1
    Big h_xu;
    pfc->random(h_xu);
    G1 h_Yu;
    h_Yu=pfc->mult(pp.g,h_xu);
    pfc->start_hash();
    pfc->add_to_hash(Yu);
    pfc->add_to_hash(h_Yu);
    Big c=pfc->finish_hash_to_group();
    Big o_xu=h_xu-c*xu;

    Pi1 pi1;
    pi1.h_Yu=h_Yu;
    pi1.o_xu=o_xu;

    pkmu.Yu=Yu;


    finish=clock();
    sum1 = (double)(finish-start)/(CLOCKS_PER_SEC*TEST_TIME);



    start=clock();

    //Step-2 CS
    //verify pi1
    pfc->start_hash();
    pfc->add_to_hash(Yu);
    pfc->add_to_hash(h_Yu);
    Big c2=pfc->finish_hash_to_group();
    if(pi1.h_Yu!=pfc->mult(Yu,c2)+pfc->mult(pp.g,o_xu)){
        printf("verification of Pi_1 failed\n");
        return -1;
    }

    //turn assignment to AS
    AS as;
    AS_tp astp0;
    as.as.push_back(astp0);
    for(int i=1;i<=ntp;i++){
        as.as.push_back(Map(assignment.task_all[i]));
    }

    // cout<<"~~~~~~AS:"<<endl;
    // show_AS(as);

    //key generation
    Big s;
    pfc->random(s);
    skmu.S0=pfc->mult(pp.h,s);
    skmu.S1=pfc->mult(pp.g,s);
    skmu.S2=pfc->mult(pfc->mult(pp.g,skcs.x1)+pfc->mult(Yu,skcs.y),s);

    for(int i=0;i<=nk;i++){
        G1 Ai=pfc->mult(pkcs.W[i],s);
        skmu.A.push_back(Ai);
    }

    G1 btmp;
    skmu.Btp.push_back(btmp);
    for(int tp=1;tp<=ntp;tp++){

        //map binary AS[tp] to char*, implicity require that nk<100
        char astp[101]={0};
        for(int i=1;i<=nk;i++) 
            astp[i-1]=as.as[tp].ass_id[i]; 
        //get H(AS[tp])
        pfc->start_hash();
        pfc->add_to_hash(astp);
        Big hastp=pfc->finish_hash_to_group();

        G1 Btp=pfc->mult(pfc->mult(pp.g,skcs.x2+skcs.z1*tp+skcs.z2*hastp)+pfc->mult(Yu,skcs.y),s);
        skmu.Btp.push_back(Btp);
    }

    mulist.push_back(pkmu);

    finish=clock();
    sum2 = (double)(finish-start)/(CLOCKS_PER_SEC*TEST_TIME);



    start=clock();

    //Step-3 MU
    GT p1,p2,p3,p4;
    skmu.xu=xu;
    G2 sigW=pkcs.w_W[0];
    for(int i=1;i<=nk;i++)
        sigW=sigW+pkcs.w_W[i];
    p1=pfc->pairing(pkcs.w_X1+pfc->mult(pkcs.w_Y,skmu.xu)+pp.w_h+sigW,skmu.S1);


    G1 sigA=skmu.A[0];
    for(int i=1;i<=nk;i++)
        sigA=sigA+skmu.A[i];
    p2=pfc->pairing(pp.w_g,skmu.S0+skmu.S2+sigA);

    if(p1!=p2){
        printf("verification of tb-SK_MU failed 1\n");
        return -1;
    }

    G2 sigZ=pfc->mult(pkcs.w_Z1,Big(1));
    for(int tp=1;tp<=ntp;tp++){

        //map binary AS[tp] to char*, implicity require that nk<=100
        char astp[101]={0};
        for(int i=1;i<=nk;i++) 
            astp[i-1]=as.as[tp].ass_id[i]; 
        //get H(AS[tp])
        pfc->start_hash();
        pfc->add_to_hash(astp);
        Big hastp=pfc->finish_hash_to_group();

        if(tp==1)
            sigZ=sigZ+pfc->mult(pkcs.w_Z2,hastp);
        else
            sigZ=sigZ+pfc->mult(pkcs.w_Z1,Big(tp))+pfc->mult(pkcs.w_Z2,hastp);
    }
    p3=pfc->pairing(pfc->mult(pkcs.w_X2+pfc->mult(pkcs.w_Y,skmu.xu),Big(ntp))+sigZ,skmu.S1);

    G1 sigB=skmu.Btp[1];
    for(int tp=2;tp<=ntp;tp++){
        sigB=sigB+skmu.Btp[tp];
    }
    p4=pfc->pairing(pp.w_g,sigB);

    if(p3!=p4){
        printf("verification of tb-SK_MU failed 2\n");
        return -1;
    }
    finish=clock();
    sum1+= (double)(finish-start)/(CLOCKS_PER_SEC*TEST_TIME);

    
    time_mu=sum1;
    time_cs=sum2;

    printf("verification of MUKeygen successful\n");
    return 0;
}

int Pri_HMTC::Update(Task_all assignment_pre,Task_all assignment_new,SKCS skcs,SKMU& skmu_pre,PKMU& pkmu_pre,SKMU& skmu_new,PKMU& pkmu_new){

    //Step-1 MU

    vector<G1> Rtp_list;
    for(int tp=1;tp<=ntp;tp++){
        Big xtp=skmu_pre.xu+Big(tp);
        Big inv_xtp=pfc->Zpinverse(xtp);
        G1 Rtp=pfc->mult(pp.g,inv_xtp);
        Rtp_list.push_back(Rtp);
    }

    Big xu_new;
    pfc->random(xu_new);
    G1 Yu;
    Yu=pfc->mult(pp.g,xu_new);

    //generate Pi2
    Big h_xu,h_xu_p;
    pfc->random(h_xu);
    pfc->random(h_xu_p);
    G1 h_Yu,h_Yu_p;
    h_Yu=pfc->mult(pp.g,h_xu);
    h_Yu_p=pfc->mult(pp.g,h_xu_p);
    pfc->start_hash();
    pfc->add_to_hash(h_xu);
    pfc->add_to_hash(h_Yu_p);
    Big c=pfc->finish_hash_to_group();
    Big o_xu=h_xu-c*skmu_pre.xu;
    Big o_xu_p=h_xu_p-c*xu_new;

    Pi2 pi2;
    pi2.h_Yu=h_Yu;
    pi2.h_Yu_p=h_Yu_p;
    pi2.o_xu=o_xu;
    pi2.o_xu_p=o_xu_p;

    pkmu_new.Yu=Yu;

    //Step-2 CS
    //verify pi2
    pfc->start_hash();
    pfc->add_to_hash(h_xu);
    pfc->add_to_hash(h_Yu_p);
    Big c2=pfc->finish_hash_to_group();
    if((pi2.h_Yu!=pfc->mult(pkmu_pre.Yu,c2)+pfc->mult(pp.g,o_xu))||(pi2.h_Yu_p!=pfc->mult(Yu,c2)+pfc->mult(pp.g,o_xu_p))){
        printf("verification of Pi_2 failed\n");
        return -1;
    }

    //isnert Rtp into Lds
    for(int i=0;i<ntp;i++){
        Lds.push_back(Rtp_list[i]);
    }

    //construct the updated AS
    AS as;
    AS_tp astp0;
    as.as.push_back(astp0);
    for(int i=1;i<=ntp;i++){
        as.as.push_back(Map(assignment_new.task_all[i]));
    }

    // cout<<"~~~~~~AS:"<<endl;
    // show_AS(as);

    //key generation
    Big s;
    pfc->random(s);
    skmu_new.S0=pfc->mult(pp.h,s);
    skmu_new.S1=pfc->mult(pp.g,s);
    skmu_new.S2=pfc->mult(pfc->mult(pp.g,skcs.x1)+pfc->mult(Yu,skcs.y),s);

    for(int i=0;i<=nk;i++){
        G1 Ai=pfc->mult(pkcs.W[i],s);
        skmu_new.A.push_back(Ai);
    }

    G1 btmp;
    skmu_new.Btp.push_back(btmp);
    for(int tp=1;tp<=ntp;tp++){

        //map binary AS[tp] to char*, implicity require that nk<100
        char astp[101]={0};
        for(int i=1;i<=nk;i++) 
            astp[i-1]=as.as[tp].ass_id[i]; 
        //get H(AS[tp])
        pfc->start_hash();
        pfc->add_to_hash(astp);
        Big hastp=pfc->finish_hash_to_group();

        G1 Btp=pfc->mult(pfc->mult(pp.g,skcs.x2+skcs.z1*tp+skcs.z2*hastp)+pfc->mult(Yu,skcs.y),s);
        skmu_new.Btp.push_back(Btp);
    }

    mulist.push_back(pkmu_new);

    //Step-3 MU
    GT p1,p2,p3,p4;
    skmu_new.xu=xu_new;
    G2 sigW=pkcs.w_W[0];
    for(int i=1;i<=nk;i++)
        sigW=sigW+pkcs.w_W[i];
    p1=pfc->pairing(pkcs.w_X1+pfc->mult(pkcs.w_Y,skmu_new.xu)+pp.w_h+sigW,skmu_new.S1);


    G1 sigA=skmu_new.A[0];
    for(int i=1;i<=nk;i++)
        sigA=sigA+skmu_new.A[i];
    p2=pfc->pairing(pp.w_g,skmu_new.S0+skmu_new.S2+sigA);

    if(p1!=p2){
        printf("verification of Update failed 1\n");
        return -1;
    }

    G2 sigZ=pfc->mult(pkcs.w_Z1,Big(1));
    for(int tp=1;tp<=ntp;tp++){

        //map binary AS[tp] to char*, implicity require that nk<=100
        char astp[101]={0};
        for(int i=1;i<=nk;i++) 
            astp[i-1]=as.as[tp].ass_id[i]; 
        //get H(AS[tp])
        pfc->start_hash();
        pfc->add_to_hash(astp);
        Big hastp=pfc->finish_hash_to_group();

        if(tp==1)
            sigZ=sigZ+pfc->mult(pkcs.w_Z2,hastp);
        else
            sigZ=sigZ+pfc->mult(pkcs.w_Z1,Big(tp))+pfc->mult(pkcs.w_Z2,hastp);
    }
    p3=pfc->pairing(pfc->mult(pkcs.w_X2+pfc->mult(pkcs.w_Y,skmu_new.xu),Big(ntp))+sigZ,skmu_new.S1);

    G1 sigB=skmu_new.Btp[1];
    for(int tp=2;tp<=ntp;tp++){
        sigB=sigB+skmu_new.Btp[tp];
    }
    p4=pfc->pairing(pp.w_g,sigB);

    if(p3!=p4){
        printf("verification of Update failed 2\n");
        return -1;
    }
    

    printf("verification of Update successful\n");
    return 0;
}



int Pri_HMTC::Sign_t(SKMU skmu,int tp,Task_tp tasktp,Dataset dataset,char* msg,Sigma& signature, Big &ot){
    Big r,t;
    pfc->random(r);
    pfc->random(t);
    signature.w_T1=pfc->mult(pp.w_g,t);
    signature.w_T2=pfc->mult(pp.w_h,t)+pfc->mult(pkcs.w_Y,skmu.xu);

    Big xtp=skmu.xu+Big(tp);
    Big inv_xtp=pfc->Zpinverse(xtp);
    signature.T3=pfc->mult(pp.g,inv_xtp);

    signature.sigma1=pfc->mult(skmu.S1,r);

    pfc->start_hash();
    pfc->add_to_hash(msg);
    Big hmsg=pfc->finish_hash_to_group();
    G1 tmp=skmu.S2+pfc->mult(skmu.S0,t)+pfc->mult(skmu.A[0],hmsg);
    for(int i=0;i<tasktp.task_id.size();i++){
        pfc->start_hash();
        pfc->add_to_hash(dataset.data[i]);
        Big hdata=pfc->finish_hash_to_group();

        tmp=tmp+pfc->mult(skmu.A[tasktp.task_id[i]],hdata);
        //current task: tasktp.task_id[i]
        //current data: dataset.data[i]
    }
    signature.sigma2=pfc->mult(tmp,r);
    signature.sigma3=pfc->mult(skmu.Btp[tp]+pfc->mult(skmu.S0,t),r);


    //generate Pi3
    Big h_xu,h_t;
    pfc->random(h_xu);
    pfc->random(h_t);
    G2 h_T1,h_T2;
    h_T1=pfc->mult(pp.w_g,h_t);
    h_T2=pfc->mult(pp.w_h,h_t)+pfc->mult(pkcs.w_Y,h_xu);
    G1 h_T3=pfc->mult(signature.T3,h_xu);

    pfc->start_hash();
    pfc->add_to_hash(signature.sigma1);
    pfc->add_to_hash(signature.sigma2);
    pfc->add_to_hash(signature.sigma3);
    pfc->add_to_hash(signature.w_T1);
    pfc->add_to_hash(signature.w_T2);
    pfc->add_to_hash(signature.T3);
    pfc->add_to_hash(h_T1);
    pfc->add_to_hash(h_T2);
    pfc->add_to_hash(h_T3);
    Big c=pfc->finish_hash_to_group();
    Big o_xu=h_xu-c*skmu.xu;
    Big o_t=h_t-c*t;

    Pi3 pi3;
    pi3.h_T1=h_T1;
    pi3.h_T2=h_T2;
    pi3.h_T3=h_T3;
    pi3.o_xu=o_xu;
    pi3.o_t=o_t;

    signature.pi3=pi3;
    ot=t;

    return 0;
}

int Pri_HMTC::Sign(SKMU skmu,int tp,Task_tp tasktp,Dataset dataset,char* msg,Sigma& signature){
    Big r,t;
    pfc->random(r);
    pfc->random(t);
    signature.w_T1=pfc->mult(pp.w_g,t);
    signature.w_T2=pfc->mult(pp.w_h,t)+pfc->mult(pkcs.w_Y,skmu.xu);

    Big xtp=skmu.xu+Big(tp);
    Big inv_xtp=pfc->Zpinverse(xtp);
    signature.T3=pfc->mult(pp.g,inv_xtp);

    signature.sigma1=pfc->mult(skmu.S1,r);

    pfc->start_hash();
    pfc->add_to_hash(msg);
    Big hmsg=pfc->finish_hash_to_group();
    G1 tmp=skmu.S2+pfc->mult(skmu.S0,t)+pfc->mult(skmu.A[0],hmsg);
    for(int i=0;i<tasktp.task_id.size();i++){
        pfc->start_hash();
        pfc->add_to_hash(dataset.data[i]);
        Big hdata=pfc->finish_hash_to_group();

        tmp=tmp+pfc->mult(skmu.A[tasktp.task_id[i]],hdata);
        //current task: tasktp.task_id[i]
        //current data: dataset.data[i]
    }
    signature.sigma2=pfc->mult(tmp,r);
    signature.sigma3=pfc->mult(skmu.Btp[tp]+pfc->mult(skmu.S0,t),r);


    //generate Pi3
    Big h_xu,h_t;
    pfc->random(h_xu);
    pfc->random(h_t);
    G2 h_T1,h_T2;
    h_T1=pfc->mult(pp.w_g,h_t);
    h_T2=pfc->mult(pp.w_h,h_t)+pfc->mult(pkcs.w_Y,h_xu);
    G1 h_T3=pfc->mult(signature.T3,h_xu);

    pfc->start_hash();
    pfc->add_to_hash(signature.sigma1);
    pfc->add_to_hash(signature.sigma2);
    pfc->add_to_hash(signature.sigma3);
    pfc->add_to_hash(signature.w_T1);
    pfc->add_to_hash(signature.w_T2);
    pfc->add_to_hash(signature.T3);
    pfc->add_to_hash(h_T1);
    pfc->add_to_hash(h_T2);
    pfc->add_to_hash(h_T3);
    Big c=pfc->finish_hash_to_group();
    Big o_xu=h_xu-c*skmu.xu;
    Big o_t=h_t-c*t;

    Pi3 pi3;
    pi3.h_T1=h_T1;
    pi3.h_T2=h_T2;
    pi3.h_T3=h_T3;
    pi3.o_xu=o_xu;
    pi3.o_t=o_t;

    signature.pi3=pi3;


    return 0;
}


int Pri_HMTC::Verify(int tp,Task_tp tasktp,Dataset dataset,char* msg,Sigma& signature){

    //Verify T3
    std::vector<G1>::iterator it = std::find(Lds.begin(), Lds.end(), signature.T3);
    
    if (it != Lds.end()) {
        std::cout << "T3 exists in Lds." << std::endl;
        return -1;
    }


    //verify Pi3
    Pi3 pi3=signature.pi3;
    pfc->start_hash();
    pfc->add_to_hash(signature.sigma1);
    pfc->add_to_hash(signature.sigma2);
    pfc->add_to_hash(signature.sigma3);
    pfc->add_to_hash(signature.w_T1);
    pfc->add_to_hash(signature.w_T2);
    pfc->add_to_hash(signature.T3);
    pfc->add_to_hash(pi3.h_T1);
    pfc->add_to_hash(pi3.h_T2);
    pfc->add_to_hash(pi3.h_T3);
    Big c=pfc->finish_hash_to_group();

    G2 A,B,C,D;
    A=pi3.h_T1;
    B=pfc->mult(signature.w_T1,c)+pfc->mult(pp.w_g,pi3.o_t);

    if(A!=B){
        printf("verification of Pi_3 failed 1\n");
        return -1;
    }
    C=pi3.h_T2;
    D=pfc->mult(signature.w_T2,c)+pfc->mult(pp.w_h,pi3.o_t)+pfc->mult(pkcs.w_Y,pi3.o_xu);
    if(C!=D){
        printf("verification of Pi_3 failed 2\n");
        return -1;
    }

    if((pi3.h_T3)!=pfc->mult(signature.T3,pi3.o_xu)+pfc->mult(-pfc->mult(signature.T3,Big(tp))+pp.g,c)){
        printf("verification of Pi_3 failed 3\n");
        return -1;
    }

    AS_tp as_tp=Map(tasktp);
    char astp[101]={0};
    for(int i=1;i<=nk;i++) 
        astp[i-1]=as_tp.ass_id[i]; 
    //get H(AS[tp])
    pfc->start_hash();
    pfc->add_to_hash(astp);
    Big hastp=pfc->finish_hash_to_group();

    pfc->start_hash();
    pfc->add_to_hash(msg);
    Big hmsg=pfc->finish_hash_to_group();


    G2 tmp=pkcs.w_X1+pkcs.w_X2+signature.w_T2+signature.w_T2+pfc->mult(pkcs.w_Z1,Big(tp))+pfc->mult(pkcs.w_Z2,hastp)+pfc->mult(pkcs.w_W[0],hmsg);
    for(int i=0;i<tasktp.task_id.size();i++){
        pfc->start_hash();
        pfc->add_to_hash(dataset.data[i]);
        Big hdata=pfc->finish_hash_to_group();
        tmp=tmp+pfc->mult(pkcs.w_W[tasktp.task_id[i]],hdata);
    }
    
    if(pfc->pairing(tmp,signature.sigma1)!=pfc->pairing(pp.w_g,signature.sigma2+signature.sigma3)){
        printf("verification of signature failed\n");
        return -1;
    }

    
    
    return 0;
}
int Pri_HMTC::Redact(SKCS skcs,Task_tp tasktp,Dataset dataset,char* msg,Task_tp disclosure,Sigma& signature,h_Sigma& redsig){
    // cout<<"tasktp:";
    // show_Tasktp(tasktp);
    // cout<<"disc:";
    // show_Tasktp(disclosure);
    
    vector<int> D,ND;
    int j=0;
    for(int i=0;i<tasktp.task_id.size();i++){
        if(tasktp.task_id[i]==disclosure.task_id[j]){
            D.push_back(i);
            j++;
        }else{
            ND.push_back(i);
        }
    }


    // Task_tp sD,sND;
    // sD.task_id=D;
    // sND.task_id=ND;
    // cout<<"sD:";
    // show_Tasktp(sD);
    // cout<<"sND:";
    // show_Tasktp(sND);


    pfc->start_hash();
    pfc->add_to_hash(msg);
    Big hmsg=pfc->finish_hash_to_group();
    Big tmp1,tmp2,tmp3;
    tmp1=tmp2=tmp3=Big(0);
    if(!tmp1.iszero()){
        printf("zero test error");
        return -1;
    }

    for(int j=0;j<ND.size();j++){
        pfc->start_hash();
        pfc->add_to_hash(dataset.data[ND[j]]);
        Big hdata=pfc->finish_hash_to_group();
        // Big wjHdata=skcs.w[tasktp.task_id[j]]*hdata;
        Big wjHdata=pfc->Zpmulti(skcs.w[tasktp.task_id[ND[j]]],hdata);


        tmp1=tmp1+wjHdata;
        // for(int i=0;i<D.size();i++){
        //     // tmp3=tmp3+skcs.w[tasktp.task_id[i]]*wjHdata;
        //     tmp3=tmp3+pfc->Zpmulti(skcs.w[tasktp.task_id[i]],wjHdata);
        // }
    }
    for(int i=0;i<D.size();i++)
        tmp2=tmp2+skcs.w[tasktp.task_id[D[i]]];
    tmp3=pfc->Zpmulti(tmp1,tmp2);

    Big aa=pfc->Zpmulti(skcs.w[0],hmsg);
    Big bb=pfc->Zpmulti(aa,tmp2);
    redsig.w_D1=pfc->mult(pp.w_g,skcs.w[0]*hmsg+tmp1);
    //redsig.w_D2=pfc->mult(signature.w_T1,tmp2)+pfc->mult(pp.w_h,skcs.w[0]*hmsg+tmp1)+pfc->mult(pp.w_g,bb+tmp3);
    redsig.w_D2=pfc->mult(pp.w_g,bb+tmp3);

    redsig.sigma1=signature.sigma1;
    redsig.sigma2=signature.sigma2;
    redsig.sigma3=signature.sigma3;
    redsig.w_T1=signature.w_T1;
    redsig.w_T2=signature.w_T2;
    redsig.T3=signature.T3;
    redsig.pi3=signature.pi3;
    return 0;
}


int Pri_HMTC::RedactVerify(Dataset dataset,Task_tp disclosure,h_Sigma redsig){
    GT p1,p2,p3,p4;

    //p1
    G2 sXDTW=pkcs.w_X1+redsig.w_D1+redsig.w_T2;
    for(int i=0;i<disclosure.task_id.size();i++){
        pfc->start_hash();
        char* data=dataset.data[i];
        pfc->add_to_hash(data);
        Big hdata=pfc->finish_hash_to_group();
        sXDTW=sXDTW+pfc->mult(pkcs.w_W[disclosure.task_id[i]],hdata);
    }

    //p2
    p2=pfc->pairing(pp.w_g,redsig.sigma2);

    //p3
    G1 sW=pkcs.W[disclosure.task_id[0]];
    if(disclosure.task_id.size()>1){
        for(int i=1;i<disclosure.task_id.size();i++){
            sW=sW+pkcs.W[disclosure.task_id[i]];
        }
    }

    p3=pfc->pairing(redsig.w_D1,sW);

    //p4
    p4=pfc->pairing(redsig.w_D2,pp.g);

    

    p1=pfc->pairing(sXDTW,redsig.sigma1);
    if(p1!=p2){
        printf("verification of redacted signature failed 1\n");
        return -1;
    }

    if(p3!=p4){
        printf("verification of redacted signature failed 2\n");
        return -1;
    }

    

    return 0;
}


// int Pri_HMTC::AggVerify(AggData aggdata,AggSigma aggsigma){

//     GT p1,p2,p3,p4;
//     //p2
//     G1 ss2=aggsigma.signatures[0].sigma2;
//     if(aggsigma.signatures.size()>1){
//         for(int k=1;k<aggsigma.signatures.size();k++)
//             ss2=ss2+aggsigma.signatures[k].sigma2;
//     }
//     p2=pfc->pairing(pp.w_g,ss2);

//     //p3
//     G1 sW=pkcs.W[aggdata.disclosure_policy.task_id[0]];
//     if(aggdata.disclosure_policy.task_id.size()>1){
//         for(int i=1;i<aggdata.disclosure_policy.task_id.size();i++){
//             sW=sW+pkcs.W[aggdata.disclosure_policy.task_id[i]];
//         }
//     }
    
//     G2 sD1=aggsigma.signatures[0].w_D1;
//     if(aggsigma.signatures.size()>1){
//         for(int k=1;k<aggsigma.signatures.size();k++)
//             sD1=sD1+aggsigma.signatures[k].w_D1;
//     }

//     p3=pfc->pairing(sD1,sW);

//     //p4
//     G2 sD2=aggsigma.signatures[0].w_D2;
//     if(aggsigma.signatures.size()>1){
//         for(int k=1;k<aggsigma.signatures.size();k++)
//             sD2=sD2+aggsigma.signatures[k].w_D2;
//     }
//     p4=pfc->pairing(sD2,pp.g);

//     //p1
//     G2 sXDTW=pkcs.w_X1+aggsigma.signatures[0].w_D1+aggsigma.signatures[0].w_T2;
//     for(int i=0;i<aggdata.disclosure_policy.task_id.size();i++){
//         pfc->start_hash();
//         char* data=aggdata.datasets[0].data[i];
//         pfc->add_to_hash(data);
//         Big hdata=pfc->finish_hash_to_group();
//         sXDTW=sXDTW+pfc->mult(pkcs.w_W[aggdata.disclosure_policy.task_id[i]],hdata);
//     }

//     p1=pfc->pairing(sXDTW,aggsigma.signatures[0].sigma1);
//     if(aggsigma.signatures.size()>1){
//         for(int k=1;k<aggsigma.signatures.size();k++){
//                 G2 sXDTWk=pkcs.w_X1+aggsigma.signatures[k].w_D1+aggsigma.signatures[k].w_T2;
//                 for(int i=0;i<aggdata.disclosure_policy.task_id.size();i++){
//                     pfc->start_hash();
//                     char* data=aggdata.datasets[k].data[i];
//                     pfc->add_to_hash(data);
//                     Big hdata=pfc->finish_hash_to_group();
//                     sXDTWk=sXDTWk+pfc->mult(pkcs.w_W[aggdata.disclosure_policy.task_id[i]],hdata);
//                 }
//             GT pk=pfc->pairing(sXDTWk,aggsigma.signatures[k].sigma1);
//             p1=p1*pk;
//         }
//     }

//     if(p1!=p2){
//         printf("verification of aggsigma failed 1\n");
//         return -1;
//     }

//     if(p3!=p4){
//         printf("verification of aggsigma failed 2\n");
//         return -1;
//     }

//     printf("verification of aggsigma successful\n");

//     return 0;
// }
int Pri_HMTC::Reveal(MSK msk,h_Sigma redsig,int& idmu,PKMU& pkmu){

    G2 w_R=pfc->mult(redsig.w_T1,msk.xta);
    w_R=-w_R+redsig.w_T2;
    GT p1=pfc->pairing(w_R,pp.g);

    for(int i=0;i<mulist.size();i++){
        GT p2=pfc->pairing(pkcs.w_Y,mulist[i].Yu);
        if(p1==p2){
            idmu=i;
            pkmu=mulist[i];
            return 0;
        }
    }

    printf("dishonest user not found");


    return 0;
}

AS_tp Pri_HMTC::Map(Task_tp& tasktp){
    AS_tp astp;
    for(int i=0;i<=nk;i++)
        astp.ass_id.push_back(0);
    for(int i=0;i<tasktp.task_id.size();i++)
        astp.ass_id[tasktp.task_id[i]]=1;
    return astp;
}

void Pri_HMTC::show_Tasktp(Task_tp tasktp){
    for(int i=0;i<tasktp.task_id.size();i++)
        cout<<tasktp.task_id[i]<<' ';
    cout<<endl;
}

void Pri_HMTC::show_AStp(AS_tp astp){
    for(int i=1;i<=nk;i++)
        cout<<astp.ass_id[i]<<' ';
    cout<<endl;
}

void Pri_HMTC::show_Taskall(Task_all taskall){
    for(int tp=1;tp<=ntp;tp++){
        cout<<"tp="<<setw(2)<<tp<<": ";
        show_Tasktp(taskall.task_all[tp]);
    }
}
void Pri_HMTC::show_AS(AS as){
    for(int tp=1;tp<=ntp;tp++){
        cout<<"tp="<<setw(2)<<tp<<": ";
        show_AStp(as.as[tp]);
    }
}
