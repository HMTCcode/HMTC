#include"priHMTC.h"
#include <tr1/unordered_map>
#include <string>
#include <ctime>
#include <time.h>
#include <iomanip>
#include "lib/aes_ctr.h"
#define TEST_TIME 1
//#define AES_SECURITY 192
#define AES_SECURITY 128
//#define AES_SECURITY 80
using namespace std;

void prihmtc_process(){

    printf("######################################        Execute PriHMTC        ######################################\n\n");
    
    srand(time(NULL));
    PFC pfc(AES_SECURITY);
    Pri_HMTC prihmtc(&pfc);

    int nk=10;
    int ntp=10;
    int ns=10;




    printf("~~~~~~~~~~~~~~~~~~          run Initialization                     ~~~~~~~~~~~~~~~~~~\n");
    PP pp;
    MSK msk;
    SKCS skcs;
    PKCS pkcs;
    prihmtc.Setup(ntp,pp,msk);
    prihmtc.CSKeygen(nk,skcs,pkcs);

    char aes_key_dr[16]={0},aes_iv_dr[8]={0};
    Big kdr;
    pfc.random(kdr);
    memcpy(aes_key_dr,kdr.fn->w,16);
    memcpy(aes_iv_dr,kdr.fn->w+2,8);

    printf("~~~~~~~~~~~~~~~~~~          run Time-bound Key Generation          ~~~~~~~~~~~~~~~~~~\n");

    SKMU skmu;
    PKMU pkmu;
    Task_all assignment;

    Task_tp nulltp;
    assignment.task_all.push_back(nulltp);
    for(int tp=1;tp<=ntp;tp++){ //random task assignment
        Task_tp tasktp;
        for(int i=1;i<=nk;i++){
            if(rand()%2||i==1)
                tasktp.task_id.push_back(i);
        }
        assignment.task_all.push_back(tasktp);
    }
    cout<<"~~~~~~assignment:"<<endl;
    prihmtc.show_Taskall(assignment);

    if(prihmtc.MUKeygen(assignment,skcs,skmu,pkmu)!=0)
        return;
    

    printf("~~~~~~~~~~~~~~~~~~          run Assignment Structure Update          ~~~~~~~~~~~~~~~~~~\n");
    SKMU skmu_pre=skmu,skmu_new;
    PKMU pkmu_pre=pkmu,pkmu_new;
    Task_all assignment_pre=assignment,assignment_new;


    Task_tp nulltp_new;
    assignment_new.task_all.push_back(nulltp_new);
    for(int tp=1;tp<=ntp;tp++){ //random task assignment
        Task_tp tasktp;
        for(int i=1;i<=nk;i++){
            if(rand()%2||i==1)
                tasktp.task_id.push_back(i);
        }
        assignment_new.task_all.push_back(tasktp);
    }
    cout<<"~~~~~~new assignment:"<<endl;
    prihmtc.show_Taskall(assignment_new);

    if(prihmtc.Update(assignment_pre,assignment_new,skcs,skmu_pre,pkmu_pre,skmu_new,pkmu_new)!=0)
        return;

    skmu=skmu_new;
    pkmu=pkmu_new;
    assignment=assignment_new;



    printf("~~~~~~~~~~~~~~~~~~          Heterogeneous Dataset Upload           ~~~~~~~~~~~~~~~~~~\n");
    int ctp=1;
    Dataset dataset;
    for(int i=0;i<assignment.task_all[ctp].task_id.size();i++){
        dataset.data.push_back(strdup("good morning"));
    }
    char* msg="hello world";
    Sigma signature;
    Big t;
    prihmtc.Sign_t(skmu,ctp,assignment.task_all[ctp],dataset,msg,signature,t);


    //encrypt data set
    char aes_key_mu[16]={0},aes_iv_mu[8]={0};
    G2 E,kmu;
    pfc.random(kmu);
    pfc.start_hash();;
    pfc.add_to_hash(kmu);
    Big hkmu=pfc.finish_hash_to_group();
    memcpy(aes_key_mu,hkmu.fn->w,16);
    memcpy(aes_iv_mu,hkmu.fn->w+2,8);

    E=pfc.mult(pkcs.w_Y,t)+kmu;
    AES_CTR aes_ctr;
    aes_ctr.init(aes_key_mu,aes_key_mu);

    Dataset edataset=dataset;
    vector<unsigned int> esize;
    for (int i = 0; i < edataset.data.size(); i++) {
        esize.push_back(strlen(edataset.data[i]));
        aes_ctr.encrypt_data(edataset.data[i],&(esize[i]));
    }

    //decrypt data set
    G2 dkmu;
    dkmu=-pfc.mult(signature.w_T1,skcs.y)+E;
    pfc.start_hash();;
    pfc.add_to_hash(dkmu);
    hkmu=pfc.finish_hash_to_group();
    memcpy(aes_key_mu,hkmu.fn->w,16);
    memcpy(aes_iv_mu,hkmu.fn->w+2,8);
    AES_CTR aes_ctr_d;
    aes_ctr_d.init(aes_key_mu,aes_key_mu);
    for (int i = 0; i < edataset.data.size(); i++) {
        aes_ctr.decrypt_data(edataset.data[i],esize[i]);
    }



    if(prihmtc.Verify(ctp,assignment.task_all[ctp],dataset,msg,signature)!=0)
        return;
    printf("verification of signature successful\n");




    printf("~~~~~~~~~~~~~~~~~~          Data and Signature Redaction           ~~~~~~~~~~~~~~~~~~\n");
    Task_tp disclosure;
    h_Sigma redsig;

    for(int i=0;i<assignment.task_all[ctp].task_id.size();i++){
        if(rand()%2||i==1)
            disclosure.task_id.push_back(assignment.task_all[ctp].task_id[i]);
    }

    if(prihmtc.Redact(skcs,assignment.task_all[ctp],dataset,msg,disclosure,signature,redsig)!=0)
        return;



    if(prihmtc.RedactVerify(dataset,disclosure,redsig)!=0)
        return;
    printf("verification of redacted signature successful\n");

    printf("~~~~~~~~~~~~~~~~~~          Dishonest User Tracing                 ~~~~~~~~~~~~~~~~~~\n");

    int idmu;
    PKMU pkmu2;
    prihmtc.Reveal(msk,redsig,idmu,pkmu2);
    if(pkmu.Yu!=pkmu2.Yu){
        printf("Tracing of dishonest user failed\n");
    }else{
        printf("Tracing of dishonest user successful\n");
    }
        


    return;
}


void prihmtc_time_test(){
    printf("\n\n\n\n######################################        Execute HMTC+ Time Test        ######################################\n\n");
    
    
    srand(time(NULL));
    PFC pfc(AES_SECURITY);
    G1 xx;
    Big yy;
    double tt=0;
    clock_t start,finish;
    for(int i=1;i<=100;i++){
            pfc.random(xx);
            pfc.random(yy);
            start=clock();
            pfc.mult(xx,yy);
            finish=clock();
            tt+=finish-start;

    }
    double sum = (double)(tt/100.0)/(CLOCKS_PER_SEC*TEST_TIME);
    printf("te: time =%f ms\n",sum*1000.0);

    
    Pri_HMTC prihmtc(&pfc);

    int nk=5;
    int ntp=20;
    int nu=1;
    int ns=0;




    printf("~~~~~~~~~~~~~~~~~~          run Initialization                     ~~~~~~~~~~~~~~~~~~\n");
    PP pp;
    MSK msk;
    SKCS skcs;
    PKCS pkcs;

    start=clock();
    prihmtc.Setup(ntp,pp,msk);

    finish=clock();
    sum = (double)(finish-start)/(CLOCKS_PER_SEC*TEST_TIME);
    printf("Setup: time =%f ms\n",sum*1000.0);



    start=clock();
    prihmtc.CSKeygen(nk,skcs,pkcs);
    finish=clock();
    sum = (double)(finish-start)/(CLOCKS_PER_SEC*TEST_TIME);
    printf("CSKeygen: time =%f ms\n",sum*1000.0);

    printf("~~~~~~~~~~~~~~~~~~          run Time-bound Key Generation          ~~~~~~~~~~~~~~~~~~\n");

    vector<SKMU> tbskmu_list;
    vector<PKMU> pkmu_list;
    Task_all assignment;
    Task_tp nulltp;
    Task_tp onetasktp;
        for(int i=1;i<=nk;i++){
            onetasktp.task_id.push_back(i);
    }
    assignment.task_all.push_back(onetasktp);
    for(int tp=1;tp<=ntp;tp++){ //the same task assignment
        assignment.task_all.push_back(onetasktp);
    }
    // cout<<"~~~~~~assignment:"<<endl;
    // prihmtc.show_Taskall(assignment);

    double time_mu=0,time_cs=0;

    start=clock();
    for(int i=0;i<nu;i++){
        SKMU skmu;
        PKMU pkmu;
        double tm=0,tc=0;
        if(prihmtc.MUKeygen_t(assignment,skcs,skmu,pkmu,tm,tc)!=0)
            return;
        time_mu+=tm;
        time_cs+=tc;
        tbskmu_list.push_back(skmu);
        pkmu_list.push_back(pkmu);
    }
    finish=clock();
    sum = (double)(finish-start)/(CLOCKS_PER_SEC*TEST_TIME);
    sum = sum / nu;
    time_mu = time_mu / nu;
    time_cs = time_cs / nu;
    printf("MUKeygen: time =%f ms\n",sum*1000.0);
    printf("MUKeygen-MU: time =%f ms\n",time_mu*1000.0);
    printf("MUKeygen-CS: time =%f ms\n",time_cs*1000.0);

    printf("~~~~~~~~~~~~~~~~~~          Heterogeneous Dataset Upload           ~~~~~~~~~~~~~~~~~~\n");
    int ctp=1;
    Dataset dataset;
    for(int i=0;i<assignment.task_all[ctp].task_id.size();i++){
        dataset.data.push_back("good morning");
    }
    char* msg="hello world";
    vector<Sigma> signature_list;

    double sum1=0,sum2=0;

    start=clock();
    for(int tp=1;tp<=ntp;tp++){
        for(int u=0;u<nu;u++){
            Sigma signature;
            
            prihmtc.Sign(tbskmu_list[u],tp,assignment.task_all[tp],dataset,msg,signature);
            

            // start=clock();
            // if(prihmtc.Verify(tp,assignment.task_all[tp],dataset,msg,signature)!=0)
            //     return;
            // finish=clock();
            // sum2+= (double)(finish-start)/(CLOCKS_PER_SEC*TEST_TIME);
            signature_list.push_back(signature);
        }
    }
    finish=clock();
    sum1= (double)(finish-start)/(CLOCKS_PER_SEC*TEST_TIME);
    ns=signature_list.size();

    sum1 = sum1 / ns;
    printf("Sign: time =%f ms\n",sum1*1000.0);



    start=clock();
    for(int i=0;i<ns;i++){
        int tp=i+1;
        if(prihmtc.Verify(tp,assignment.task_all[tp],dataset,msg,signature_list[i])!=0)
                return;
    }
    finish=clock();
    sum2= (double)(finish-start)/(CLOCKS_PER_SEC*TEST_TIME);



    sum2 = sum2 / ns;
    printf("Verify: time =%f ms\n",sum2*1000.0);





    printf("~~~~~~~~~~~~~~~~~~          Data and Signature Redaction           ~~~~~~~~~~~~~~~~~~\n");
    Task_tp disclosure;
    vector<h_Sigma> redsig_list;

    for(int i=0;i<assignment.task_all[ctp].task_id.size();i++){
        if(i%2==0)
            disclosure.task_id.push_back(assignment.task_all[ctp].task_id[i]);
    }


    start=clock();
    for(int k=0;k<signature_list.size();k++){

        h_Sigma redsig;
        if(prihmtc.Redact(skcs,assignment.task_all[ctp],dataset,msg,disclosure,signature_list[k],redsig)!=0)
            return;
        redsig_list.push_back(redsig);
        
    }
    finish=clock();
    sum = (double)(finish-start)/(CLOCKS_PER_SEC*TEST_TIME);
    sum = sum / ns;
    printf("Redact: time =%f ms\n",sum*1000.0);



    start=clock();
    for (size_t i = 1; i < redsig_list.size(); ++i) {
        h_Sigma redsig = redsig_list[i];
        if(prihmtc.RedactVerify(dataset,disclosure,redsig)!=0)
            return;  
    }
    finish=clock();
    sum = (double)(finish-start)/(CLOCKS_PER_SEC*TEST_TIME);
    sum = sum / ns;
    printf("RedactVerify time =%f ms\n",sum*1000.0);

    printf("~~~~~~~~~~~~~~~~~~          Dishonest User Tracing                 ~~~~~~~~~~~~~~~~~~\n");

    int idmu;
    PKMU pkmu2;

    //int trace_id=rand()%ns;
    int trace_id=0;
    int uid=trace_id%10;

    cout<<"Trace the signature of ID "<<trace_id<<", created by userID "<<uid<<endl;

    start = clock();
    prihmtc.Reveal(msk,redsig_list[trace_id],idmu,pkmu2);
    finish=clock();
    sum = (double)(finish-start)/(CLOCKS_PER_SEC*TEST_TIME);
    printf("Reveal: time =%f ms\n",sum*1000.0);


        


}




int main()
{

    prihmtc_process();
    prihmtc_time_test();
    return 0;
}
