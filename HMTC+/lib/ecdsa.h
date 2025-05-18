#ifndef ECDSA_H
#define ECDSA_H
#include "acme.h"
#include "macddh.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <cstring>
using namespace std;

struct Keypair
{
    int vk_ep;
    Big vk_x;
    Big sk_d;
};

struct Signature
{
    Big r;
    Big s;
};

class ECDSA
{

public:
    PFC *pfc;
    MACddh mac_ddh;
    ACME acme;
    FAC fac;

    ECDSA(PFC *p);
    ~ECDSA();
    Keypair keygen();
    Signature sign(Big msg_hash, Keypair keypair);
    bool verify(Big msg_hash, Signature signature,Keypair keypair);

};

#endif // PRISRV_H
