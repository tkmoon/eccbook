//
//  polarcode -- encoding and decoding for polar codes
//
//
//  Todd K. Moon
//  Utah State University
//
//  Date: July 31, 2018
//
//  Derived from code provided by Harish Vangala (Monash University)
//    See http://www.ecse.monash.edu.au/eviterbo/polarcodes.html,
//  A translation into C++ from the pseudocode provided in
//     "List Decoding of Polar Codes," by Ido Tal and
//     Alexander Vardy, arXiv :1206.0050v1 31 May 2012
//  For the list decoding, also from Saurabh Tavildar
//    See https://github.com/tavildar/Polar/blob/master/PolarC/PolarCode.cpp
//  who clarifies the getArraypointer functions,
//  and provides LLR decoding

#ifndef POLARCODE_H
#define POLARCODE_H

#include <math.h>
#include "matalloc.h"
#include <string>
#include <iostream>
#include "polarstack.h"
using namespace std;

// types:
//   Encoders:
//   regular encoder    
//   systematic encoder  

//   Decoders
//   SC decoder          

//   list decoder
//   list decoder, llr
//   list decoder, CRC
//   list decoder, systematic
//   list decoder, llr, CRC
//   list decoder, llr, systematic
//   list decoder, CRC, systematic
//   list decoder, llr, systematic, CRC

typedef unsigned char encodertype;
constexpr encodertype nonsystematicenc = 1;
constexpr encodertype systematicenc = 2;
constexpr encodertype withcrcenc = 4;

typedef unsigned char decodertype;
constexpr decodertype SCdec = 1;
constexpr decodertype listdec = 2;
constexpr decodertype listllrdec = 4;
constexpr decodertype withcrcdec = 8;

#include "bittype.h"

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

class polarcode
{
public:
   int N;						// code length
   int K;						// message length
   int n;						// number of levels in tree
   int L;						// number of elements in list
   int crcsize;					// number of bits of CRC parity


   encodertype enctype;
   decodertype dectype;
   
   double channelfact; // -4*srt(Ec)/N0.  Set using setchannelfact(  )

   BITTYPE *uwork;		   // working data for encoder
   BITTYPE *uworksave;	   // saved working data for encoder (for debug)
   
   // bit reverse data:  bitreverser, level1array, level0array:
   //    These could be computed on the fly, but for
   //    speed we'll compute and save as part of the constructor
   int *bitreverser;

   // Data for Systematic encoding
   BITTYPE **Bsys;


   // Data for SC decoding
   int *level1array, *level0array;
   double *LLR;

   // Polar Code Design Data
   SIGNEDBITTYPE *polarcodedesign; // polar code construction:
   // -1 = used bit.   0,1 = frozen bit value

   // SC decoder data
   BITTYPE **BITS;				// propagated bits
   BITTYPE *bitsout;			// output bits
   BITTYPE *messagebitsout;
   double *LLRfinal;			   // save final LLRs (for test/debug)

   // Data for list decoding
   polarstack inactivePathIndices;
   bool *activePath;
   double ****arrayPointer_P; // (non-llr decoding)
   double ***arrayPointer_llr; // llr decoding
   BITTYPE ****arrayPointer_C;
   int **pathIndexToArrayIndex;
   polarstack *inactiveArrayIndices;
   int **arrayReferenceCount;
   BITTYPE *decout;
   double** probForks;
   bool** contForks;
   int *idx;					// used for sorting/selecting
   double *pathmetric;			// LLR decoding

   // Data for CRC decoding (with list decoding)
   BITTYPE **arrayPointer_message;	// Holds the message bits
   // for each element in the list.  This is used
   // for CRC-based list decoding.  If the CRC stuff is not used,
   // this variable can be entirely removed from the code.
   // (In which case, the decoded codeword can be un-encoded
   // to obtain the message bits.)
   bool do_message;
   // if this is set, then this keeps track of the message bits.
   // do_message is set automatically if crc_check > 0.  It can also be
   // set by calling set_do_message().
   BITTYPE *message;  // K - crcsize - bit message

   bool crc_check(BITTYPE* messageblock);

   // channel parameters
   double R, Ec, Ecsqrt, Eb, EbN0, N0, sigma, sigma2;

   // CRC parameters
   int crcp;   // degree of CRC polynomial (16)
   unsigned int crcg; // CRC polynomial
   unsigned int crcmask;  // (1<< crcp) - 1;


   polarcode(int N_in, int Kp_in,
			 encodertype enctype_in=nonsystematicenc,
			 decodertype dectype_in=SCdec,
			 int L_in = 0, int crcsize_in = 0);
// Constructor:
// Kp: size of code without the CRC (if any).  K = Kp + crcsize
// If L_in > 1, then list decoding is assumed.
// If crcsize_in != 0, then crc decoding is assumed (relevant only for list decoding)

   ~polarcode();					 // destructor

   int bitreverse(int x, int n, int& level1,int& level0);

   // Functions for nonsystematic encoding
   BITTYPE* encode(BITTYPE* u);
   void vectFn(BITTYPE *uwork);

   BITTYPE* unencode(BITTYPE *x);



   // Functions for SC decoder
   bool SCdecode(double *y, BITTYPE* &codewordout, BITTYPE* & messagebitsout,
					 encodertype no_op);
   void updateLLR(int i, int level1);
   void updateBITS(BITTYPE latestbit, int i, int level0);
   double LLRupper(double LLR1, double LLR2);
   double LLRlower(BITTYPE bit, double LLR1, double LLR2);
   double logdomain_sum(double x, double y);
   double logdomain_diff(double x, double y);
   void LLRtobit(int j, int i);

   // Functions for systematic encoding and decoding
   BITTYPE *encodesyst(BITTYPE *u); 
   bool decodesyst(double *y, BITTYPE* &codewordout,
					   BITTYPE* &messagebitsout, encodertype nop);
   void initencodesyst(void);  	// initialize data space for systematic encoding
   

   // Functions for List decoding for likelihood and LLR-based formulations
   bool listdecode(double *y, BITTYPE* &cw, BITTYPE* &info, encodertype thisenctype);
   bool listdecodeLLR(double *y, BITTYPE* &cw, BITTYPE* &info,
				 encodertype thisenctype);
   void buildListDecodeDataStructures();
   void initializeListDecodeDataStructures();
   int assignInitialPath();
   int clonePath(int l);
   int clonePathllr(int l);
   void killPath(int l);
   double ** getArrayPointer_P(int lambda, int l);
   double * getArrayPointer_llr(int lambda, int l);
   BITTYPE ** getArrayPointer_C(int lambda, int l);
   BITTYPE ** getArrayPointer_Cllr(int lambda, int l);
   void recursivelyCalcP(int lambda, int phi);
   void recursivelyCalcllr(int lambda, int phi);
   void recursivelyUpdateC(int lambda, int phi);
   void recursivelyUpdateCllr(int lambda, int phi);
   void continuePaths_UnFrozenBit(int phi);
   void continuePaths_UnFrozenBitllr(int phi);
   int findMostProbablePath(bool check_crc, encodertype thisenctype);
   int findMostProbablePathllr(bool check_crc, encodertype thisenctype);
   double W(double y, BITTYPE bit);

   // Functions for CRC checking (for list decoding)
   BITTYPE* encodewithCRC(BITTYPE* u);
   unsigned int CRCenc(BITTYPE *data, int len);
   unsigned int CRCcheck(BITTYPE *data, int len);
   unsigned int CRCpolydiv(BITTYPE *data,int len);
   BITTYPE *encodesystwithCRC(BITTYPE *u);
   
   // miscelleaneous functions
   void set_do_message();
   void extractmessage(BITTYPE *messageout, BITTYPE *messagein);
   void extractmessagefromsyst(BITTYPE *messageout, BITTYPE *messagein);

   // Functions for channel parameters
   void setchannelfact(double cf);
   void computeLLR(double *y, double *llrout);
   void setAWGNchannelparams(double EbN0dB);

   // Functions for designing polar codes
   void designBhatt(double designSNRdB);
   void designBhatt2(double designSNRdB);
   void designMC(double designSNRdB, int M);
   void setdesign(SIGNEDBITTYPE *polarcodedesign);

   // Dump out the auxilliary data
   void dumpdat(void);

};


#endif

/*
Local Variables:
compile-command: "clang++ -c polarcode.cc -std=c++11 -stdlib=libc++"
End:
*/



