// interleave.h -- A random interleaver
// Todd K. Moon

#ifndef INTERLEAVE_H
#define INTERLEAVE_H

class interleave {
   int size;
   int *pi;
   int *piinv;
public:
   interleave(int in_size, unsigned int seed=0);
   ~interleave() { delete[] pi;  delete[] piinv; };
   void Pi(const double *in, double *out);
   void Pi(const unsigned char *in, unsigned char *out);
   void Pi(const double **in, double **out, int nrow);
   void Piinv(const double *in, double *out);
   void Piinv(const unsigned char *in, unsigned char *out);
   void Piinv(const double **in, double **out, int nrow);
   void PiinvTimesoverlay(const double **in, double **out, int nrow);
};


#endif

/*
Local Variables:
compile-command: "g++ -c interleave.cc"
End:
*/
