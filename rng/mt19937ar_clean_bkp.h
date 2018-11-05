#ifndef MT19937AR_CLEAN_BKP_H
#define MT19937AR_CLEAN_BKP_H

#define MT19937_N (624)

typedef struct
{
  unsigned long mt[MT19937_N]; 
  int mti;
  unsigned long int mag01[2];
} rnd_strct;
void get_rnd_params(rnd_strct *p);
void set_rnd_params(rnd_strct *p);

#endif
  
  
