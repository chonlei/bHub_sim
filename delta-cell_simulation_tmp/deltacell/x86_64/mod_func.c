#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _caldel_reg(void);
extern void _candel_reg(void);
extern void _catdel_reg(void);
extern void _kadel_reg(void);
extern void _kdrdel_reg(void);
extern void _nadel_reg(void);
extern void _release_reg(void);
extern void _sst_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," caldel.mod");
    fprintf(stderr," candel.mod");
    fprintf(stderr," catdel.mod");
    fprintf(stderr," kadel.mod");
    fprintf(stderr," kdrdel.mod");
    fprintf(stderr," nadel.mod");
    fprintf(stderr," release.mod");
    fprintf(stderr," sst.mod");
    fprintf(stderr, "\n");
  }
  _caldel_reg();
  _candel_reg();
  _catdel_reg();
  _kadel_reg();
  _kdrdel_reg();
  _nadel_reg();
  _release_reg();
  _sst_reg();
}
