#include <string.h>
#include "AtomData.hpp"
#include <hdf5.h>
#include "writeSingleAtomData.hpp"

int writeSingleAtomData_hdf5(hid_t loc_id, AtomData &atom)
{
  printf("WARNING: Writing of hdf5 potentials not implemented yet!\n");

  return -1;
}
