#include "calculateChargesPotential.hpp"


void calculateChargesPotential(LSMSCommunication &comm, LSMSSystemParameters &lsms, LocalTypeInfo &local, CrystalParameters &crystal, int chargeSwitch)
{

  Real *qsub;
  Array3d<Real> rhoTemp;

  //qsub = new Real[crystal.num_types];
  qsub = (Real*)shmalloc(crystal.num_types*sizeof(Real));

  for (int i=0; i<crystal.num_types; i++) qsub[i] = 0.0;

  rhoTemp.resize(lsms.global.iprpts+1, 2, local.num_local);
  rhoTemp = 0.0;

  printf("%d:calculateCharges\n",comm.comm.rank);
  calculateCharges(comm, lsms, local, crystal, qsub, rhoTemp, chargeSwitch);  

  // for (int i=0; i<crystal.num_types; i++) printf("i, qsub = %5d %25.15f\n", i, qsub[i]);
  printf("%d:calculatePotential\n",comm.comm.rank);
  calculatePotential(comm, lsms, local, crystal, qsub, rhoTemp, chargeSwitch);
  printf("%d:end of calculatePotential\n",comm.comm.rank);

  //delete[] qsub;
  shfree(qsub);

  return;
}


void calculateLocalCharges(LSMSSystemParameters &lsms, LocalTypeInfo &local, int chargeSwitch)
{

  Array3d<Real> rhoTemp;
  rhoTemp.resize(lsms.global.iprpts+1, 2, local.num_local);

  // Compute integrated densities of states and store in xval**
  // (from mufind_c.f)
  for(int i=0; i<local.num_local; i++)
  {
    if (lsms.n_spin_cant == 2)         // nspin >=3
    {
      local.atom[i].qvalmt = local.atom[i].dosckint[0];
      local.atom[i].qvalws = local.atom[i].dosint[0];
      local.atom[i].mvalws = local.atom[i].dosint[1] * local.atom[i].evecNew[0] + \
                             local.atom[i].dosint[2] * local.atom[i].evecNew[1] + \
                             local.atom[i].dosint[3] * local.atom[i].evecNew[2];
      local.atom[i].mvalmt = local.atom[i].dosckint[1] * local.atom[i].evecNew[0] + \
                             local.atom[i].dosckint[2] * local.atom[i].evecNew[1] + \
                             local.atom[i].dosckint[3] * local.atom[i].evecNew[2];
      local.atom[i].xvalmt[0]    = 0.5 * (local.atom[i].qvalmt + local.atom[i].mvalmt);
      local.atom[i].xvalwsNew[0] = 0.5 * (local.atom[i].qvalws + local.atom[i].mvalws);
      local.atom[i].xvalmt[1]    = 0.5 * (local.atom[i].qvalmt - local.atom[i].mvalmt);
      local.atom[i].xvalwsNew[1] = 0.5 * (local.atom[i].qvalws - local.atom[i].mvalws);
    }
    else if (lsms.n_spin_pola == 2)    // nspin = 2
    {
      local.atom[i].qvalmt = local.atom[i].dosckint[0] + local.atom[i].dosckint[1];
      local.atom[i].qvalws = local.atom[i].dosint[0] + local.atom[i].dosint[1];
      local.atom[i].mvalmt = local.atom[i].dosckint[0] - local.atom[i].dosckint[1];
      local.atom[i].mvalws = local.atom[i].dosint[0] - local.atom[i].dosint[1];

      local.atom[i].xvalmt[0] = local.atom[i].dosckint[0];
      local.atom[i].xvalwsNew[0] = local.atom[i].dosint[0];
      local.atom[i].xvalmt[1] = local.atom[i].dosckint[1];
      local.atom[i].xvalwsNew[1] = local.atom[i].dosint[1];
    }
    else                               // nspin = 1
    {
      local.atom[i].qvalmt = local.atom[i].dosckint[0];
      local.atom[i].qvalws = local.atom[i].dosint[0];
      local.atom[i].mvalmt = 0.0;
      local.atom[i].mvalws = 0.0;

      local.atom[i].xvalmt[0] = local.atom[i].dosckint[0];
      local.atom[i].xvalwsNew[0] = local.atom[i].dosint[0];
    }

    if (lsms.global.iprint > 0)
    {
      printf(" LOCAL WS Int[n(e)] = %18.11f\n", local.atom[i].qvalws);
      printf(" LOCAL MT Int[n(e)] = %18.11f\n", local.atom[i].qvalmt);
      printf(" LOCAL interstial Q = %18.11f\n", local.atom[i].qvalws - local.atom[i].qvalmt);
      
      if (lsms.n_spin_pola == 2)
      {
        printf(" LOCAL WS Int[m(e)] = %18.11f\n", local.atom[i].mvalws);
        printf(" LOCAL MT Int[m(e)] = %18.11f\n", local.atom[i].mvalmt);
      }

      printf(" LOCAL interstial M = %18.11f\n", local.atom[i].mvalws - local.atom[i].mvalmt);
      printf(" Spin = 1 LOCAL WS Int[n(e)] = %18.11f\n", local.atom[i].xvalwsNew[0]);
      printf(" Spin = 2 LOCAL WS Int[n(e)] = %18.11f\n", local.atom[i].xvalwsNew[1]);
      printf(" Spin = 1 LOCAL MT Int[n(e)] = %18.11f\n", local.atom[i].xvalmt[0]);
      printf(" Spin = 2 LOCAL MT Int[n(e)] = %18.11f\n", local.atom[i].xvalmt[1]);
      printf(" Spin = 1 LOCAL interstial Q = %18.11f\n", local.atom[i].xvalwsNew[0] - local.atom[i].xvalmt[0]);
      printf(" Spin = 2 LOCAL interstial Q = %18.11f\n", local.atom[i].xvalwsNew[1] - local.atom[i].xvalmt[1]);
      printf(" LOCAL Moment orientation = (%18.11f, %18.11f, %18.11f)\n", local.atom[i].evecNew[0], local.atom[i].evecNew[1], local.atom[i].evecNew[2]);
    }

 // Calculate qtotws and mtotws

    switch (chargeSwitch)
    {
      case 1:
      {
        local.atom[i].qtotws = local.atom[i].xvalws[0] + \
                               (lsms.n_spin_pola-1) * local.atom[i].xvalws[lsms.n_spin_pola-1] + \
                               local.atom[i].zsemss + \
                               local.atom[i].zcorss;
        local.atom[i].mtotws = local.atom[i].xvalws[0] - local.atom[i].xvalws[lsms.n_spin_pola-1];
        break;
      }
      default:
      {
        local.atom[i].qtotws = local.atom[i].xvalwsNew[0] + \
                               (lsms.n_spin_pola-1) * local.atom[i].xvalwsNew[lsms.n_spin_pola-1] + \
                               local.atom[i].zsemss + \
                               local.atom[i].zcorss;
        local.atom[i].mtotws = local.atom[i].xvalwsNew[0] - local.atom[i].xvalwsNew[lsms.n_spin_pola-1];
      }
    }

//    printf("qtotws = %12.8f\n", local.atom[i].qtotws);
//    printf("mtotws = %12.8f\n", local.atom[i].mtotws);

// from genpot_c.f
//   ================================================================
//   calculate qtotmt and mtotmt.....................................
//   ----------------------------------------------------------------

    Real rSphere;

    switch (lsms.mtasa)
    {
      case 1:
        rSphere = local.atom[i].rws;
        break;
      case 2:
        rSphere = local.atom[i].rws;
        break;
      default:
        rSphere = local.atom[i].rInscribed;
    }

    Real *rTemp;
    rTemp = new Real[local.atom[i].jmt+3];

    rTemp[0] = 0.0;
    for (int j=0; j<local.atom[i].jmt+2; j++)
    {
      //rTemp's indices need to be shifted by 1 for passing into getqm_mt!
      rTemp[j+1] = std::sqrt(local.atom[i].r_mesh[j]);
    }

    switch (chargeSwitch) 
    {
      case 1:
      {
        if (local.atom[i].rhotot(local.atom[i].jmt,0) == 0.0) 
          local.atom[i].rhotot(local.atom[i].jmt,0) = local.atom[i].rhotot(local.atom[i].jmt-1,0);
        if (local.atom[i].rhotot(local.atom[i].jmt+1,0) == 0.0) 
          local.atom[i].rhotot(local.atom[i].jmt+1,0) = local.atom[i].rhotot(local.atom[i].jmt-1,0);
        if (local.atom[i].rhotot(local.atom[i].jmt,lsms.n_spin_pola-1) == 0.0)
          local.atom[i].rhotot(local.atom[i].jmt,lsms.n_spin_pola-1) = local.atom[i].rhotot(local.atom[i].jmt-1,lsms.n_spin_pola-1);
        if (local.atom[i].rhotot(local.atom[i].jmt+1,lsms.n_spin_pola-1) == 0.0) 
          local.atom[i].rhotot(local.atom[i].jmt+1,lsms.n_spin_pola-1) = local.atom[i].rhotot(local.atom[i].jmt-1,lsms.n_spin_pola-1);
    
        if (lsms.global.iprint > 0)
        {
          printf("Spin = 1 rhotot[jmt-1], [jmt], [jmt+1] = %20.16f %20.16f %20.16f\n", local.atom[i].rhotot(local.atom[i].jmt-1,0), local.atom[i].rhotot(local.atom[i].jmt,0), local.atom[i].rhotot(local.atom[i].jmt+1,0));
          printf("Spin = 2 rhotot[jmt-1], [jmt], [jmt+1] = %20.16f %20.16f %20.16f\n", local.atom[i].rhotot(local.atom[i].jmt-1,lsms.n_spin_pola-1), local.atom[i].rhotot(local.atom[i].jmt,lsms.n_spin_pola-1), local.atom[i].rhotot(local.atom[i].jmt+1,lsms.n_spin_pola-1));
        }
    
        getqm_mt_(&lsms.n_spin_pola, &local.atom[i].jmt, &local.atom[i].rInscribed, rTemp, &local.atom[i].rhotot(0,0), &lsms.global.iprpts, &rhoTemp(0,0,i), &lsms.mtasa, &local.atom[i].qtotmt, &local.atom[i].mtotmt, &rSphere, &lsms.global.iprint);
    
        break;
      }
      default:
      {
        if (local.atom[i].rhoNew(local.atom[i].jmt,0) == 0.0) 
          local.atom[i].rhoNew(local.atom[i].jmt,0) = local.atom[i].rhoNew(local.atom[i].jmt-1,0);
        if (local.atom[i].rhoNew(local.atom[i].jmt+1,0) == 0.0) 
          local.atom[i].rhoNew(local.atom[i].jmt+1,0) = local.atom[i].rhoNew(local.atom[i].jmt-1,0);
        if (local.atom[i].rhoNew(local.atom[i].jmt,lsms.n_spin_pola-1) == 0.0)
          local.atom[i].rhoNew(local.atom[i].jmt,lsms.n_spin_pola-1) = local.atom[i].rhoNew(local.atom[i].jmt-1,lsms.n_spin_pola-1);
        if (local.atom[i].rhoNew(local.atom[i].jmt+1,lsms.n_spin_pola-1) == 0.0) 
          local.atom[i].rhoNew(local.atom[i].jmt+1,lsms.n_spin_pola-1) = local.atom[i].rhoNew(local.atom[i].jmt-1,lsms.n_spin_pola-1);
    
        if (lsms.global.iprint > 0)
        {   
          printf("Spin = 1 rhoNew[jmt-1], [jmt], [jmt+1] = %20.16f %20.16f %20.16f\n", local.atom[i].rhoNew(local.atom[i].jmt-1,0), local.atom[i].rhoNew(local.atom[i].jmt,0), local.atom[i].rhoNew(local.atom[i].jmt+1,0));
          printf("Spin = 2 rhoNew[jmt-1], [jmt], [jmt+1] = %20.16f %20.16f %20.16f\n", local.atom[i].rhoNew(local.atom[i].jmt-1,lsms.n_spin_pola-1), local.atom[i].rhoNew(local.atom[i].jmt,lsms.n_spin_pola-1), local.atom[i].rhoNew(local.atom[i].jmt+1,lsms.n_spin_pola-1));
        }
    
        getqm_mt_(&lsms.n_spin_pola, &local.atom[i].jmt, &local.atom[i].rInscribed, rTemp, &local.atom[i].rhoNew(0,0), &lsms.global.iprpts, &rhoTemp(0,0,i), &lsms.mtasa, &local.atom[i].qtotmt, &local.atom[i].mtotmt, &rSphere, &lsms.global.iprint);
      }
    }

    if (lsms.global.iprint >= 0)
    {
      printf("\n");
      printf(" GENPOT / calculateLocalCharges: \n");
      printf(" Total charge and moment in W-S cell:\n");
      printf(" qtotws = %18.11f\n", local.atom[i].qtotws);
      printf(" mtotws = %18.11f\n", local.atom[i].mtotws);
    }

    delete[] rTemp;

  }
}


void calculateCharges(LSMSCommunication &comm, LSMSSystemParameters &lsms, LocalTypeInfo &local, CrystalParameters &crystal, Real *qsub, Array3d<Real> &rhoTemp, int chargeSwitch)
{

  // Compute integrated densities of states and store in xval**
  // (from mufind_c.f)
  for(int i=0; i<local.num_local; i++)
  {
    if (lsms.n_spin_cant == 2)         // nspin >=3
    {
      local.atom[i].qvalmt = local.atom[i].dosckint[0];
      local.atom[i].qvalws = local.atom[i].dosint[0];
      local.atom[i].mvalws = local.atom[i].dosint[1] * local.atom[i].evecNew[0] + \
                             local.atom[i].dosint[2] * local.atom[i].evecNew[1] + \
                             local.atom[i].dosint[3] * local.atom[i].evecNew[2];
      local.atom[i].mvalmt = local.atom[i].dosckint[1] * local.atom[i].evecNew[0] + \
                             local.atom[i].dosckint[2] * local.atom[i].evecNew[1] + \
                             local.atom[i].dosckint[3] * local.atom[i].evecNew[2];
      local.atom[i].xvalmt[0]    = 0.5 * (local.atom[i].qvalmt + local.atom[i].mvalmt);
      local.atom[i].xvalwsNew[0] = 0.5 * (local.atom[i].qvalws + local.atom[i].mvalws);
      local.atom[i].xvalmt[1]    = 0.5 * (local.atom[i].qvalmt - local.atom[i].mvalmt);
      local.atom[i].xvalwsNew[1] = 0.5 * (local.atom[i].qvalws - local.atom[i].mvalws);
    }
    else if (lsms.n_spin_pola == 2)    // nspin = 2
    {
      local.atom[i].qvalmt = local.atom[i].dosckint[0] + local.atom[i].dosckint[1];
      local.atom[i].qvalws = local.atom[i].dosint[0] + local.atom[i].dosint[1];
      local.atom[i].mvalmt = local.atom[i].dosckint[0] - local.atom[i].dosckint[1];
      local.atom[i].mvalws = local.atom[i].dosint[0] - local.atom[i].dosint[1];

      local.atom[i].xvalmt[0] = local.atom[i].dosckint[0];
      local.atom[i].xvalwsNew[0] = local.atom[i].dosint[0];
      local.atom[i].xvalmt[1] = local.atom[i].dosckint[1];
      local.atom[i].xvalwsNew[1] = local.atom[i].dosint[1];
    }
    else                               // nspin = 1
    {
      local.atom[i].qvalmt = local.atom[i].dosckint[0];
      local.atom[i].qvalws = local.atom[i].dosint[0];
      local.atom[i].mvalmt = 0.0;
      local.atom[i].mvalws = 0.0;

      local.atom[i].xvalmt[0] = local.atom[i].dosckint[0];
      local.atom[i].xvalwsNew[0] = local.atom[i].dosint[0];
    }

    if (lsms.global.iprint > 0)
    {
      printf(" LOCAL WS Int[n(e)] = %18.11f\n", local.atom[i].qvalws);
      printf(" LOCAL MT Int[n(e)] = %18.11f\n", local.atom[i].qvalmt);
      printf(" LOCAL interstial Q = %18.11f\n", local.atom[i].qvalws - local.atom[i].qvalmt);
      
      if (lsms.n_spin_pola == 2)
      {
        printf(" LOCAL WS Int[m(e)] = %18.11f\n", local.atom[i].mvalws);
        printf(" LOCAL MT Int[m(e)] = %18.11f\n", local.atom[i].mvalmt);
      }

      printf(" LOCAL interstial M = %18.11f\n", local.atom[i].mvalws - local.atom[i].mvalmt);
      printf(" Spin = 1 LOCAL WS Int[n(e)] = %18.11f\n", local.atom[i].xvalwsNew[0]);
      printf(" Spin = 2 LOCAL WS Int[n(e)] = %18.11f\n", local.atom[i].xvalwsNew[1]);
      printf(" Spin = 1 LOCAL MT Int[n(e)] = %18.11f\n", local.atom[i].xvalmt[0]);
      printf(" Spin = 2 LOCAL MT Int[n(e)] = %18.11f\n", local.atom[i].xvalmt[1]);
      printf(" Spin = 1 LOCAL interstial Q = %18.11f\n", local.atom[i].xvalwsNew[0] - local.atom[i].xvalmt[0]);
      printf(" Spin = 2 LOCAL interstial Q = %18.11f\n", local.atom[i].xvalwsNew[1] - local.atom[i].xvalmt[1]);
      printf(" LOCAL Moment orientation = (%18.11f, %18.11f, %18.11f)\n", local.atom[i].evecNew[0], local.atom[i].evecNew[1], local.atom[i].evecNew[2]);
    }

 // Calculate qtotws and mtotws

    switch (chargeSwitch)
    {
      case 1:
      {
        local.atom[i].qtotws = local.atom[i].xvalws[0] + \
                               (lsms.n_spin_pola-1) * local.atom[i].xvalws[lsms.n_spin_pola-1] + \
                               local.atom[i].zsemss + \
                               local.atom[i].zcorss;
        local.atom[i].mtotws = local.atom[i].xvalws[0] - local.atom[i].xvalws[lsms.n_spin_pola-1];
        break;
      }
      default:
      {
        local.atom[i].qtotws = local.atom[i].xvalwsNew[0] + \
                               (lsms.n_spin_pola-1) * local.atom[i].xvalwsNew[lsms.n_spin_pola-1] + \
                               local.atom[i].zsemss + \
                               local.atom[i].zcorss;
        local.atom[i].mtotws = local.atom[i].xvalwsNew[0] - local.atom[i].xvalwsNew[lsms.n_spin_pola-1];
      }
    }

    if (lsms.global.iprint > 0)
    {
      printf("qtotws = %12.8f\n", local.atom[i].qtotws);
      printf("mtotws = %12.8f\n", local.atom[i].mtotws);
    }

// from genpot_c.f
/*
     ================================================================
     calculate qtotmt and mtotmt.....................................
     ----------------------------------------------------------------
*/
    Real rSphere;

    switch (lsms.mtasa)
    {
      case 1:
        rSphere = local.atom[i].rws;
        break;
      case 2:
        rSphere = local.atom[i].rws;
        break;
      default:
        rSphere = local.atom[i].rInscribed;
    }

    Real *rTemp;
    rTemp = new Real[local.atom[i].jmt+3];

    rTemp[0] = 0.0;
    for (int j=0; j<local.atom[i].jmt+2; j++)
    {
      //rTemp's indices need to be shifted by 1 for passing into getqm_mt!
      rTemp[j+1] = std::sqrt(local.atom[i].r_mesh[j]);
    }

    switch (chargeSwitch) 
    {
      case 1:
      {
        if (local.atom[i].rhotot(local.atom[i].jmt,0) == 0.0) 
          local.atom[i].rhotot(local.atom[i].jmt,0) = local.atom[i].rhotot(local.atom[i].jmt-1,0);
        if (local.atom[i].rhotot(local.atom[i].jmt+1,0) == 0.0) 
          local.atom[i].rhotot(local.atom[i].jmt+1,0) = local.atom[i].rhotot(local.atom[i].jmt-1,0);
        if (local.atom[i].rhotot(local.atom[i].jmt,lsms.n_spin_pola-1) == 0.0)
          local.atom[i].rhotot(local.atom[i].jmt,lsms.n_spin_pola-1) = local.atom[i].rhotot(local.atom[i].jmt-1,lsms.n_spin_pola-1);
        if (local.atom[i].rhotot(local.atom[i].jmt+1,lsms.n_spin_pola-1) == 0.0) 
          local.atom[i].rhotot(local.atom[i].jmt+1,lsms.n_spin_pola-1) = local.atom[i].rhotot(local.atom[i].jmt-1,lsms.n_spin_pola-1);
    
        if (lsms.global.iprint > 0)
        {
          printf("Spin = 1 rhotot[jmt-1], [jmt], [jmt+1] = %20.16f %20.16f %20.16f\n", local.atom[i].rhotot(local.atom[i].jmt-1,0), local.atom[i].rhotot(local.atom[i].jmt,0), local.atom[i].rhotot(local.atom[i].jmt+1,0));
          printf("Spin = 2 rhotot[jmt-1], [jmt], [jmt+1] = %20.16f %20.16f %20.16f\n", local.atom[i].rhotot(local.atom[i].jmt-1,lsms.n_spin_pola-1), local.atom[i].rhotot(local.atom[i].jmt,lsms.n_spin_pola-1), local.atom[i].rhotot(local.atom[i].jmt+1,lsms.n_spin_pola-1));
        }
    
        getqm_mt_(&lsms.n_spin_pola, &local.atom[i].jmt, &local.atom[i].rInscribed, rTemp, &local.atom[i].rhotot(0,0), &lsms.global.iprpts, &rhoTemp(0,0,i), &lsms.mtasa, &local.atom[i].qtotmt, &local.atom[i].mtotmt, &rSphere, &lsms.global.iprint);
    
        break;
      }
      default:
      {
        if (local.atom[i].rhoNew(local.atom[i].jmt,0) == 0.0) 
          local.atom[i].rhoNew(local.atom[i].jmt,0) = local.atom[i].rhoNew(local.atom[i].jmt-1,0);
        if (local.atom[i].rhoNew(local.atom[i].jmt+1,0) == 0.0) 
          local.atom[i].rhoNew(local.atom[i].jmt+1,0) = local.atom[i].rhoNew(local.atom[i].jmt-1,0);
        if (local.atom[i].rhoNew(local.atom[i].jmt,lsms.n_spin_pola-1) == 0.0)
          local.atom[i].rhoNew(local.atom[i].jmt,lsms.n_spin_pola-1) = local.atom[i].rhoNew(local.atom[i].jmt-1,lsms.n_spin_pola-1);
        if (local.atom[i].rhoNew(local.atom[i].jmt+1,lsms.n_spin_pola-1) == 0.0) 
          local.atom[i].rhoNew(local.atom[i].jmt+1,lsms.n_spin_pola-1) = local.atom[i].rhoNew(local.atom[i].jmt-1,lsms.n_spin_pola-1);

        if (lsms.global.iprint > 0)
        {   
          printf("Spin = 1 rhoNew[jmt-1], [jmt], [jmt+1] = %20.16f %20.16f %20.16f\n", local.atom[i].rhoNew(local.atom[i].jmt-1,0), local.atom[i].rhoNew(local.atom[i].jmt,0), local.atom[i].rhoNew(local.atom[i].jmt+1,0));
          printf("Spin = 2 rhoNew[jmt-1], [jmt], [jmt+1] = %20.16f %20.16f %20.16f\n", local.atom[i].rhoNew(local.atom[i].jmt-1,lsms.n_spin_pola-1), local.atom[i].rhoNew(local.atom[i].jmt,lsms.n_spin_pola-1), local.atom[i].rhoNew(local.atom[i].jmt+1,lsms.n_spin_pola-1));
        }

        getqm_mt_(&lsms.n_spin_pola, &local.atom[i].jmt, &local.atom[i].rInscribed, &rTemp[0], &local.atom[i].rhoNew(0,0), &lsms.global.iprpts, &rhoTemp(0,0,i), &lsms.mtasa, &local.atom[i].qtotmt, &local.atom[i].mtotmt, &rSphere, &lsms.global.iprint);

      }
    }

    if (lsms.global.iprint >= 0)
    {
      printf("\n");
      printf(" %d:GENPOT / calculateCharges: \n",comm.comm.rank);
      printf(" Total charge and moment in W-S cell:\n");
      printf(" qtotws = %18.11f\n", local.atom[i].qtotws);
      printf(" mtotws = %18.11f\n", local.atom[i].mtotws);
    }

    delete[] rTemp;

  }

/*
    ================================================================
    calculate interstitial charge, moment & charge density..........
    ================================================================
*/
  // Total interstitial charge or moment components
  static Real qmIntTotal[4];
  for (int i=0; i<4; i++) qmIntTotal[i] = 0.0;

  for(int i=0; i<local.num_local; i++)
    local.atom[i].qInt = local.atom[i].qtotws - local.atom[i].qtotmt;

  switch (lsms.mtasa)
  {
    case 1:             // ASA (not implemented)
      break;

    default:            // Muffin-tin or ASA-MT cases
      for(int i=0; i<local.num_local; i++)
      {
        int i_vdif = 0;       // Temp. fix, should have been read-in and defined somewhere else

        //local.atom[i].mInt = (local.atom[i].mtotws - local.atom[i].mtotmt) * i_vdif;
        if (i_vdif != 0) local.atom[i].mInt = local.atom[i].mtotws - local.atom[i].mtotmt;

        local.atom[i].mIntComponent[0] = local.atom[i].mInt * local.atom[i].evecNew[0];
        local.atom[i].mIntComponent[1] = local.atom[i].mInt * local.atom[i].evecNew[1];
        local.atom[i].mIntComponent[2] = local.atom[i].mInt * local.atom[i].evecNew[2];

        qmIntTotal[0] += local.atom[i].qInt * Real(local.n_per_type[i]);
        qmIntTotal[1] += local.atom[i].mIntComponent[0] * Real(local.n_per_type[i]);
        qmIntTotal[2] += local.atom[i].mIntComponent[1] * Real(local.n_per_type[i]);
        qmIntTotal[3] += local.atom[i].mIntComponent[2] * Real(local.n_per_type[i]);
      }

      // YingWai: should we communicate only qInterstitial if mInterstitial = 0?
      //globalSum(comm, qmIntTotal, 4);
      globalSum_real(comm, qmIntTotal, 4);

      for (int i=0; i<local.num_local; i++)
      {
        local.atom[i].qInt = qmIntTotal[0];

        if (local.atom[i].qInt <= 0.0)
        {
          printf("GENPOT / calculateCharges:: negative interstitial charges: %20.14f\n", local.atom[i].qInt);
          exitLSMS(comm, 99);
        }

        local.atom[i].rhoInt = local.atom[i].qInt / lsms.volumeInterstitial;

        local.atom[i].mIntComponent[0] = qmIntTotal[1];
        local.atom[i].mIntComponent[1] = qmIntTotal[2];
        local.atom[i].mIntComponent[2] = qmIntTotal[3];

        if (lsms.global.iprint >= 0)
        {
          printf("Total interstitial charge = %18.11f\n", local.atom[i].qInt);
          printf("Total interstitial moment = %18.11f\n", local.atom[i].mInt);
        }
      }
  }
 

  for (int i=0; i<local.num_local; i++)
  {
    switch (lsms.n_spin_cant)
    {
      case 2:
        local.atom[i].mInt = std::sqrt(local.atom[i].mIntComponent[0] * local.atom[i].mIntComponent[0] +
                                       local.atom[i].mIntComponent[1] * local.atom[i].mIntComponent[1] +
                                       local.atom[i].mIntComponent[2] * local.atom[i].mIntComponent[2] );
        break;
      default:
        local.atom[i].mInt = local.atom[i].mIntComponent[2];
    }

    if (std::abs(local.atom[i].mInt) < 1.0e-5)
    {
      local.atom[i].mInt = 0.0;
      local.atom[i].mIntComponent[0] = 0.0;
      local.atom[i].mIntComponent[1] = 0.0;
      local.atom[i].mIntComponent[2] = 0.0;
    }

    // YingWai: Need to see if dqSite, momentSite, momentSiteComponent be kept local or moved to Atom class
    Real dqSite;
    Real momentSite;
    Real momentSiteComponent[3];

    if (lsms.mtasa >= 2)       // ASA-MT case (not implemented)
    {

    }
    else                       // Muffin-tin or ASA cases
    {
/*
      ================================================================
      calculate qsub for this atom....................................
      ================================================================
*/
      qsub[local.global_id[i]] = local.atom[i].ztotss - local.atom[i].qtotmt + \
                                 local.atom[i].rhoInt * local.atom[i].omegaMT;
/*
      ================================================================
      calculate the no. of excess electrons and the moment on the site
      ================================================================
*/
      dqSite = local.atom[i].qtotmt + local.atom[i].qInt / Real(crystal.num_atoms) - \
               local.atom[i].ztotss;
      momentSiteComponent[0] = local.atom[i].mtotmt * local.atom[i].evecNew[0] + \
                               local.atom[i].mIntComponent[0] / Real(crystal.num_atoms);
      momentSiteComponent[1] = local.atom[i].mtotmt * local.atom[i].evecNew[1] + \
                               local.atom[i].mIntComponent[1] / Real(crystal.num_atoms);
      momentSiteComponent[2] = local.atom[i].mtotmt * local.atom[i].evecNew[2] + \
                               local.atom[i].mIntComponent[2] / Real(crystal.num_atoms);
    }
  
    momentSite = std::sqrt(momentSiteComponent[0] * momentSiteComponent[0] + \
                           momentSiteComponent[1] * momentSiteComponent[1] + \
                           momentSiteComponent[2] * momentSiteComponent[2]);

    if(local.atom[i].mtotws < 0) momentSite = -momentSite;

  }

/*
  ================================================================
  obtain the qsub from all other nodes............................
  ----------------------------------------------------------------
*/
  //globalSum(comm, qsub, crystal.num_types);
  globalSum_real(comm, qsub, crystal.num_types);

  if (lsms.global.iprint >= 1)
  {
    for (int i=0; i<local.num_local; i++)
    {
      Real meis_h = 0.0;
      Real meis_hh = 0.0;

      printf("qsub for all atoms:\n");
      printf("j, qsub, madmat, qsub*madmat\n");

      for (int j=0; j<crystal.num_types; j++)
      {
        printf("%5d %25.15f %25.15f %25.15f\n", j, qsub[j], local.atom[i].madelungMatrix[j], qsub[i]*local.atom[i].madelungMatrix[j]);

        meis_h += qsub[j];
        meis_hh += qsub[j] * local.atom[i].madelungMatrix[j];
      }

      printf("sum qsub = %25.15f\n", meis_h);
      printf("sum qsub*madmat = %25.15f\n", meis_hh);
    }
  }
}


void calculatePotential(LSMSCommunication &comm, LSMSSystemParameters &lsms, LocalTypeInfo &local, CrystalParameters &crystal, Real *qsub, Array3d<Real> &rhoTemp, int chargeSwitch)
{
/*
  ================================================================
  calculate muffin-tin zero potential and its contribution to the
  Coulomb energy..................................................
  note: lattice constant factor is included in madmat.............
  ================================================================
*/
  static Real vmtSum = 0.0;
  Real vmt = 0.0;
  Real vmt1 = 0.0;
  static Real u0Sum = 0.0;
  Real u0 = 0.0;

  for (int i=0; i<local.num_local; i++)
  {
    // These few lines are from set_u0v0_const.f (may need to make them stand-alone)
    //Real alpha_mad = 0.0;
    //for (int j=0; j<crystal.num_atoms; j++)
    //  alpha_mad += local.atom[i].madelungMatrix[j];
    //alpha_mad *= 2.0 * local.atom[i].omegaWS;

    getvmt(lsms, local.atom[i], crystal, qsub, local.global_id[i], vmt, vmt1, u0);

    //Real madterm = -(vmt1 - alpha_mad * local.atom[i].rhoInt);
    vmtSum += vmt * local.n_per_type[i];
    u0Sum += u0 * local.n_per_type[i];
  }

/*
  =============================================================
  calculate vmt, the electro-static contribution to the muffin-
  tin zero potential through a global sum......................
  vmt1 is the site-dependent constant potential................
  =============================================================
*/
  Real *ro3, *dz;
  ro3 = new Real[local.num_local];
  dz = new Real[local.num_local];

  switch (lsms.mtasa)
  {
    case 1:                            // ASA case
      // not implemented
      break;

    case 2:                            // MT-ASA case
      // not implemented
      break;

    default:                           // Muffin-tin case
      //globalSum(comm, vmtSum);
      globalSum_double(comm, vmtSum);
      //globalSum(comm, u0Sum);
      globalSum_double(comm, u0Sum);

      for (int i=0; i<local.num_local; i++)
      {
        vmt = vmtSum / lsms.volumeInterstitial;
        lsms.u0 = u0Sum;
/*
        ===============================================================
        calculate the exchange-correlation potential related parameters
        ===============================================================
*/
        ro3[i] = std::pow( (4.0*M_PI/3.0) * local.atom[i].rhoInt, -1.0/3.0);
        dz[i] = local.atom[i].mInt / local.atom[i].qInt;

        if (dz[i] < -1.0)
          printf("GENPOT :: dz, mint, qint = %20.12f %20.12f %20.12f\n", \
                 dz[i], local.atom[i].mInt, local.atom[i].qInt);
      }
  }

  for (int i=0; i<local.num_local; i++)
  {

    Real *vrms;
    vrms = new Real[lsms.n_spin_pola];

    int iexch = 0;                  //tmp. fix, should have been passed from input

    Real rSphere;

    switch (lsms.mtasa)
    {
      case 1:
        rSphere = local.atom[i].rws;
        break;
      case 2:
        rSphere = local.atom[i].rws;
        break;
      default:
        rSphere = local.atom[i].rInscribed;
    }

    Real *rTemp;
    rTemp = new Real[local.atom[i].jmt+3];

    rTemp[0] = 0.0;
    for (int j=0; j<local.atom[i].jmt+2; j++)
    {
      //rTemp's indices need to be shifted by 1 for passing into getqm_mt!
      rTemp[j+1] = std::sqrt(local.atom[i].r_mesh[j]);
    }

    local.atom[i].vrNew = 0.0;

    for (int is=0; is<lsms.n_spin_pola; is++)
    {
/*
      =============================================================
      calculate vxcout, the exchange-correlation potential corres-
      ponding to the interstial constant charge density, and excout,
      the exchange-correlation energy.............................. 
      vmtz is the muffin-tin zero potential.
      emad will be used in the total energy calculation.
      emadp will be used in the pressure calculation.
      =============================================================
*/
      Real spin = 1.0 - Real(is) * 2.0;
      local.atom[i].exchangeCorrelationV[is] = 0.0;       // vxcout
      local.atom[i].exchangeCorrelationE = 0.0;           // excout
  
      switch (chargeSwitch) 
      {
        case 1:
        {
          newexchg_(&lsms.n_spin_pola, &spin, &local.atom[i].rhotot(0,0), &local.atom[i].rhotot(0,lsms.n_spin_pola-1), &local.atom[i].exchangeCorrelationPotential(0,is), &local.atom[i].exchangeCorrelationEnergy(0,is), &local.atom[i].exchangeCorrelationV[is], &local.atom[i].exchangeCorrelationE, &ro3[i], &dz[i], &local.atom[i].r_mesh[0], &local.atom[i].jmt, &iexch);
          break;
        }
        default:
        {
          newexchg_(&lsms.n_spin_pola, &spin, &local.atom[i].rhoNew(0,0), &local.atom[i].rhoNew(0,lsms.n_spin_pola-1), &local.atom[i].exchangeCorrelationPotential(0,is), &local.atom[i].exchangeCorrelationEnergy(0,is), &local.atom[i].exchangeCorrelationV[is], &local.atom[i].exchangeCorrelationE, &ro3[i], &dz[i], &local.atom[i].r_mesh[0], &local.atom[i].jmt, &iexch);
        }
      }

      // Energy calculations for the moment
      //if (lsms.mtasa < 1)
      //{
      //   emad(is) = sfac*(qint+sp*mint)*excout
      //   emadp(is)=-sfac*(qint+sp*mint)*three*(excout-vxout(is))
      //}

      // Skipped electron-positron correlation here:
      // (note: if this is put back in, have to get out of the "is" loop first)
      //      if(iexch.ge.100)then
      //        epcorrave=n_per_type*
      //     >  epcorr(rhojmt,0.d0,r_mesh(jmt),0,alpgga,50.d0)
      //        call GlobalSum(epcorrave)
      //        epcorrave=epcorrave/dble(num_atoms)
      //        vxout(1)=epcorrave
      //        vxout(n_spin_pola)=epcorrave
      //      endif

      int isOld;

      local.atom[i].afm = 0;          // tmp. fix, need to calculate it somewhere else 
      switch (local.atom[i].afm)
      {
        case 1:
          isOld = 1 - is;
          break;
        default:
          isOld = is;
      }

/*    =============================================================
      generate new potential.......................................
      -------------------------------------------------------------
*/
      switch (chargeSwitch)
      {
        case 1:
        {
          newpot_(&lsms.n_spin_pola, &local.atom[i].ztotss, &local.atom[i].rhotot(0,0), 
                  &local.atom[i].rhotot(0,lsms.n_spin_pola-1), &rhoTemp(0,0,i), 
                  &local.atom[i].vr(0,isOld), &local.atom[i].vrNew(0,is), &vrms[is], 
                  &local.atom[i].exchangeCorrelationPotential(0,is), &vmt1, &vmt, &local.atom[i].exchangeCorrelationV[is], rTemp, 
                  &local.atom[i].jmt, &local.atom[i].rInscribed, &rSphere, 
                  &lsms.mtasa, &iexch);
          break;
        }
        default:
        {
          newpot_(&lsms.n_spin_pola, &local.atom[i].ztotss, &local.atom[i].rhotot(0,0), 
                  &local.atom[i].rhotot(0,lsms.n_spin_pola-1), &rhoTemp(0,0,i), 
                  &local.atom[i].vr(0,isOld), &local.atom[i].vrNew(0,is), &vrms[is], 
                  &local.atom[i].exchangeCorrelationPotential(0,is), &vmt1, &vmt, &local.atom[i].exchangeCorrelationV[is], rTemp, 
                  &local.atom[i].jmt, &local.atom[i].rInscribed, &rSphere, 
                  &lsms.mtasa, &iexch);
        }
      }

    }

/*
    ================================================================
    vdif is the difference between the spin-down muffin-tin zero and
    the spin-up muffin-tin zero potentials..........................
    Note that vr(ir,1) and vr(ir,2) are all relative to their own
    muffin-tin zero.................................................
    ================================================================
*/
    local.atom[i].vdifNew = local.atom[i].exchangeCorrelationV[lsms.n_spin_pola] - local.atom[i].exchangeCorrelationV[0];

    // shift minority potential by an amount coresponding to -i_vdif Tesla

    int i_vdif = 0;		// tmp. fix, should have been passed in from inputs
    if ( (i_vdif < 0) && (lsms.n_spin_cant != 2) )
      local.atom[i].vdifNew -= Real(i_vdif) * 4.256e-6;

    delete[] rTemp;
    delete[] vrms;
  }

  delete[] ro3;
  delete[] dz;

  return;

}
