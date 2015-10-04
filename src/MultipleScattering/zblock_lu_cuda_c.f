! Notes
! 

      subroutine zblock_lu(a,lda,blk_sz,nblk,ipvt,mp,idcol,k)

      external zblock_lu_cuda_c

!      write (*,*), "Inside zblock_lu_cuda_c.f !!! "

      call zblock_lu_cuda_c
     &     ( a, lda, blk_sz, nblk, ipvt, mp, idcol, k)

      return
      
      end

