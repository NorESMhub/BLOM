      subroutine initime
c
c --- ------------------------------------------------------------------
c --- Initializes variables used for timing. Must be called before use
c --- of 'gettime'
c --- ------------------------------------------------------------------
c
      implicit none
c
#include "timing.h"
c
      call system_clock(end_count,count_rate,count_max)
c
      total_time       =0.
      total_xio_time   =0.
      auxil_total_time =0.
      getflx_total_time=0.
      seaice_total_time=0.
      tmsmt1_total_time=0.
      advdif_total_time=0.
      sfcstr_total_time=0.
      momtum_total_time=0.
      pgforc_total_time=0.
      barotp_total_time=0.
      pbtcor_total_time=0.
      convec_total_time=0.
      diapfl_total_time=0.
      thermf_total_time=0.
      mxlayr_total_time=0.
      tmsmt2_total_time=0.
      diaacc_total_time=0.
      io_total_time    =0.
c
      return
      end
