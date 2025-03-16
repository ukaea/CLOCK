module Lib_RandomSeed
!---^^^^^^^^^^^^^^^^^^^^^
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    Lib_RandomSeed from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
!*    Copyright (C) 2024  Daniel Mason

!*    This program is free software: you can redistribute it and/or modify
!*    it under the terms of the GNU General Public License as published by
!*    the Free Software Foundation, either version 3 of the License, or
!*    (at your option) any later version.

!*    This program is distributed in the hope that it will be useful,
!*    but WITHOUT ANY WARRANTY; without even the implied warranty of
!*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!*    GNU General Public License for more details.

!*    You should have received a copy of the GNU General Public License
!*    along with this program.  If not, see <https://www.gnu.org/licenses/>.
!*-----------------------------------------------------------------------------------------------------------------------------------
!*      A simple module which sets a random seed.
!*      Also has a Gaussian variate - this module could be extended for more variates.
!*
   use iso_fortran_env
   implicit none
   private

   public      ::      init_random_seed            !   set a seed for built-in pseudorandom generator
   public      ::      get_random_seed             !   get a seed for built-in pseudorandom generator from clock
   public      ::      gaussianVariate             !   returns a Gaussian variate
   public      ::      ran0                        !   minimum standard pseudorandom number
   public      ::      ran0array                   !   many minimum standard pseudorandom numbers

   !---

   interface gaussianVariate
      module procedure gaussianVariate0
      module procedure gaussianVariate1
   end interface

   interface ran0array
      module procedure ran0arrayr
      module procedure ran0arrayi
   end interface

contains
!---^^^^^^^^

   subroutine get_random_seed(seed_out)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      returns the current value of the random seed(s)
      !*      note that the actual length of the seed may vary by processor
      !*      make sure there is sufficient length in the output array to return the answer

      integer, dimension(:), intent(out)                 ::      seed_out

      integer, dimension(:), allocatable    ::      seed
      integer                             ::      nn

      call random_seed(size=nn)
      allocate (seed(nn))
      call random_seed(get=seed)

      nn = min(nn, size(seed_out))
      seed_out = 0
      seed_out(1:nn) = seed(1:nn)

      return
   end subroutine get_random_seed

   subroutine init_random_seed(seed_in)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      set the random seed to be the input value, or the current time on the system clock
      !*
      integer, intent(in), optional         ::      seed_in
      integer                             ::      ii, nn, clock
      integer, dimension(:), allocatable  ::      seed

      call random_seed(size=nn)
      allocate (seed(nn))

      if (present(seed_in)) then
         seed = seed_in + 37*(/(ii - 1, ii=1, nn)/)
      else
         call system_clock(count=clock)
         seed = clock + 37*(/(ii - 1, ii=1, nn)/)
      end if

      call random_seed(put=seed)

      deallocate (seed)
      return
   end subroutine init_random_seed

   !---

   function gaussianVariate0() result(v)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      use Marsaglia polar method to produce a random gaussian variate
      !*      mean = 0, std dev = 1
      real(kind=real64)       ::      v
      real(kind=real64)       ::      g1, g2
      do
         call random_number(g1)
         call random_number(g2)
         g1 = 2*g1 - 1
         g2 = 2*g2 - 1
         v = g1*g1 + g2*g2
         if (v*(1 - v) > 0) exit
      end do
      v = sqrt(-2*log(v)/v)
      v = g1*v
      return
   end function gaussianVariate0

   function gaussianVariate1(n) result(v)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      use Marsaglia polar method to produce n random gaussian variates
      integer, intent(in)      ::      n
      real(kind=real64), dimension(n)       ::      v
      real(kind=real64)       ::      g1, g2, ss
      integer     ::      ii
      if (n == 1) then
         v(1) = gaussianVariate0()
         return
      end if
      do ii = 1, n, 2
         do
            call random_number(g1)
            call random_number(g2)
            g1 = 2*g1 - 1
            g2 = 2*g2 - 1
            ss = g1*g1 + g2*g2
            if (ss*(1 - ss) > 0) exit
         end do
         ss = sqrt(-2*log(ss)/ss)
         v(ii) = g1*ss
         v(ii + 1) = g2*ss
      end do
      if (mod(n, 2) == 1) then
         v(n) = gaussianVariate0()
      end if
      return
   end function gaussianVariate1

   real(kind=real32) function ran0(seed)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      implementation of Park and Miller's "minimal" random number generator

      integer, intent(inout)       ::      seed
      integer, parameter           ::      IA = 16807
      integer, parameter           ::      IM = 2147483647
      real(kind=real32), parameter ::      AM = 1.0/real(IM, kind=real32)   !   just can't shake the conversion warning. But it's fine.
      integer, parameter           ::      IQ = 127773
      integer, parameter           ::      IR = 2836
      integer, parameter           ::      MASK = 123459876

      integer         ::      kk

      seed = ieor(seed, MASK)            !   using a simple mask on input allows for 'accidentally' calling with seed = 0.

      kk = seed/IQ
      seed = IA*(seed - kk*IQ) - IR*kk
      if (seed < 0) seed = seed + IM
      ran0 = AM*seed
      seed = ieor(seed, MASK)

      return
   end function ran0

   subroutine ran0arrayr(seed, ra)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      implementation of Park and Miller's "minimal" random number generator
      integer, intent(inout)       ::      seed
      real(kind=real32), dimension(:), intent(out)      ::      ra
      integer, parameter           ::      IA = 16807
      integer, parameter           ::      IM = 2147483647
      real(kind=real32), parameter ::      AM = 1.0/real(IM, kind=real32)   !   just can't shake the conversion warning. But it's fine.
      integer, parameter           ::      IQ = 127773
      integer, parameter           ::      IR = 2836
      integer, parameter           ::      MASK = 123459876

      integer         ::      ii, kk

      seed = ieor(seed, MASK)            !   using a simple mask on input allows for 'accidentally' calling with seed = 0.

      do ii = 1, size(ra)
         kk = seed/IQ
         seed = IA*(seed - kk*IQ) - IR*kk
         if (seed < 0) seed = seed + IM
         ra(ii) = AM*seed
      end do
      seed = ieor(seed, MASK)

      return
   end subroutine ran0arrayr

   subroutine ran0arrayi(seed, ra)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      implementation of Park and Miller's "minimal" random number generator
      !*      returns 31 bit random integer.
      integer, intent(inout)       ::      seed
      integer, dimension(:), intent(out)      ::      ra
      integer, parameter           ::      IA = 16807
      integer, parameter           ::      IM = 2147483647
      integer, parameter           ::      IQ = 127773
      integer, parameter           ::      IR = 2836
      integer, parameter           ::      MASK = 123459876

      integer         ::      ii, kk

      seed = ieor(seed, MASK)            !   using a simple mask on input allows for 'accidentally' calling with seed = 0.

      do ii = 1, size(ra)
         kk = seed/IQ
         seed = IA*(seed - kk*IQ) - IR*kk
         if (seed < 0) seed = seed + IM
         ra(ii) = seed
      end do
      seed = ieor(seed, MASK)

      return
   end subroutine ran0arrayi

end module Lib_RandomSeed
