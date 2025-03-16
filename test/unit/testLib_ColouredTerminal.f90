program testLib_ColouredTerminal
!   Tests that Lib_ColouredTerminal.f90 is producing the correct text decorators,
!   Failure to do so may cause other tests to fail

   use Lib_ColouredTerminal
   use Lib_UtilsForTests
   implicit none

   character(len=4)         ::      text = "test"
   logical             ::      ok
   character(4 + 9)    ::      rtext, gtext, rstart, gstart, rend, gend, arrow!(len=len_trim(text)+9)
   integer             ::      correctcount

   correctcount = 0
   gtext = colour(LIGHT_GREEN, text)
   rtext = colour(RED, text)
   !arrow=rtext(:1)//gtext(:1)
   rstart = rtext(2:5)
   gstart = gtext(2:5)

   rend = rtext(11:13)
   gend = gtext(11:13)

   gtext = gtext(6:9)
   rtext = rtext(6:9)

   if (rstart == "[31m") correctcount = correctcount + 1
   if (rend == "[0m") correctcount = correctcount + 1
   if (gstart == "[92m") correctcount = correctcount + 1
   if (gend == "[0m") correctcount = correctcount + 1
   if (gtext == "test") correctcount = correctcount + 1
   if (rtext == "test") correctcount = correctcount + 1

   print *, gstart, " ", gtext, " ", gend, " ", rstart, " ", rtext, " ", rend
   !print*,correctcount
   ok = haveAllSubTestsPassed(correctcount, 6)

   call announcePassOrFail(ok)
end program testLib_ColouredTerminal
