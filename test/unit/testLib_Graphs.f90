!------This tests multiple subroutines/function in Lib_Primitives.
!------If these test fails then the user must rerun the failed
!------test verbosely to diagnose the problem.
program testLib_Graphs
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    testLib_Graphs from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
!*    Copyright (C) 2024  James Heath

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
!*      A simple test program to show correct working of Lib_CommandLineArguments
   !   Describe test
   !   Any other details

   use Lib_ColouredTerminal !needed
   use Lib_UtilsForTests
   use Lib_Graphs

   implicit none

   logical             ::      ok
   integer             ::      correctcount
   integer             ::      noofTests = 26 !
   character(len=32)   ::      libName = "testLib_Graphs"
   logical             ::      tempCheck

   type(Line)          ::      singleLineCheck1, singleLineCheck2, conncheckL1, conncheckL2
   type(LineList)      ::      lineArrayCheck1
   type(Point)         ::      singlePointCheck1, singlePointCheck2, conncheckP1, conncheckP2, conncheckP3
   type(PointList)     ::      PointArrayCheck1
   type(graph)         ::      GraphCheck1, GraphCheck2, GraphCheck3, GraphCheck4
   integer             ::      ii, intCheck, tempStore !index, and store results

   correctcount = 0

   !test 1 - Can we make a line with null properties?
   tempCheck = .false.
   singleLineCheck1 = Line_ctor() !make line
   if ((singleLineCheck1%from + singleLineCheck1%to & !fields should sum to zero
        + singleLineCheck1%id + singleLineCheck1%subgraph) == 0) tempCheck = .true.
   tempCheck = tempCheck .and. (.not. (singleLineCheck1%active)) !null line should not be active
   call announceSubTest(libName, "Line_null", 1, noofTests, tempCheck, correctcount)

   !test 2 - Can we make a line with defined properties?
   tempCheck = .false.
   singleLineCheck2 = Line_ctor(1, 2, 3, 4)! should sum to 10
   if ((singleLineCheck2%from + singleLineCheck2%to &
        + singleLineCheck2%id + singleLineCheck2%subgraph) == 10) tempCheck = .true.
   tempCheck = tempCheck .and. singleLineCheck2%active !defined line should be active
   call announceSubTest(libName, "Line_ctor0", 2, noofTests, tempCheck, correctcount)

   !test 3 - Can we make an array of pointers to line dervied types all with null properties?
   tempCheck = .false.
   lineArrayCheck1 = LineList_Ctor(10)
   !print*,"DBG arrayOfLinesCtor0 1",(size(lineArrayCheck1%Lines,dim=1))
   !print*,"DBG arrayOfLinesCtor0 2",(storage_size(lineArrayCheck1%Lines))
   if ((size(lineArrayCheck1%Lines, dim=1) + storage_size(lineArrayCheck1%Lines) &
        + lineArrayCheck1%noofLinesActive + lineArrayCheck1%noofLinesMax) == 180) tempCheck = .true. !10 elemnts of 180 bits each + 10 lines max
   call announceSubTest(libName, "LineList_null", 3, noofTests, tempCheck, correctcount)

   !test 4 - Can we make a point with null properties?
   tempCheck = .false.
   singlePointCheck1 = point_Ctor()
   if ((singlePointCheck1%position(1) + singlePointCheck1%position(2) & !fields + size() should sum to 0
        + singlePointCheck1%id + singlePointCheck1%noofConnex &
        + singlePointCheck1%maxNoofConnex) == 0) tempCheck = .true.
   tempCheck = tempCheck .and. (.not. (singlePointCheck1%active)) !null point should not be active
   call announceSubTest(libName, "Point_null", 4, noofTests, tempCheck, correctcount)

   !test 5 - Can we make a point with defined properties?
   tempCheck = .false.
   singlePointCheck2 = point_Ctor((/1.0d0, 2.0d0/), 3)
   if ((singlePointCheck2%position(1) + singlePointCheck2%position(2) & !fields + size() should sum to 18
        + singlePointCheck2%id + singlePointCheck2%noofConnex &
        + singlePointCheck2%maxNoofConnex + size(singlePointCheck2%connex, dim=1)) == 18) tempCheck = .true.
   tempCheck = tempCheck .and. (singlePointCheck2%active) !point should be active
   call announceSubTest(libName, "Point_ctor0", 5, noofTests, tempCheck, correctcount)

   !test 6 - Can we make an array of pointers to Point dervied types all with null properties?
   tempCheck = .false.
   PointArrayCheck1 = PointList_Ctor(10)
   !print*,"DBG PointList_null 1",(size(PointArrayCheck1%Points,dim=1))
   !print*,"DBG PointList_null 2",(storage_size(PointArrayCheck1%Points))
   if ((size(PointArrayCheck1%Points, dim=1) + storage_size(PointArrayCheck1%Points) &
        + PointArrayCheck1%noofPointsActive + PointArrayCheck1%noofPointsMax) == 852) tempCheck = .true. !10 elemnts of 832 bits each + 10 points max
   call announceSubTest(libName, "PointList_null", 6, noofTests, tempCheck, correctcount)         !There is some padding being added is this compiler depedent?

   !test 7 - can we clone a point?
   tempCheck = .false.
   call clone_point(singlePointCheck2, singlePointCheck1)
   if ((singlePointCheck1%position(1) + singlePointCheck1%position(2) & !fields + size() should sum to 18
        + singlePointCheck1%id + singlePointCheck1%noofConnex &
        + singlePointCheck1%maxNoofConnex + size(singlePointCheck1%connex, dim=1)) == 18) tempCheck = .true.
   tempCheck = tempCheck .and. (singlePointCheck1%active) !point should be active
   call announceSubTest(libName, "clone_point", 7, noofTests, tempCheck, correctcount)

   !test 8 - Can we add a point to an existing point list?
   tempCheck = .false.
   do ii = 1, 12
      call AddPoint(PointArrayCheck1, (/1.0d0, 1.0d0/)) !10 points origionally, so more space will need
   end do                                                  !to be allocated.
   if ((PointArrayCheck1%noofPointsActive + PointArrayCheck1%noofPointsMax) == 32) tempCheck = .true.!maxpoints=20,noofpoints=12
   call announceSubTest(libName, "AddPoint", 8, noofTests, tempCheck, correctcount)

   !test 9 - Can we add a line to an existing line list?
   tempCheck = .false.
   intcheck = 0
   do ii = 1, 11 !one more than max so the reallocation is tested
      call AddLine(lineArrayCheck1, PointArrayCheck1, ii, (ii + 1))!AddLine(LineListIn, vertexIn1, vertexIn2)
      associate (line => lineArrayCheck1%Lines)
         intCheck = intcheck + line(ii)%to + line(ii)%from + line(ii)%subgraph
         if (line(ii)%active .eqv. .true.) intCheck = intCheck + 1
      end associate
   end do
   if (intcheck == 220) tempCheck = .true.
   call announceSubTest(libName, "AddLine", 9, noofTests, tempCheck, correctcount)

   !test 10 - Can we delete an existing point?
   tempCheck = .false.
   call deletePoint(singlePointCheck1)
   if ((singlePointCheck1%position(1) + singlePointCheck1%position(2) & !fields + size() should sum to 0
        + singlePointCheck1%id + singlePointCheck1%noofConnex &
        + singlePointCheck1%maxNoofConnex) == 0) tempCheck = .true.
   tempCheck = tempCheck .and. (.not. (singlePointCheck1%active)) !null point should not be active
   call announceSubTest(libName, "deletePoint", 10, noofTests, tempCheck, correctcount)

   !test 11 - Can we remove an existing point?
   tempCheck = .false.
   call removePoint(PointArrayCheck1, 11)
   associate (delpoint => PointArrayCheck1%Points(11))
      if (PointArrayCheck1%noofPointsActive + delpoint%position(1) + delpoint%position(2) + delpoint%id &
          + delpoint%noofConnex + delpoint%maxNoofConnex == 11) tempCheck = .true.
   end associate
   call announceSubTest(libName, "removePoint", 11, noofTests, tempCheck, correctcount)

   !test 12 - Can we allocate a null graph?
   tempCheck = .false.
   GraphCheck1 = graph_null(10, 10)
   intcheck = 0
   associate (g => GraphCheck1)
      intcheck = size(g%PointList%Points, dim=1) + size(g%LineList%Lines, dim=1) + g%noofSubGraphs + storage_size(g)!should sum to 1237
      !print*,"DBG graph_null test, total is", intCheck
      if (intCheck == 1237) tempCheck = .true.
   end associate
   call announceSubTest(libName, "graph_null", 12, noofTests, tempCheck, correctcount)

   !test 13 - can we find how many lines are active in a graph derived type?
   tempCheck = .false.
   call AddPoint(GraphCheck1%PointList, (/0.0d0, 0.0d0/))!unit square of four points
   call AddPoint(GraphCheck1%PointList, (/0.0d0, 1.0d0/))
   call AddPoint(GraphCheck1%PointList, (/1.0d0, 1.0d0/))
   call AddPoint(GraphCheck1%PointList, (/1.0d0, 0.0d0/))

   call AddLine(GraphCheck1%LineList, GraphCheck1%PointList, 1, 2) !join the dots up
   call AddLine(GraphCheck1%LineList, GraphCheck1%PointList, 2, 3)
   call AddLine(GraphCheck1%LineList, GraphCheck1%PointList, 3, 4)
   call AddLine(GraphCheck1%LineList, GraphCheck1%PointList, 4, 1)
   if (reportNoofActiveLines(GraphCheck1) == 4) tempCheck = .true. !there should be four active lines
   call announceSubTest(libName, "reportNoofActiveLines", 13, noofTests, tempCheck, correctcount)

   !test 14 - can we find how many points are active in a graph derived type?
   tempCheck = .false.
   if (reportNoofActiveVerts(GraphCheck1) == 4) tempCheck = .true. !there should be four active points
   call announceSubTest(libName, "reportNoofActiveVerts", 14, noofTests, tempCheck, correctcount)

   !test 15 - can we find how total active line length in a graph derived type?
   tempCheck = .false.
   if (reportLength(GraphCheck1) == 4.0d0) tempCheck = .true. !unit square should have total line liength of four
   call announceSubTest(libName, "reportLength", 15, noofTests, tempCheck, correctcount)

   !test 16 - can we toggle a line's active attribute?
   tempCheck = .false.
   call toggleLineActive(GraphCheck1%LineList, 5)
   if (GraphCheck1%LineList%Lines(5)%active) tempCheck = .true. !this line should now be active
   call announceSubTest(libName, "toggleLineActive", 16, noofTests, tempCheck, correctcount)
   call toggleLineActive(GraphCheck1%LineList, 5)

   !test 17 - can we toggle a point's active attribute?
   tempCheck = .false.
   call togglePointActive(GraphCheck1%PointList, 5)
   if (GraphCheck1%PointList%Points(5)%active) tempCheck = .true. !this point should now be active
   call announceSubTest(libName, "toggleLineActive", 17, noofTests, tempCheck, correctcount)

   !test 18 - can we connect a point to a line?
   tempCheck = .false.
   conncheckP3 = Point_ctor((/2.0d0, 2.0d0/), 3)
   conncheckL2 = Line_ctor(3, 4, 1, 1)
   call addConnection(conncheckL2, conncheckP3)
   if (conncheckP3%connex(1) == 1) tempCheck = .true. !point 3 should now be connected to line 1
   call announceSubTest(libName, "addConnection", 18, noofTests, tempCheck, correctcount)

   !test 19 - can we check a point is connected to a line?
   tempCheck = .false.
   call addConnection(GraphCheck1%LineList%Lines(1), GraphCheck1%PointList%Points(2))
   tempCheck = checkConnectivity(GraphCheck1%LineList%Lines(1), GraphCheck1%PointList%Points(1), GraphCheck1%PointList%Points(2))
   call announceSubTest(libName, "checkConnectivity", 19, noofTests, tempCheck, correctcount)

   !test 20 - can we define the connectivity between a line and its two vertices,
   !checking that the connection does not already exist?
   tempCheck = .false.
   intCheck = 0
   conncheckP1 = Point_ctor((/0.0d0, 0.0d0/), 1)
   conncheckP2 = Point_ctor((/1.0d0, 1.0d0/), 2)
   conncheckL1 = Line_ctor(1, 2, 1, 1)
   intCheck = conncheckP1%connex(1) + conncheckP1%connex(1) !should sum to 0
   call setConnectivity(conncheckL1, conncheckP2, conncheckP1)
   intCheck = intCheck + conncheckP1%connex(1) + conncheckP1%connex(1) !should now sum to 2
   if (intCheck == 2) tempCheck = .true.
   call announceSubTest(libName, "setConnectivity", 20, noofTests, tempCheck, correctcount)

   !test 21 - can we find the smallest line ID connected to a point?
   tempCheck = .false.
   intCheck = 0
   call findSmallestLineIdOnPoint(GraphCheck1%PointList%Points(6), tempStore)
   intCheck = intCheck + tempStore!should be no connection so 0
   call findSmallestLineIdOnPoint(GraphCheck1%PointList%Points(2), tempStore)
   intCheck = intCheck + tempStore! point 2 is connected to lines 1 and 2, so smallest line ID should be 1
   if (intCheck == 1) tempCheck = .true. !hence intcheck should be 1
   call announceSubTest(libName, "findSmallestLineIdOnPoint", 21, noofTests, tempCheck, correctcount)

   !test 22 - can we find the smallest subgraph ID connected to a point?
   tempCheck = .false.
   intCheck = 0
   call findSmallestSubGraphIdOnPoint(GraphCheck1%LineList, GraphCheck1%PointList%Points(6), tempStore)
   intCheck = intCheck + tempStore!should be no connection so 0
   call findSmallestSubGraphIdOnPoint(GraphCheck1%LineList, GraphCheck1%PointList%Points(2), tempStore)
   intCheck = intCheck + tempStore! point 2 is connected to lines 1 and 2, so smallest subgraph ID should be 1
   if (intCheck == 1) tempCheck = .true. !hence intcheck should be 1
   call announceSubTest(libName, "findSmallestSubGraphIdOnPoint", 22, noofTests, tempCheck, correctcount)

   !test 23 - can we find the subgraphs on a graph?                                                (7)
   tempCheck = .false.!                                                                               |
   intcheck = 0!                                                                                      6
   !add a unconnected subgraph to GraphCheck1                                                       |
   call AddPoint(GraphCheck1%PointList, (/0.0d0, 3.0d0/))!P5 !the points                     (5)--5--(6)--7--(8)
   call AddPoint(GraphCheck1%PointList, (/1.0d0, 3.0d0/))!P6
   call AddPoint(GraphCheck1%PointList, (/1.0d0, 4.0d0/))!P7
   call AddPoint(GraphCheck1%PointList, (/2.0d0, 3.0d0/))!P8                                 (2)--2--(3)
   !        |       |
   call AddLine(GraphCheck1%LineList, GraphCheck1%PointList, 5, 6)!L5 !join the dots up        1       3
   call AddLine(GraphCheck1%LineList, GraphCheck1%PointList, 6, 7)!L6                          |       |
   call AddLine(GraphCheck1%LineList, GraphCheck1%PointList, 6, 8)!L7                         (1)--4--(4)

   associate (Ln => GraphCheck1%LineList%Lines)!alias
      if ((Ln(1)%subgraph + Ln(2)%subgraph + Ln(3)%subgraph + Ln(4)%subgraph) == 10) intCheck = 1! sub graphs should = line ids
      if ((Ln(5)%subgraph + Ln(6)%subgraph + Ln(7)%subgraph) == 18) intCheck = intCheck + 1

      call findSubgraphs(GraphCheck1)

      if ((Ln(1)%subgraph + Ln(2)%subgraph + Ln(3)%subgraph + Ln(4)%subgraph) == 4) intCheck = intCheck + 1!4 lines in sub graph, lowest line in is 1
      if ((Ln(5)%subgraph + Ln(6)%subgraph + Ln(7)%subgraph) == 15) intCheck = intCheck + 1!3 lines in sub graph, lowest line in is 5
      !print*,"DBG TEST findSubgraph A intcheck is:",intCheck
   end associate
   if (intCheck == 4) tempCheck = .true.
   call announceSubTest(libName, "findSubgraphs A", 23, noofTests, tempCheck, correctcount)

   !test 24 - can we find the subgraphs on a graph where the lines are not in a nice order?
   tempCheck = .false.
   intCheck = 0
   GraphCheck2 = graph_null(4, 5)

   do ii = 1, 5
      call AddPoint(GraphCheck2%PointList, (/(ii*1.0d0), 0.0d0/))!points in a line, coords dont actually matter here just connectivity
   end do

   call AddLine(GraphCheck2%LineList, GraphCheck2%PointList, 1, 2) !connect the points so the line positions are not in order of ID
   call AddLine(GraphCheck2%LineList, GraphCheck2%PointList, 3, 4) ! (1)--1--(2)--4--(3)--2--(4)--3--(5)
   call AddLine(GraphCheck2%LineList, GraphCheck2%PointList, 4, 5)
   call AddLine(GraphCheck2%LineList, GraphCheck2%PointList, 2, 3)

   associate (ln => GraphCheck2%LineList%Lines)
      if ((Ln(1)%subgraph + Ln(2)%subgraph + Ln(3)%subgraph + Ln(4)%subgraph) == 10) intCheck = intCheck + 1
      call findSubgraphs(GraphCheck2)
      if ((Ln(1)%subgraph + Ln(2)%subgraph + Ln(3)%subgraph + Ln(4)%subgraph) == 4) intCheck = intCheck + 1

   end associate
   if (intCheck == 2) tempCheck = .true.

   call announceSubTest(libName, "findSubgraphs B", 24, noofTests, tempCheck, correctcount)

   !test 25 - can we find the id of a point on a line given the other.
   tempCheck = .false.

   GraphCheck3 = graph_null(3, 3) !define triangle graph
   call AddPoint(GraphCheck3%PointList, (/0.0d0, 0.0d0/))
   call AddPoint(GraphCheck3%PointList, (/0.0d0, 1.0d0/))
   call AddPoint(GraphCheck3%PointList, (/1.0d0, 0.0d0/))
   call AddLine(GraphCheck3%LineList, GraphCheck3%PointList, 1, 2)
   call AddLine(GraphCheck3%LineList, GraphCheck3%PointList, 2, 3)
   call AddLine(GraphCheck3%LineList, GraphCheck3%PointList, 3, 1)

   if (getLineOtherEnd(GraphCheck3, 3, 3) == 1) tempCheck = .true. !line 3 is connected to points 3 & 1, so inpuit 3 get 1
   call announceSubTest(libName, "getLineOtherEnd", 25, noofTests, tempCheck, correctcount)

   !test 26 - can we merge two points
   tempCheck = .false.
   GraphCheck4 = graph_null(2, 4)
   call AddPoint(GraphCheck4%PointList, (/0.0d0, 0.0d0/))
   call AddPoint(GraphCheck4%PointList, (/1.0d0, 0.25d0/))
   call AddPoint(GraphCheck4%PointList, (/2.0d0, 0.0d0/))
   call AddPoint(GraphCheck4%PointList, (/1.0d0, -0.25d0/))
   call AddLine(GraphCheck4%LineList, GraphCheck4%PointList, 1, 2)
   call AddLine(GraphCheck4%LineList, GraphCheck4%PointList, 4, 3)

   associate (testPoint => GraphCheck4%PointList%Points)
      call mergePoints(testPoint(2), testPoint(4), GraphCheck4, (/1.0d0, 0.0d0/))
      if ((testPoint(2)%connex(2) + testPoint(3)%connex(2) + testPoint(2)%position(1)) == 7.0d0) tempCheck = .true.
   end associate

   call announceSubTest(libName, "mergePoints", 26, noofTests, tempCheck, correctcount)
   !-----------------------------------------------
   ok = haveAllSubTestsPassed(correctcount, noofTests)

   call announcePassOrFail(ok)

end program testLib_Graphs
