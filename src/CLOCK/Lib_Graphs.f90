module Lib_Graphs
   !---^^^^^^^^^^^^^^^^^^
   !*-----------------------------------------------------------------------------------------------------------------------------------
   !*    Lib_Graphs from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
   !*      Module containing derived data types for various geometric primitives to make a graph and constructors to set them up.

   use iso_fortran_env
   implicit none
   private

   public      ::      Line_ctor
   public      ::      LineList_Ctor !may need iterface later
   public      ::      Point_ctor
   public      ::      PointList_Ctor
   public      ::      clone_point
   public      ::      AddPoint
   public      ::      AddLine
   public      ::      deletePoint
   public      ::      removePoint
   public      ::      graph_null
   public      ::      reportNoofActiveLines
   public      ::      reportNoofActiveVerts
   public      ::      reportLength
   public      ::      toggleLineActive
   public      ::      togglePointActive
   public      ::      findSmallestLineIdOnPoint
   public      ::      checkConnectivity
   public      ::      addConnection
   public      ::      setConnectivity
   public      ::      findSubgraphs
   public      ::      findSmallestSubGraphIdOnPoint
   public      ::      getLineOtherEnd
   public      ::      mergePoints

   integer, private, parameter       ::      LIB_GRAPH_CUT = 0     !   indicating a primitive has been cut from the list.

   !DERIVED TYPES
   type, public     ::      Point
      !---    defines a point in space connected to one or more lines
      integer                             ::      id
      real(kind=real64), dimension(2)      ::      position                   !   position of the vertices. Stored including offset, ie (1:Nx,1:Ny)
      integer                             ::      noofConnex          !   number of lines connected
      integer                             ::      maxNoofConnex       !   max number of lines space allocated
      integer, dimension(:), pointer        ::      connex              !   (1:nConnex) lines connected
      logical                             ::      active              !   Is this point currently used?
   end type

   type, public     ::      Line
      !---    defines a line between two vertices which can be drawn to reconstruct the image
      !private !keep an eye on this, does it need to be private?
      integer                     ::      id, subgraph
      integer                     ::      from, to         !   from and to vertices
      logical                     ::      active          !   Is this line currently used?
   end type

   type, public     ::      LineList
      !---    defines a list of line derived types
      type(Line), dimension(:), pointer                     ::      Lines               !   (0:nLine0) origionally, will start from 1
      integer                                             ::      noofLinesMax        !   Maximum nuber of lines allocated.
      integer                                             ::      noofLinesActive     !   number of lines currently used
   end type

   type, public     ::      PointList
      !---    defines a list of line derived types
      type(Point), dimension(:), pointer                    ::      Points              !   (0:nLine0) origionally, will start from 1
      integer                                             ::      noofPointsMax       !   Maximum nuber of lines allocated.
      integer                                             ::      noofPointsActive    !   number of lines currently used
   end type

   type, public     ::      graph
      !---    Defines and undirected graph
      type(LineList)                                      ::      LineList            !   holds line connectivty
      type(PointList)                                     ::      PointList           !   Hold graph vertices
      integer                                             ::      noofSubGraphs       !   Number of subgraphs
   end type

   interface Line_ctor
      module procedure Line_null
      module procedure Line_ctor0
   end interface

   interface LineList_Ctor
      module procedure LineList_null
      !module procedure LineListPopulate
      !module procedure LineList_ctor0
   end interface

   interface Point_ctor
      module procedure Point_null
      module procedure Point_ctor0
   end interface

   interface PointList_Ctor
      module procedure PointList_null
      !module procedure LineListPopulate
      !module procedure LineList_ctor0
   end interface

contains
   !---^^^^^^^^
   !FUNCTIONS

   function Line_null() result(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !Initialises a line derived type with null properties or marks it as inactive (no deallocation), better to have 0 for nul -1 for ignore?
      type(Line)           ::      this! input line
      this%from = LIB_GRAPH_CUT !initalise as 0
      this%to = LIB_GRAPH_CUT
      this%id = LIB_GRAPH_CUT
      this%subgraph = LIB_GRAPH_CUT
      this%active = .false.
      return
   end function Line_null

   function Line_ctor0(from, to, id, sub) result(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !Sets up a line derived type
      integer, intent(in)              ::      from, to !vertex IDs
      integer, intent(in)              ::      id, sub      !Line ID & sub graph
      type(Line)                 ::      this         !Line object

      this%from = from    !assing members
      this%to = to
      this%id = id
      this%subgraph = sub
      this%active = .true.
      return
   end function Line_ctor0

   function LineList_null(noofLines) result(LineListOut)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !Allocates LineList derived type (essentially array of pointers to line derived types, nested derived types needed
      !to make this work in fortran) with null properties
      integer, intent(in)                              ::      noofLines !number of lines
      type(LineList)               ::        LineListOut
      integer                                         ::      ii !line index
      !allocate(lineArrayOut(noofLines)) !old one started at 0, why? (0:this%nLine0)
      allocate (LineListOut%Lines(noofLines))
      do ii = 1, noofLines                         !for noofLines Lines
         LineListOut%Lines(ii) = Line_ctor()     !Contruct null line
      end do
      LineListOut%noofLinesMax = noofLines          !update non-array fields
      LineListOut%noofLinesActive = 0

      return
   end function LineList_null

   function PointList_null(noofPoints) result(pointListOut)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !Allocates PointList derived type (essentially array of pointers to point derived types, nested derived types needed
      !to make this work in fortran) with null properties
      integer, intent(in)                              ::      noofPoints !number of lines
      type(PointList)               ::        pointListOut
      integer                                         ::      ii !point index
      allocate (pointListOut%Points(noofPoints))
      do ii = 1, noofPoints                        !for noofPoints points
         pointListOut%Points(ii) = Point_ctor()  !construct null point
      end do
      pointListOut%noofPointsMax = noofPoints       !update non-array members
      pointListOut%noofPointsActive = 0

      return
   end function PointList_null

   function Point_null() result(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !sets up a point derived type with null properties
      type(Point)           ::      this !point derived type
      this%position = 0
      this%id = LIB_GRAPH_CUT     !id = 0 for unused
      this%noofConnex = 0
      this%maxNoofConnex = 0
      nullify (this%connex)
      this%active = .false.

      return
   end function Point_null

   function Point_ctor0(p, id) result(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !Defines a point at postion p, with given ID and subimage coords
      real(kind=real64), dimension(2), intent(in)   ::      p       !exact position of point
      integer, intent(in)                          ::      id      !Point ID
      type(Point)                 ::      this                   !point derived type
      this%position = p       !assign members
      this%id = id
      this%noofConnex = 0
      this%maxNoofConnex = 6
      allocate (this%connex(1:this%maxNoofConnex))
      this%connex = 0 !initalise all element to 0
      this%active = .true.
      return
   end function Point_ctor0

   function graph_null(noofLines, noofPoints) result(GraphOut)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !Allocates Graph derived type with a null lineList and null pointList.
      !number of subgraphs initialised as 1.
      integer, intent(in)      ::      noofLines !number of lines
      integer, intent(in)      ::      noofPoints !number of points
      type(Graph)             ::      GraphOut

      GraphOut%noofSubGraphs = 1 !1 before calculation.
      GraphOut%PointList = PointList_null(noofPoints)   !allocate point list
      GraphOut%LineList = LineList_null(noofLines)      !allocate line list
   end function graph_null

   pure integer function reportNoofActiveLines(graphIn)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !returns the number of active lines in a graph derived type.
      type(graph), intent(in)      ::  graphIn
      reportNoofActiveLines = graphIn%LineList%noofLinesActive
   end function reportNoofActiveLines

   pure integer function reportNoofActiveVerts(graphIn)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !returns the number of lines in a graph derived type.
      type(graph), intent(in)      ::  graphIn
      reportNoofActiveVerts = graphIn%PointList%noofPointsActive
   end function reportNoofActiveVerts

   real(kind=real64) function reportLength(graphIn)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !return the total length of all active lines in a graph derived type
      type(graph), intent(in)      ::      graphIn

      integer                     ::      ii  !line index
      !real(kind=real64)           ::      lengthTotal
      type(Line), pointer          ::      currentLine
      type(Point), pointer         ::      pointTo, pointFrom !both ends of the line

      do ii = 1, reportNoofActiveLines(graphIn)
         currentLine => graphIn%LineList%Lines(ii)!shorten expression with pointers
         pointTo => graphIn%PointList%Points(currentLine%to)
         pointFrom => graphIn%PointList%Points(currentLine%from)

         reportLength = reportLength + norm2(pointTo%position - pointFrom%position)!sum up lengths between each line
      end do

   end function reportLength

   logical function checkConnectivity(lineIn, point1In, point2In)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !Returns true if points 1 & 2 are connected by a line, prints warning if
      !the connection is only one way
      type(Line), intent(inout)       ::      lineIn
      type(Point), intent(inout)      ::      point1In
      type(Point), intent(inout)      ::      point2In

      logical                         ::      OK1, OK2 !is point n connected to the line
      integer                         ::      conIdx !connection index

      !check if point 1 is connected to the line already
      OK1 = .false.
      do conIdx = 1, point1In%noofConnex
         if (point1In%connex(conIdx) == 0) CYCLE !0 connection is unset, ignore
         if (lineIn%id == point1In%connex(conIdx)) OK1 = .true.
      end do

      !check if point 2 is connected to the line already
      OK2 = .false.
      do conIdx = 1, point2In%noofConnex
         if (point2In%connex(conIdx) == 0) CYCLE !0 connection is unset, ignore
         if (lineIn%id == point2In%connex(conIdx)) OK2 = .true.
      end do
      !one way connections shouold not happen, warn if it does.
   if (XOR(OK1, OK2)) print *, "WARN| checkConnectivity, one way connection in Lid, Pid1, Pid2", lineIn%id, point1In%id, point2In%id
      checkConnectivity = (OK1 .and. OK2) !proper two way connection detected
   end function checkConnectivity

   integer function getLineOtherEnd(graphIn, lineIdIn, pointIdIn)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !For a given end point (with id pointIdIn) on a line (with id lineIdIn), returns
      !The id of the other end point
      type(Graph), intent(in)     ::      graphIn
      integer, intent(in)         ::      lineIdIn, pointIdIn

      associate (lineIn => graphIn%LineList%Lines(lineIdIn), &        !aliases for clarity
                 pointKnown => graphIn%PointList%Points(pointIdIn))
         if (lineIn%from == pointKnown%id) then  !If the point id in is the same as from
            getLineOtherEnd = lineIn%to         !then the other end is to.
         else if (lineIn%to == pointKnown%id) then   !If the point id in is the same as to
            getLineOtherEnd = lineIn%from           !then the other end is from.
         else    !Otherwise the ponit ID in does not exist on this line
            print *, "WARN | getLineOtherEnd, pointIdIn", pointIdIn, "is not on lineIdIn", lineIdIn
         end if
      end associate
   end function getLineOtherEnd

   !SUBROUTINES

   subroutine clone_point(pointIn, cloneInout)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !Clones pointIn, copying all its field to cloneInout
      type(Point), intent(in)              ::      pointIn !Point dervied types
      type(Point), intent(inout)           ::      cloneInout
      cloneInout%active = pointIn%active
      cloneInout%id = pointIn%id
      cloneInout%position = pointIn%position
      cloneInout%noofConnex = pointIn%noofConnex
      cloneInout%maxNoofConnex = pointIn%maxNoofConnex
      allocate (cloneInout%connex(1:pointIn%maxNoofConnex))
      cloneInout%connex(1:pointIn%maxNoofConnex) = pointIn%connex(1:pointIn%maxNoofConnex)!copy all connections
   end subroutine clone_point

   subroutine AddPoint(PointListIn, p)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !Adds a point to a point list, giving it the next avaiable ID, and checking that ther is sufficient memory to hold it
      !The point list is reallocated if additional memory is required.
      type(PointList), intent(inout)                  ::      PointListIn
      real(kind=real64), dimension(2), intent(in)   ::      p       !exact position of the added point

      type(Point), dimension(:), pointer            ::      tempPoints   !temporary array of pointers to point derived types
      integer                         ::      id, ii      !Point ID, point index

      id = PointListIn%noofPointsActive + 1                   !one more point
      if (id > PointListIn%noofPointsMax) then              !if that exceeds maximum allowed noof points
         !need to allocate more space
         allocate (tempPoints(PointListIn%noofPointsMax*2))   !allocate temporary points list with twice the size
         do ii = 1, (PointListIn%noofPointsActive)              !For all existing ponits
            call clone_point(PointListIn%Points(ii), tempPoints(ii)) !copy over to the temp list
         end do
         deallocate (PointListIn%Points)  !don't need orional any more
         PointListIn%Points => tempPoints  !point to the temporary to replace
         PointListIn%noofPointsMax = PointListIn%noofPointsMax*2 !update maximum number of allowed points
      end if
      PointListIn%Points(id) = Point_ctor(p, id) !add the new point at next available index. !may put some code to find what subimage in here.
      PointListIn%noofPointsActive = PointListIn%noofPointsActive + 1 !update number of active points.

      return
   end subroutine AddPoint

   subroutine AddLine(LineListIn, pointListIn, vertexIn1, vertexIn2)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !Adds a line to a line list, giving it the next avaiable ID, and checking that ther is sufficient memory to hold it
      !The point list is reallocated if additional memory is required, this assumes your points are already allocated and defined.
      type(LineList), intent(inout)                    ::      LineListIn
      type(PointList), intent(inout)                   ::      pointListIn
      integer, intent(in)                              ::      vertexIn1, vertexIn2     !from and to point IDs

      type(Line), dimension(:), pointer            ::      tempLines   !temporary array of pointers to line derived types
      integer                                         ::      lineId, ii ! line id, index

      if (.not. (pointListIn%Points(vertexIn1)%active)) print *, "WARN | AddLine in Lib_Graphs, point", vertexIn1, "is not active"
      if (.not. (pointListIn%Points(vertexIn2)%active)) print *, "WARN | AddLine in Lib_Graphs, point", vertexIn2, "is not active"

      lineId = LineListIn%noofLinesActive + 1                 !one more line
      if (lineId > LineListIn%noofLinesMax) then             !if adding one more line exceeds the maximum allowed number of lines
         allocate (tempLines(LineListIn%noofLinesMax*2))  !allocate a temporary line list twice as large as the origional.
         do ii = 1, (LineListIn%noofLinesMax)               !For the origional lines
            tempLines(ii) = LineListIn%Lines(ii)          !copy across to temporoary line list
         end do
         do ii = (LineListIn%noofLinesMax + 1), (LineListIn%noofLinesMax*2) !for all new lines after the newly added one
            tempLines(ii) = Line_ctor()                                 !construct a line with null properties
         end do
         deallocate (LineListIn%Lines)                        !Dont need origonal, deallocate it.
         LineListIn%Lines => tempLines                         !Point to the temporay compy to replace.
         LineListIn%noofLinesMax = LineListIn%noofLinesMax*2   !update the maximum allowed number of lines
      end if

      LineListIn%noofLinesActive = lineId !update the number of acitve lines
      LineListIn%Lines(lineId) = Line_ctor0(vertexIn1, vertexIn2, lineId, lineId) !Assing the new lines at the nex available index.

      !call addConnection(LineListIn%Lines(lineId),pointListIn%Points(vertexIn1))! connect the line with both vertices
      !call addConnection(LineListIn%Lines(lineId),pointListIn%Points(vertexIn2))
      call setConnectivity(LineListIn%Lines(lineId), pointListIn%Points(vertexIn1), pointListIn%Points(vertexIn2))! connect the line with both vertices

   end subroutine AddLine

   subroutine removePoint(pointListIn, pIdIn)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !Removes a point from a point list
      !and calls delete to free any allocated memory
      type(PointList), intent(inout)       ::      pointListIn !Point list derived type.
      integer, intent(in)                  ::      pIdIn !ID of the point to be removed from the list.

      pointListIn%noofPointsActive = pointListIn%noofPointsActive - 1 !one less active point.
      call deletePoint(pointListIn%Points(pIdIn))

   end subroutine removePoint

   subroutine deletePoint(pointIn) !DEPRECIATED
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !Frees any allocated memory associate with a point derived type.
      type(Point), intent(inout)       ::      pointIn

      if (pointIn%maxNoofConnex == 0) return
      deallocate (pointIn%connex)
      pointIn = Point_null()
      return
   end subroutine deletePoint

   subroutine findSmallestLineIdOnPoint(pointIn, LidOut)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !For a given point finds the smallest line ID connected to it
      type(point), intent(in)      ::  pointIn !Point dervied type
      integer, intent(inout)       ::  LidOut

      integer                     ::  conIdx !connection index

      if (pointIn%noofConnex > 0) then !any lines connected?
         LidOut = pointIn%connex(1)
         do conIdx = 1, pointIn%noofConnex !for all connections
            if (pointIn%connex(conIdx) < LidOut) then
               LidOut = pointIn%connex(conIdx)   !if the line ID is smallest then set LidOut to it
            end if
         end do
      else
         LidOut = 0 !if there are no connections
      end if

   end subroutine findSmallestLineIdOnPoint

   subroutine findSmallestSubGraphIdOnPoint(lineListIn, pointIn, SGidOut)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !For a given point finds the smallest subgraph (of a line) ID connected to it
      type(LineList)              ::  lineListIn
      type(point), intent(in)      ::  pointIn !Point dervied type
      integer, intent(inout)       ::  SGidOut !subgraph id out

      integer                     ::  conIdx !connection index

      if (pointIn%noofConnex > 0) then !any lines connected?
         SGidOut = lineListIn%Lines(pointIn%connex(1))%subgraph
         do conIdx = 1, pointIn%noofConnex !for all connections
            if (lineListIn%Lines(pointIn%connex(conIdx))%subgraph < SGidOut) then
               SGidOut = lineListIn%Lines(pointIn%connex(conIdx))%subgraph   !if the subgraph ID is smallest then set SGidOut to it
            end if
         end do
      else
         SGidOut = 0 !if there are no connections
      end if

   end subroutine findSmallestSubGraphIdOnPoint

   subroutine toggleLineActive(LinelistIn, lineId)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !toggles the active attribute of the line within the linelist selected by lineId
      type(LineList), intent(inout)    ::  LinelistIn
      integer, intent(in)             ::  lineId

      LinelistIn%Lines(lineId)%active = .not. (LinelistIn%Lines(lineId)%active)!toggle active
   end subroutine toggleLineActive

   subroutine togglePointActive(PointlistIn, PointId)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !toggles the active attribute of the point within the pointlist selected by pointId
      type(PointList), intent(inout)   ::  PointlistIn
      integer, intent(in)             ::  PointId

      PointlistIn%Points(PointId)%active = .not. (PointlistIn%Points(PointId)%active)!toggle active
   end subroutine togglePointActive

   subroutine addConnection(lineIn, pointIn)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !adds the line to the list of connected lines to a point. A check for existing connectivity
      !(checkConnectivity) should be called before this. If more space is required to store the
      !connection, the connection array is reallocated and size fields updated.
      !NB: assigning to and from point IDs to the line is handled by addLine
      !modifies point%connex, point%noofConnex and (potentially) point%maxNoofConnex
      type(Line), intent(in)          ::      LineIn
      type(Point), intent(inout)      ::      pointIn

      integer                         ::      tempConCount, conIdx !temporary connection counter, connection index
      integer, dimension(:), pointer    ::      tempConnex   !temporary (1:nConnex) lines connected

      tempConCount = pointIn%noofConnex + 1                   !one more connection

      if (tempConCount > pointIn%maxNoofConnex) then         !is there space for it
         allocate (tempConnex(pointIn%maxNoofConnex*2))   !if not expand

         do conIdx = 1, pointIn%noofConnex
            tempConnex(conIdx) = pointIn%connex(conIdx)   !copy the origional values
         end do

         deallocate (pointIn%connex)                      !deallocate origional
         pointIn%connex => tempConnex                    !replace origional with temporary
         pointIn%maxNoofConnex = pointIn%maxNoofConnex*2   !update max allowed connections
      end if

      pointIn%connex(tempConCount) = LineIn%id              !add connection
      pointIn%noofConnex = tempConCount                     !update noof connections

   end subroutine addConnection

   subroutine setConnectivity(lineIn, point1In, point2In)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !sets the connectivity of points 1 & 2 when they are connected by a line,
      !called at the end of AddLine, allocates more space for connections if needed
      type(Line), intent(inout)       ::      LineIn
      type(Point), intent(inout)      ::      point1In
      type(Point), intent(inout)      ::      point2In

      !check if points are connected to this line already, warn if connection is one way
      if (checkConnectivity(lineIn, point1In, point2In)) return
      !point1
      call addConnection(lineIn, point1In)!add connection at next aviaiable index
      !and update storage field
      !point2
      call addConnection(lineIn, point2In)

   end subroutine setConnectivity

   subroutine findSubgraphs(graphIn)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !finds what subgraph each line is connected to and updates datastructure accordingly.
      type(Graph), intent(inout)       ::      graphIn

      integer                     ::      actLineCount, lineCount, iterationCount!, changeCount
      integer                     ::      currentLowestSGId, trialId
      logical                     ::      change !has a chane happened

      actLineCount = 0
      lineCount = 0
      iterationCount = 0
      currentLowestSGId = 1
      iterationCount = 0
      change = .true.

      !iterate through graphIn%lineList%line
      do while ((change .eqv. .true.) .or. (iterationCount < graphIn%LineList%noofLinesActive)) !Shouldnt have to do more iterations
         iterationCount = iterationCount + 1 !1,2,3,4 slowly getting faster...               !than active lines. Stop when no changes.
         !changeCount=0
         lineCount = 0
         actLineCount = 0
         do while (actLineCount < graphIn%LineList%noofLinesActive)    !new iteration when all active lines have been checked.
            lineCount = lineCount + 1                                   !need to keep track of this if there are inactive lines
            change = .false.
            !only considers active lines (graphIn%lineList%line(:)%active==true)
            if (graphIn%LineList%Lines(lineCount)%active .eqv. .true.) then !only look at active lines
               actLineCount = actLineCount + 1
               associate (lineTmp => graphIn%LineList%Lines(lineCount))
                  currentLowestSGId = lineTmp%subgraph
                  !check from point
                  associate (pointTmp => graphIn%PointList%Points(lineTmp%from))   !alias for readability
                     if (pointTmp%noofConnex > 0) then !only look at points that have connections
                        call findSmallestSubGraphIdOnPoint(graphIn%LineList, pointTmp, trialId)   !look at all lines connected
                        if (trialId < currentLowestSGId) then                                     !to a point, which has the lowest
                           currentLowestSGId = trialId                                           !subgraph id, return min value.
                           change = .true.
                        end if
                     end if
                  end associate
                  !check to point
                  associate (pointTmp => graphIn%PointList%Points(lineTmp%to))     !same again for the line's "to" point.
                     if (pointTmp%noofConnex > 0) then
                        call findSmallestSubGraphIdOnPoint(graphIn%LineList, pointTmp, trialId)
                        if (trialId < currentLowestSGId) then
                           currentLowestSGId = trialId
                           change = .true.
                        end if
                     end if
                  end associate
                  graphIn%LineList%Lines(lineCount)%subgraph = currentLowestSGId    !set the lines subgraph to it own ID or
                  !the smallest line ID it is connected to.
               end associate
            end if
         end do
      end do

   end subroutine findSubgraphs

   subroutine mergePoints(keepPoint, hidePoint, graphIn, newPosition)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !merges two points by moving keepPoint (KP) to (newPosition) and setting hidePoint to
      !inactive. Lines connected to hidePoint (HP) are duplicated, swapping their connectivity with
      !keepoint. The original lines are then set to inactive
      ! Before :(P1)--1--(KP)             |After:  (P1)--1--(KP)--3--(P2)
      !                  (HP)--2--(P2)    |                           /
      !                                   |                          2
      !                                   |                         /
      !                                   |                      {HP}
      !where: -- active line, / inactive line, (PN) active point, {PN} inactive point
      !More than one line can be conected to the input points. Point/line active fields are updated
      !NB: this is a simplistic but memory inefficient solution, a method with deallocation
      !may be required for more complex systems. To merger or not is decided elsewhere
      type(Point), intent(inout)      ::  keepPoint, hidePoint
      type(Graph), intent(inout)      ::  graphIn
      real(kind=real64), dimension(2), intent(in)    ::  newPosition

      integer                         ::  conIdx !connection index

      keepPoint%position = newPosition !set new position

      call togglePointActive(graphIn%PointList, hidePoint%id) !set hide point inactive
      graphIn%PointList%noofPointsActive = graphIn%PointList%noofPointsActive - 1 !one less active point

      do conIdx = 1, hidePoint%noofConnex !for all connections to hidPoint
         associate (currentLine => graphIn%LineList%Lines(hidePoint%connex(conIdx))) !alias for clarity
            call AddLine(graphIn%LineList, & !add new lines between keep point and the outer points connected to hide point
                         graphIn%PointList, &
                         keepPoint%id, &
                         getLineOtherEnd(graphIn, currentLine%id, hidePoint%id))
            if (currentLine%active .eqv. .true.) then !if the connection is active
               call toggleLineActive(graphIn%LineList, hidePoint%connex(conIdx))!set it innactive
               graphIn%LineList%noofLinesActive = graphIn%LineList%noofLinesActive - 1 !one less active line
            end if
         end associate

      end do

   end subroutine mergePoints

end module Lib_Graphs
