module pipes_module
    ! BSD 3-Clause License

    ! Copyright (c) 2014-2023, Jacob Williams
    
    ! Redistribution and use in source and binary forms, with or without
    ! modification, are permitted provided that the following conditions are met:
    
    ! 1. Redistributions of source code must retain the above copyright notice, this
    !    list of conditions and the following disclaimer.
    
    ! 2. Redistributions in binary form must reproduce the above copyright notice,
    !    this list of conditions and the following disclaimer in the documentation
    !    and/or other materials provided with the distribution.
    
    ! 3. Neither the name of the copyright holder nor the names of its
    !    contributors may be used to endorse or promote products derived from
    !    this software without specific prior written permission.
    
    ! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    ! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    ! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    ! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
    ! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    ! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    ! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    ! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    ! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    ! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

    !obtained from this repo: https://github.com/jacobwilliams/popen-fortran/tree/master

    use,intrinsic :: iso_c_binding
    
    implicit none
    
    private
    
    interface
    
        function popen(command, mode) bind(C,name='popen')
        import :: c_char, c_ptr
        character(kind=c_char),dimension(*) :: command
        character(kind=c_char),dimension(*) :: mode
        type(c_ptr) :: popen
        end function popen
    
        function fgets(s, siz, stream) bind(C,name='fgets')
        import :: c_char, c_ptr, c_int
        type (c_ptr) :: fgets
        character(kind=c_char),dimension(*) :: s
        integer(kind=c_int),value :: siz
        type(c_ptr),value :: stream
        end function fgets
    
        function pclose(stream) bind(C,name='pclose')
        import :: c_ptr, c_int
        integer(c_int) :: pclose
        type(c_ptr),value :: stream
        end function pclose
    
    end interface
    
    public :: c2f_string, get_command_as_string
    
    contains
    
    !**********************************************
    ! convert a C string to a Fortran string
    !**********************************************
    function c2f_string(c) result(f)
    
        implicit none
    
        character(len=*),intent(in) :: c
        character(len=:),allocatable :: f
    
        integer :: i
    
        i = index(c,c_null_char)
    
        if (i<=0) then
            f = c
        else if (i==1) then
            f = ''
        else if (i>1) then
            f = c(1:i-1)
        end if
    
    end function c2f_string
    
    !**********************************************
    ! return the result of the command as a string
    !**********************************************
    function get_command_as_string(command) result(str)
    
        implicit none
    
        character(len=*),intent(in) :: command
        character(len=:),allocatable :: str
    
        integer,parameter :: buffer_length = 1000
    
        type(c_ptr) :: h
        integer(c_int) :: istat
        character(kind=c_char,len=buffer_length) :: line
    
        str = ''
        h = c_null_ptr
        h = popen(command//c_null_char,'r'//c_null_char)
    
        if (c_associated(h)) then
            do while (c_associated(fgets(line,buffer_length,h)))
                str = str//c2f_string(line)
            end do
            istat = pclose(h)
        end if
    
    end function get_command_as_string
    
    end module pipes_module
    