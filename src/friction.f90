! MIT License
!
! Copyright (c) 2026 Jason Christopherson
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
module friction
    !! Provides a collection of routines for modeling frictional behaviors
    !! of contacting bodies.
    use friction_core
    use friction_coulomb
    use friction_lugre
    use friction_maxwell
    use friction_gmsm
    use friction_stribeck
    use friction_modified_stribeck
    implicit none
    private
    public :: friction_model
    public :: friction_evaluation
    public :: friction_logical_query
    public :: friction_state_model
    public :: friction_model_to_array
    public :: friction_model_from_array
    public :: friction_integer_query
    public :: regression_statistics
    public :: coulomb_model
    public :: lugre_model
    public :: maxwell_model
    public :: generalized_maxwell_slip_model
    public :: stribeck_model
    public :: modified_stribeck_model

end module